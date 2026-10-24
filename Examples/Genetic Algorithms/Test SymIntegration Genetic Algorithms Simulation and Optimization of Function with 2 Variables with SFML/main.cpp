#include <SFML/Graphics.hpp>
#include "symintegrationc++.h"
#include <ctime>
#include <cstdlib>

const int WIDTH = 1400; 
const int HEIGHT = 900;
const int GENERATIONS_PER_FRAME = 1; 

FunctionType currentFunc = ROSENBROCK;
float cameraZoom = 600.0f; 

std::vector<float> fitnessHistory;
const size_t MAX_HISTORY_POINTS = 150;

// Graphical projection engine utilizing math constructs
sf::Vector2f project(Point3D p, float angleX, float angleY) {
    float radY = angleY * PI_VAL / 180.0f;
    float x1 = p.x * std::cos(radY) - p.z * std::sin(radY);
    float z1 = p.x * std::sin(radY) + p.z * std::cos(radY);

    float radX = angleX * PI_VAL / 180.0f;
    float y2 = p.y * std::cos(radX) - z1 * std::sin(radX);
    float z2 = p.y * std::sin(radX) + z1 * std::cos(radX);

    float distance = 4.0f; 
    float sz = distance + z2;
    if (sz <= 0.1f) sz = 0.1f;

    Domain d = getDomain(currentFunc);
    float normX = ((x1 - (d.xMin + d.xMax) / 2.0f) / (d.xMax - d.xMin)) * 2.0f;

    float screenX = (normX / sz) * cameraZoom + (WIDTH * 0.42f);
    float screenY = -(y2 / sz) * cameraZoom + (HEIGHT * 0.42f); 

    return sf::Vector2f(screenX, screenY);
}

// GUI UI Components
class Slider {
public:
    sf::RectangleShape track;
    sf::CircleShape handle;
    float minValue, maxValue, currentValue;
    bool isDragging;
    std::string label;

    Slider(float x, float y, float w, float minV, float maxV, float startV, std::string lbl) {
        minValue = minV; maxValue = maxV; currentValue = startV;
        isDragging = false; label = lbl;

        track.setSize(sf::Vector2f(w, 4.0f));
        track.setFillColor(sf::Color(180, 180, 190));
        track.setPosition(x, y);

        handle.setRadius(8.0f);
        handle.setFillColor(sf::Color(0, 120, 240));
        handle.setOrigin(8.0f, 8.0f);
        
        float fraction = (currentValue - minValue) / (maxValue - minValue);
        handle.setPosition(x + fraction * w, y + 2.0f);
    }

    void handleEvent(const sf::Event& event, const sf::RenderWindow& window) {
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
            if (handle.getGlobalBounds().contains(mousePos.x, mousePos.y)) isDragging = true;
        }
        if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) isDragging = false;
        
        if (isDragging && event.type == sf::Event::MouseMoved) {
            float trackLeft = track.getPosition().x;
            float trackWidth = track.getSize().x;
            float newX = std::max(trackLeft, std::min(static_cast<float>(mousePos.x), trackLeft + trackWidth));
            handle.setPosition(newX, handle.getPosition().y);
            float fraction = (newX - trackLeft) / trackWidth;
            currentValue = minValue + fraction * (maxValue - minValue);
        }
    }

    void draw(sf::RenderWindow& window, sf::Font& font) {
        window.draw(track);
        window.draw(handle);
        sf::Text text; text.setFont(font); text.setCharacterSize(14); text.setFillColor(sf::Color::Black);
        std::stringstream ss; ss << label << ": " << std::fixed << std::setprecision(2) << currentValue;
        text.setString(ss.str()); text.setPosition(track.getPosition().x, track.getPosition().y - 22.0f);
        window.draw(text);
    }
};

class Button {
public:
    sf::RectangleShape shape;
    sf::Text text;

    Button(float x, float y, float w, float h, std::string label, sf::Font& font) {
        shape.setSize(sf::Vector2f(w, h));
        shape.setFillColor(sf::Color(230, 235, 245));
        shape.setOutlineThickness(1.5f);
        shape.setOutlineColor(sf::Color(160, 170, 190));
        shape.setPosition(x, y);

        text.setFont(font);
        text.setString(label);
        text.setCharacterSize(14);
        text.setFillColor(sf::Color::Black);
        
        sf::FloatRect textBounds = text.getLocalBounds();
        text.setOrigin(textBounds.left + textBounds.width / 2.0f, textBounds.top + textBounds.height / 2.0f);
        text.setPosition(x + w / 2.0f, y + h / 2.0f);
    }

    bool isClicked(const sf::Event& event, const sf::RenderWindow& window) {
        if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
            sf::Vector2i mousePos = sf::Mouse::getPosition(window);
            return shape.getGlobalBounds().contains(mousePos.x, mousePos.y);
        }
        return false;
    }

    void updateHover(const sf::RenderWindow& window) {
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        if (shape.getGlobalBounds().contains(mousePos.x, mousePos.y)) {
            shape.setFillColor(sf::Color(210, 220, 235));
        } else {
            shape.setFillColor(sf::Color(230, 235, 245));
        }
    }

    void draw(sf::RenderWindow& window) {
        window.draw(shape);
        window.draw(text);
    }
};

void rebuildMesh(std::vector<std::vector<Point3D>>& mesh, int size) {
    Domain d = getDomain(currentFunc);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            float x = d.xMin + (d.xMax - d.xMin) * i / (size - 1);
            float y = d.yMin + (d.yMax - d.yMin) * j / (size - 1);
            mesh[i][j] = { x, getZHeight(x, y, currentFunc), y };
        }
    }
}

void clearGraphHistory() {
    fitnessHistory.clear();
}

int main() {
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Modular 3D GA Optimization Playground");
    window.setFramerateLimit(60);

    sf::Font font;
    if (!font.loadFromFile("oswald.ttf")) {
        if (!font.loadFromFile("arial.ttf")) {
            if (!font.loadFromFile("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")) {
                std::cout << "Warning: Missing visual UI typography components!" << std::endl;
            }
        }
    }

    sf::Text leftHUD; leftHUD.setFont(font); leftHUD.setCharacterSize(16);
    leftHUD.setFillColor(sf::Color::Black); leftHUD.setPosition(20.f, 20.f);

    sf::Text rightHUD; rightHUD.setFont(font); rightHUD.setCharacterSize(15);
    rightHUD.setFillColor(sf::Color(20, 80, 40)); rightHUD.setPosition(WIDTH - 360.f, 20.f);

    Slider mutationSlider(WIDTH - 360.f, 540.f, 250.f, 0.0f, 1.0f, 0.15f, "Mutation Rate");
    Slider crossoverSlider(WIDTH - 360.f, 610.f, 250.f, 0.0f, 1.0f, 0.5f, "Crossover Bias (Alpha)");
    Button resetButton(WIDTH - 360.f, 660.f, 250.f, 35.f, "Reset Population (R)", font);

    GeneticAlgorithm ga;
    ga.reset(currentFunc);
    clearGraphHistory();

    float cameraAngleX = 35.0f;
    float cameraAngleY = -50.0f;
    
    bool isRotatingCamera = false;
    sf::Vector2i lastMousePos;

    const int gridSize = 40;
    std::vector<std::vector<Point3D>> mesh(gridSize, std::vector<Point3D>(gridSize));
    rebuildMesh(mesh, gridSize);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Tab) {
                currentFunc = static_cast<FunctionType>((currentFunc + 1) % 4);
                rebuildMesh(mesh, gridSize);
                ga.reset(currentFunc); 
                clearGraphHistory();
            }

            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::R) {
                ga.reset(currentFunc);
                clearGraphHistory();
            }

            if (event.type == sf::Event::MouseWheelScrolled) {
                if (event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel) {
                    cameraZoom += event.mouseWheelScroll.delta * 40.0f;
                    if (cameraZoom < 100.0f) cameraZoom = 100.0f;
                    if (cameraZoom > 2500.0f) cameraZoom = 2500.0f;
                }
            }

            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                sf::Vector2i mousePos = sf::Mouse::getPosition(window);
                if (!mutationSlider.handle.getGlobalBounds().contains(mousePos.x, mousePos.y) &&
                    !crossoverSlider.handle.getGlobalBounds().contains(mousePos.x, mousePos.y) &&
                    !resetButton.shape.getGlobalBounds().contains(mousePos.x, mousePos.y) &&
                    mousePos.x < (WIDTH - 380) && mousePos.y < 650) { 
                    isRotatingCamera = true;
                    lastMousePos = mousePos;
                }
            }

            if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
                isRotatingCamera = false;
            }

            if (event.type == sf::Event::MouseMoved && isRotatingCamera) {
                sf::Vector2i currentMousePos = sf::Mouse::getPosition(window);
                float deltaX = static_cast<float>(currentMousePos.x - lastMousePos.x);
                float deltaY = static_cast<float>(currentMousePos.y - lastMousePos.y);
                cameraAngleY += deltaX * 0.4f;
                cameraAngleX += deltaY * 0.4f;
                lastMousePos = currentMousePos;
            }

            if (resetButton.isClicked(event, window)) {
                ga.reset(currentFunc);
                clearGraphHistory();
            }

            mutationSlider.handleEvent(event, window);
            crossoverSlider.handleEvent(event, window);
        }

        resetButton.updateHover(window);

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))  cameraAngleY -= 1.5f;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) cameraAngleY += 1.5f;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))    cameraAngleX += 1.5f;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))  cameraAngleX -= 1.5f;

        if (!mutationSlider.isDragging && !crossoverSlider.isDragging) {
            for (int g = 0; g < GENERATIONS_PER_FRAME; ++g) {
                ga.evolve(currentFunc, mutationSlider.currentValue, crossoverSlider.currentValue);
            }
        }

        float bestX = 0.0f, bestY = 0.0f, bestCost = 0.0f;
        if (!ga.population.empty()) {
            bestX = ga.population[0].x;
            bestY = ga.population[0].y;
            bestCost = evaluateObjective(bestX, bestY, currentFunc);
        }

        if (fitnessHistory.empty() || ga.generation > static_cast<int>(fitnessHistory.size())) {
            fitnessHistory.push_back(bestCost);
            if (fitnessHistory.size() > MAX_HISTORY_POINTS) {
                fitnessHistory.erase(fitnessHistory.begin());
            }
        }

        window.clear(sf::Color::White); 

        // Draw Surface Mesh Map Landscape
        for (int i = 0; i < gridSize; ++i) {
            sf::VertexArray line(sf::LinesStrip, gridSize);
            for (int j = 0; j < gridSize; ++j) {
                line[j].position = project(mesh[i][j], cameraAngleX, cameraAngleY);
                line[j].color = sf::Color(0, 90, 200, 45);
            }
            window.draw(line);
        }
        for (int j = 0; j < gridSize; ++j) {
            sf::VertexArray line(sf::LinesStrip, gridSize);
            for (int i = 0; i < gridSize; ++i) {
                line[i].position = project(mesh[i][j], cameraAngleX, cameraAngleY);
                line[i].color = sf::Color(0, 90, 200, 45);
            }
            window.draw(line);
        }

        // --- RENDER 3D AXIS MARKERS (X, Y, Z) AT ORIGIN ---
        Domain d = getDomain(currentFunc);
        float axisLengthX = (d.xMax - d.xMin) * 0.3f;
        float axisLengthY = (d.yMax - d.yMin) * 0.3f;
        float axisLengthZ = 1.5f; 

        Point3D origin = { 0.0f, 0.0f, 0.0f };
        sf::Vector2f screenOrigin = project(origin, cameraAngleX, cameraAngleY);

        sf::Vector2f screenX = project({ axisLengthX, 0.0f, 0.0f }, cameraAngleX, cameraAngleY);
        sf::Vector2f screenY = project({ 0.0f, 0.0f, axisLengthY }, cameraAngleX, cameraAngleY); 
        sf::Vector2f screenZ = project({ 0.0f, axisLengthZ, 0.0f }, cameraAngleX, cameraAngleY);

        sf::VertexArray axisLines(sf::Lines, 6);
        // X Axis -> Red
        axisLines[0] = sf::Vertex(screenOrigin, sf::Color(230, 50, 50, 200));
        axisLines[1] = sf::Vertex(screenX,      sf::Color(230, 50, 50, 200));
        // Y Axis -> Green
        axisLines[2] = sf::Vertex(screenOrigin, sf::Color(40, 180, 40, 200));
        axisLines[3] = sf::Vertex(screenY,      sf::Color(40, 180, 40, 200));
        // Z Axis -> Blue
        axisLines[4] = sf::Vertex(screenOrigin, sf::Color(40, 100, 230, 200));
        axisLines[5] = sf::Vertex(screenZ,      sf::Color(40, 100, 230, 200));
        
        window.draw(axisLines);

        // Axis Text Labels
        sf::Text labelX("X", font, 14); labelX.setFillColor(sf::Color(200, 30, 30));   labelX.setPosition(screenX);
        sf::Text labelY("Y", font, 14); labelY.setFillColor(sf::Color(30, 150, 30));   labelY.setPosition(screenY);
        sf::Text labelZ("Z", font, 14); labelZ.setFillColor(sf::Color(30, 80, 200));   labelZ.setPosition(screenZ);
        window.draw(labelX);
        window.draw(labelY);
        window.draw(labelZ);

        // Draw Agents Optimization Nodes
        for (const auto& ind : ga.population) {
            Point3D agentPos = { ind.x, getZHeight(ind.x, ind.y, currentFunc), ind.y };
            sf::Vector2f screenPos = project(agentPos, cameraAngleX, cameraAngleY);

            sf::CircleShape agentNode(3.5f);
            if (ind.origin == "Mutated")      agentNode.setFillColor(sf::Color(255, 120, 0)); 
            else if (ind.origin == "Elite")   agentNode.setFillColor(sf::Color(0, 180, 60)); 
            else                              agentNode.setFillColor(sf::Color(230, 30, 30));  
                
            agentNode.setOrigin(3.5f, 3.5f);
            agentNode.setPosition(screenPos);
            window.draw(agentNode);
        }

        // --- RENDER REAL-TIME FITNESS DECAY GRAPH OVERLAY ---
        float gX = 50.0f, gY = 660.0f, gW = 950.0f, gH = 180.0f;
        
        sf::RectangleShape graphBG(sf::Vector2f(gW, gH));
        graphBG.setPosition(gX, gY);
        graphBG.setFillColor(sf::Color(245, 247, 250, 230));
        graphBG.setOutlineThickness(1.5f);
        graphBG.setOutlineColor(sf::Color(200, 210, 225));
        window.draw(graphBG);

        float maxValInHistory = 1.0f;
        for (float val : fitnessHistory) {
            if (val > maxValInHistory) maxValInHistory = val;
        }

        if (fitnessHistory.size() > 1) {
            sf::VertexArray graphLine(sf::LinesStrip, fitnessHistory.size());
            for (size_t i = 0; i < fitnessHistory.size(); ++i) {
                float progressFraction = static_cast<float>(i) / (MAX_HISTORY_POINTS - 1);
                float plotX = gX + progressFraction * gW;
                
                float heightFraction = fitnessHistory[i] / maxValInHistory;
                float plotY = (gY + gH) - (heightFraction * (gH - 20.0f)) - 10.0f;

                graphLine[i].position = sf::Vector2f(plotX, plotY);
                graphLine[i].color = sf::Color(0, 120, 245);
            }
            window.draw(graphLine);
        }

        sf::Text txtGraphTitle; txtGraphTitle.setFont(font); txtGraphTitle.setCharacterSize(13);
        txtGraphTitle.setFillColor(sf::Color(60, 70, 90)); txtGraphTitle.setPosition(gX + 15.0f, gY + 8.0f);
        txtGraphTitle.setString("REAL-TIME FITNESS DECAY PROFILE (OBJECTIVE MINIMIZATION)");
        window.draw(txtGraphTitle);

        sf::Text txtGraphScale; txtGraphScale.setFont(font); txtGraphScale.setCharacterSize(12);
        txtGraphScale.setFillColor(sf::Color(120, 130, 150)); txtGraphScale.setPosition(gX + 15.0f, gY + gH - 22.0f);
        std::stringstream ssScale; ssScale << "Peak Vertical Scale bounds: [0.00 -> " << std::scientific << std::setprecision(2) << maxValInHistory << "]";
        txtGraphScale.setString(ssScale.str());
        window.draw(txtGraphScale);

        // Display Dashboard Text Metrics
        std::stringstream ssLeft;
        ssLeft << "=== SIMULATION ENGINE ===\n"
               << "Active Landscape: \n>> " << getFunctionName(currentFunc) << " <<\n"
               << "Generation: " << ga.generation << "\n"
               << "Population Size: " << POP_SIZE << "\n\n"
               << "CURRENT BEST CHROMOSOME:\n"
               << "X: " << std::fixed << std::setprecision(4) << bestX << "\n"
               << "Y: " << std::fixed << std::setprecision(4) << bestY << "\n"
               << "Objective Cost: " << bestCost << "\n\n"
               << "[Tab]         : Switch Landscape\n"
               << "[R]           : Reset Population\n"
               << "[Left-Drag]   : Rotate Camera View\n"
               << "[Mouse Scroll]: Zoom In/Out View";
        leftHUD.setString(ssLeft.str());

        std::stringstream ssRight;
        ssRight << "=== POPULATION TELEMETRY ===\n"
                << "Lineage Operations (This Frame):\n"
                << "  Elite Selected   : " << ga.countElite << "\n"
                << "  Crossovers Built  : " << ga.countCrossover << "\n"
                << "  Mutations Applied : " << ga.countMutation << "\n\n"
                << "=== TOP 3 HISTORICAL MINIMA ===\n";
        for (size_t i = 0; i < ga.historicalMinima.size(); ++i) {
            float hX = ga.historicalMinima[i].x;
            float hY = ga.historicalMinima[i].y;
            ssRight << "Rank #" << (i + 1) << ": f(" << std::setprecision(2) << hX << ", " << hY << ") = " 
                    << std::scientific << std::setprecision(2) << evaluateObjective(hX, hY, currentFunc) << "\n";
        }
        ssRight << "\n\n\n\n\n=== INTERACTIVE PARAMETERS ===";
        rightHUD.setString(ssRight.str());

        window.draw(leftHUD);
        window.draw(rightHUD);

        mutationSlider.draw(window, font);
        crossoverSlider.draw(window, font);
        resetButton.draw(window);

        window.display();
    }

    return 0;
}
