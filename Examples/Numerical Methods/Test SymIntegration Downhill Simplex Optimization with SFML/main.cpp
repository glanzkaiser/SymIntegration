
#include "symintegrationc++.h"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <sstream>

const int WIDTH = 1000;
const int HEIGHT = 700;

struct Vector3 {
    float x, y, z;
    Vector3 operator+(const Vector3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vector3 operator-(const Vector3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vector3 operator*(float s) const { return {x * s, y * s, z * s}; }
    Vector3 operator/(float s) const { return {x / s, y / s, z / s}; }
};

struct Vertex {
    Vector3 pos;
    float val;
    bool operator<(const Vertex& other) const { return val < other.val; }
};

enum ObjectiveType {
    RASTRIGIN = 0,
    BOOTH = 1,
    ACKLEY = 2,
    ROSENBROCK = 3,
    COUNTDS = 4
};

ObjectiveType currentType = RASTRIGIN;
float cameraZoom = 50.0f; 

struct LeaderboardEntry {
    ObjectiveType functionType;
    float x, y;
    float absoluteScore;
};
std::vector<LeaderboardEntry> leaderboard;

float getTrueValue(float x, float y, ObjectiveType type) {
    if (type == RASTRIGIN) {
        return 20.0f + (x*x - 10.0f * cos(2.0f * 3.14159265f * x)) 
                     + (y*y - 10.0f * cos(2.0f * 3.14159265f * y));
    } 
    else if (type == BOOTH) {
        float term1 = x + 2.0f * y - 7.0f;
        float term2 = 2.0f * x + y - 5.0f;
        return (term1 * term1) + (term2 * term2);
    } 
    else if (type == ACKLEY) {
        float sum1 = x*x + y*y;
        float sum2 = cos(2.0f * 3.14159265f * x) + cos(2.0f * 3.14159265f * y);
        return -20.0f * exp(-0.2f * sqrt(0.5f * sum1)) - exp(0.5f * sum2) + 20.0f + 2.71828182f;
    }
    else { // ROSENBROCK
        float term1 = 1.0f - x;
        float term2 = y - (x * x);
        return (term1 * term1) + 100.0f * (term2 * term2);
    }
}

float objectiveFunction(float x, float y) {
    float z = getTrueValue(x, y, currentType);
    if (currentType == RASTRIGIN)   return z * 0.2f; 
    else if (currentType == BOOTH)  return z * 0.02f; 
    else if (currentType == ACKLEY) return z * 1.2f;
    else                            return z * 0.02f; // ROSENBROCK
}

void getFunctionMeta(float& minX, float& maxX, float& minY, float& maxY, float& gridStep, Vector3& optimum) {
    if (currentType == RASTRIGIN) {
        minX = -4.0f; maxX = 4.0f; minY = -4.0f; maxY = 4.0f; gridStep = 0.25f;
        optimum = {0.0f, 0.0f, objectiveFunction(0.0f, 0.0f)};
    } else if (currentType == BOOTH) {
        minX = -6.0f; maxX = 8.0f; minY = -4.0f; maxY = 8.0f; gridStep = 0.5f;
        optimum = {1.0f, 3.0f, objectiveFunction(1.0f, 3.0f)};
    } else if (currentType == ACKLEY) {
        minX = -4.0f; maxX = 4.0f; minY = -4.0f; maxY = 4.0f; gridStep = 0.25f;
        optimum = {0.0f, 0.0f, objectiveFunction(0.0f, 0.0f)};
    } else { // ROSENBROCK
        minX = -2.0f; maxX = 2.0f; minY = -1.0f; maxY = 3.0f; gridStep = 0.2f;
        optimum = {1.0f, 1.0f, objectiveFunction(1.0f, 1.0f)};
    }
}

void submitToLeaderboard(const std::vector<Vertex>& simplex) {
    if (simplex.empty()) return;
    float trueScore = getTrueValue(simplex[0].pos.x, simplex[0].pos.y, currentType);
    LeaderboardEntry entry = { currentType, simplex[0].pos.x, simplex[0].pos.y, trueScore };
    leaderboard.push_back(entry);
    
    std::sort(leaderboard.begin(), leaderboard.end(), [](const LeaderboardEntry& a, const LeaderboardEntry& b) {
        return a.absoluteScore < b.absoluteScore;
    });
    
    if (leaderboard.size() > 5) {
        leaderboard.resize(5);
    }
}

void resetSimplex(std::vector<Vertex>& simplex, int& iterations, int& r, int& e, int& c, int& s, std::string& currentOp) {
    if (!simplex.empty()) {
        std::sort(simplex.begin(), simplex.end());
        submitToLeaderboard(simplex);
    }

    iterations = 0; r = 0; e = 0; c = 0; s = 0;
    currentOp = "Context Reset";

    float minX, maxX, minY, maxY, step; Vector3 opt;
    getFunctionMeta(minX, maxX, minY, maxY, step, opt);

    float spanX = (maxX - minX) * 0.4f;
    float spanY = (maxY - minY) * 0.4f;
    float baseX = minX + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / spanX));
    float baseY = minY + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / spanY));

    if (std::abs(baseX - opt.x) < 0.5f) baseX += 1.0f;
    if (std::abs(baseY - opt.y) < 0.5f) baseY += 1.0f;

    float size = (maxX - minX) * 0.15f;
    simplex = {
        {{baseX,        baseY,        0}, 0},
        {{baseX + size, baseY + size, 0}, 0},
        {{baseX + size*0.5f, baseY - size*0.5f, 0}, 0},
        {{baseX + size*1.2f, baseY + size*0.3f, 0}, 0} 
    };

    for (auto& vertex : simplex) {
        vertex.pos.z = objectiveFunction(vertex.pos.x, vertex.pos.y);
        vertex.val = vertex.pos.z; 
    }

    cameraZoom = (currentType == BOOTH) ? 35.0f : 50.0f;
}

sf::Vector2f project(Vector3 p, float angleX, float angleY) {
    float radX = angleX * 3.14159265f / 180.0f;
    float radY = angleY * 3.14159265f / 180.0f;

    float x1 = p.x * cos(radY) - p.z * sin(radY);
    float z1 = p.x * sin(radY) + p.z * cos(radY);
    float y2 = p.y * cos(radX) - z1 * sin(radX);
    float z2 = p.y * sin(radX) + z1 * cos(radX);

    float distance = 15.0f;
    float scale = cameraZoom / (1.0f + z2 / distance); 

    return sf::Vector2f((WIDTH - 220) / 2.0f + x1 * scale, HEIGHT / 2.0f - y2 * scale);
}

int main() {
    srand(static_cast<unsigned>(time(nullptr)));
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Nelder-Mead 3D Surface Visualizer");
    window.setFramerateLimit(60);

    sf::Font font;
    if (!font.loadFromFile("oswald.ttf")) {
        std::cout << "Warning: Could not load oswald.ttf from project folder.\n";
    }

    float sliderX = 800.0f, sliderY = 120.0f, sliderWidth = 160.0f, handleX = sliderX + 80.0f; 
    bool isDragging = false;
    float cameraAngleX = 25.0f, cameraAngleY = -50.0f;
    bool autoResetEnabled = true;

    int iterationCount = 0;
    int reflectionCount = 0, expansionCount = 0, contractionCount = 0, shrinkCount = 0;
    std::string currentOperation = "Initialization";

    std::vector<Vertex> simplex;
    resetSimplex(simplex, iterationCount, reflectionCount, expansionCount, contractionCount, shrinkCount, currentOperation);

    sf::Clock timer;
    bool rKeyWasPressed = false;
    bool tabKeyWasPressed = false;

    while (window.isOpen()) {
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            
            if (event.type == sf::Event::MouseWheelScrolled) {
                if (event.mouseWheelScroll.delta > 0) {
                    cameraZoom += 4.0f;
                } else {
                    cameraZoom -= 4.0f;
                }
                if (cameraZoom < 5.0f)   cameraZoom = 5.0f;
                if (cameraZoom > 250.0f) cameraZoom = 250.0f;
            }
            
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                if (mousePos.x >= handleX - 8 && mousePos.x <= handleX + 8 &&
                    mousePos.y >= sliderY - 12 && mousePos.y <= sliderY + 12) {
                    isDragging = true;
                }
            }
            if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
                isDragging = false;
            }
        }

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Tab)) {
            if (!tabKeyWasPressed) {
                currentType = static_cast<ObjectiveType>((static_cast<int>(currentType) + 1) % static_cast<int>(COUNTDS));
                resetSimplex(simplex, iterationCount, reflectionCount, expansionCount, contractionCount, shrinkCount, currentOperation);
                tabKeyWasPressed = true;
            }
        } else {
            tabKeyWasPressed = false;
        }

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::R)) {
            if (!rKeyWasPressed) {
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::LShift) || sf::Keyboard::isKeyPressed(sf::Keyboard::RShift)) {
                    resetSimplex(simplex, iterationCount, reflectionCount, expansionCount, contractionCount, shrinkCount, currentOperation);
                } else {
                    autoResetEnabled = !autoResetEnabled;
                }
                rKeyWasPressed = true;
            }
        } else {
            rKeyWasPressed = false;
        }

        if (isDragging) {
            handleX = static_cast<float>(mousePos.x);
            if (handleX < sliderX) handleX = sliderX;
            if (handleX > sliderX + sliderWidth) handleX = sliderX + sliderWidth;
        }

        float sliderPercentage = (handleX - sliderX) / sliderWidth;
        float updateInterval = 1.05f - (sliderPercentage * 1.04f); 

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))  cameraAngleY -= 1.5f;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) cameraAngleY += 1.5f;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))    cameraAngleX += 1.5f;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))  cameraAngleX -= 1.5f;

        if (timer.getElapsedTime().asSeconds() > updateInterval) {
            timer.restart();
            iterationCount++;

            std::sort(simplex.begin(), simplex.end());

            float maxDistance = 0.0f;
            for (size_t i = 1; i < simplex.size(); ++i) {
                float dist = std::sqrt(std::pow(simplex[0].pos.x - simplex[i].pos.x, 2) + 
                                       std::pow(simplex[0].pos.y - simplex[i].pos.y, 2));
                if (dist > maxDistance) maxDistance = dist;
            }

            float convergenceLimit = (currentType == BOOTH) ? 0.015f : 0.005f;

            if (autoResetEnabled && maxDistance < convergenceLimit) {
                resetSimplex(simplex, iterationCount, reflectionCount, expansionCount, contractionCount, shrinkCount, currentOperation);
                currentOperation = "Auto-Reset Triggered";
            } else {
                Vector3 centroid = (simplex[0].pos + simplex[1].pos + simplex[2].pos) / 3.0f;
                Vector3 worst = simplex[3].pos;
                Vector3 reflected = centroid + (centroid - worst) * 1.0f;
                float reflectedVal = objectiveFunction(reflected.x, reflected.y);
                reflected.z = reflectedVal;

                if (reflectedVal >= simplex[0].val && reflectedVal < simplex[2].val) {
                    simplex[3] = {reflected, reflectedVal};
                    currentOperation = "Reflection"; reflectionCount++;
                } 
                else if (reflectedVal < simplex[0].val) {
                    Vector3 expanded = centroid + (reflected - centroid) * 2.0f;
                    float expandedVal = objectiveFunction(expanded.x, expanded.y);
                    expanded.z = expandedVal;

                    if (expandedVal < reflectedVal) {
                        simplex[3] = {expanded, expandedVal};
                        currentOperation = "Expansion"; expansionCount++;
                    } else {
                        simplex[3] = {reflected, reflectedVal};
                        currentOperation = "Reflection"; reflectionCount++;
                    }
                } 
                else {
                    Vector3 contracted = centroid + (worst - centroid) * 0.5f;
                    float contractedVal = objectiveFunction(contracted.x, contracted.y);
                    contracted.z = contractedVal;

                    if (contractedVal < simplex[3].val) {
                        simplex[3] = {contracted, contractedVal};
                        currentOperation = "Contraction"; contractionCount++;
                    } 
                    else {
                        currentOperation = "Shrink Step"; shrinkCount++;
                        for (size_t i = 1; i < simplex.size(); ++i) {
                            simplex[i].pos = simplex[0].pos + (simplex[i].pos - simplex[0].pos) * 0.5f;
                            simplex[i].pos.z = objectiveFunction(simplex[i].pos.x, simplex[i].pos.y);
                            simplex[i].val = simplex[i].pos.z;
                        }
                    }
                }
            }
        }

        window.clear(sf::Color(255, 255, 255));

        float minX, maxX, minY, maxY, gridStep; Vector3 optimum;
        getFunctionMeta(minX, maxX, minY, maxY, gridStep, optimum);

        sf::VertexArray grid(sf::Lines, 0);
        for (float x = minX; x <= maxX; x += gridStep) {
            for (float y = minY; y <= maxY; y += gridStep) {
                Vector3 p1 = {x, y, objectiveFunction(x, y)};
                Vector3 p2 = {x + gridStep, y, objectiveFunction(x + gridStep, y)};
                Vector3 p3 = {x, y + gridStep, objectiveFunction(x, y + gridStep)};

                if (x + gridStep <= maxX) {
                    grid.append(sf::Vertex(project(p1, cameraAngleX, cameraAngleY), sf::Color(60, 65, 85)));
                    grid.append(sf::Vertex(project(p2, cameraAngleX, cameraAngleY), sf::Color(60, 65, 85)));
                }
                if (y + gridStep <= maxY) {
                    grid.append(sf::Vertex(project(p1, cameraAngleX, cameraAngleY), sf::Color(60, 65, 85)));
                    grid.append(sf::Vertex(project(p3, cameraAngleX, cameraAngleY), sf::Color(60, 65, 85)));
                }
            }
        }
        window.draw(grid);

        sf::Vector2f targetScreenPos = project(optimum, cameraAngleX, cameraAngleY);
        sf::CircleShape targetMarker(8.0f);
        targetMarker.setFillColor(sf::Color::Magenta);
        targetMarker.setOutlineThickness(2.0f);
        targetMarker.setOutlineColor(sf::Color::White);
        targetMarker.setPosition(targetScreenPos.x - 8.0f, targetScreenPos.y - 8.0f);
        window.draw(targetMarker);

        sf::VertexArray simplexEdges(sf::Lines, 0);
        for (size_t i = 0; i < simplex.size(); ++i) {
            for (size_t j = i + 1; j < simplex.size(); ++j) {
                simplexEdges.append(sf::Vertex(project(simplex[i].pos, cameraAngleX, cameraAngleY), sf::Color::Cyan));
                simplexEdges.append(sf::Vertex(project(simplex[j].pos, cameraAngleX, cameraAngleY), sf::Color::Cyan));
            }
        }
        window.draw(simplexEdges);

        for (const auto& v : simplex) {
            sf::CircleShape dot(5.0f);
            dot.setFillColor(sf::Color::Yellow);
            sf::Vector2f screenPos = project(v.pos, cameraAngleX, cameraAngleY);
            dot.setPosition(screenPos.x - 5.0f, screenPos.y - 5.0f);
            window.draw(dot);
        }

        sf::RectangleShape uiBg(sf::Vector2f(220.0f, HEIGHT));
        uiBg.setPosition(WIDTH - 220.0f, 0);
        uiBg.setFillColor(sf::Color(30, 30, 40));
        window.draw(uiBg);

        sf::Text txtFunctionTitle, txtIteration, txtOp, txtSpeed, txtTarget, txtAutoReset, txtCountersHeader, txtRef, txtCon, txtInstructions;
        std::vector<sf::Text*> allTexts = { &txtFunctionTitle, &txtIteration, &txtOp, &txtSpeed, &txtTarget, &txtAutoReset, &txtCountersHeader, &txtRef, &txtCon, &txtInstructions };
        for (auto* t : allTexts) {
            t->setFont(font); t->setCharacterSize(13); t->setFillColor(sf::Color::White);
        }

        std::string fName = "UNKNOWN";
        if(currentType == RASTRIGIN) fName = "RASTRIGIN";
        else if(currentType == BOOTH) fName = "BOOTH";
        else if(currentType == ACKLEY) fName = "ACKLEY";
        else if(currentType == ROSENBROCK) fName = "ROSENBROCK";

        txtFunctionTitle.setString("Landscape: " + fName);
        txtFunctionTitle.setPosition(WIDTH - 200, 15);
        txtFunctionTitle.setStyle(sf::Text::Bold);
        txtFunctionTitle.setFillColor(sf::Color::Yellow);
        window.draw(txtFunctionTitle);

        txtIteration.setString("Iteration: " + std::to_string(iterationCount));
        txtIteration.setPosition(WIDTH - 200, 40);
        window.draw(txtIteration);

        txtOp.setString("Last Action:\n" + currentOperation);
        txtOp.setPosition(WIDTH - 200, 60);
        txtOp.setFillColor(sf::Color::Cyan);
        window.draw(txtOp);

        txtSpeed.setString("Algorithm Speed:");
        txtSpeed.setPosition(WIDTH - 200, 100);
        window.draw(txtSpeed);

        sf::RectangleShape sliderTrack(sf::Vector2f(sliderWidth, 4.0f));
        sliderTrack.setPosition(sliderX, sliderY);
        sliderTrack.setFillColor(sf::Color(100, 100, 120));
        window.draw(sliderTrack);

        sf::CircleShape sliderHandle(8.0f);
        sliderHandle.setFillColor(sf::Color::Green);
        sliderHandle.setPosition(handleX - 8.0f, sliderY - 6.0f);
        window.draw(sliderHandle);

        std::string fMinCoord = "(0.0, 0.0)";
        if(currentType == BOOTH) fMinCoord = "(1.0, 3.0)";
        else if(currentType == ROSENBROCK) fMinCoord = "(1.0, 1.0)";

        txtTarget.setString("Target Min: " + fMinCoord);
        txtTarget.setPosition(WIDTH - 200, 140);
        txtTarget.setFillColor(sf::Color::Magenta);
        window.draw(txtTarget);

        txtAutoReset.setString("Auto-Reset: " + std::string(autoResetEnabled ? "ENABLED" : "DISABLED"));
        txtAutoReset.setPosition(WIDTH - 200, 165);
        txtAutoReset.setFillColor(autoResetEnabled ? sf::Color::Green : sf::Color::Red);
        window.draw(txtAutoReset);

        txtCountersHeader.setString("Operation Totals:");
        txtCountersHeader.setPosition(WIDTH - 200, 195);
        txtCountersHeader.setStyle(sf::Text::Underlined);
        window.draw(txtCountersHeader);

        txtRef.setString("Ref: " + std::to_string(reflectionCount) + " | Exp: " + std::to_string(expansionCount));
        txtRef.setPosition(WIDTH - 200, 215); window.draw(txtRef);

        txtCon.setString("Con: " + std::to_string(contractionCount) + " | Shr: " + std::to_string(shrinkCount));
        txtCon.setPosition(WIDTH - 200, 235); window.draw(txtCon);

        sf::Text txtBoardHeader;
        txtBoardHeader.setFont(font); txtBoardHeader.setCharacterSize(14);
        txtBoardHeader.setFillColor(sf::Color::Yellow);
        txtBoardHeader.setStyle(sf::Text::Underlined | sf::Text::Bold);
        txtBoardHeader.setString("Top 5 Abs Min Points:");
        txtBoardHeader.setPosition(WIDTH - 200, 270);
        window.draw(txtBoardHeader);

        for (size_t i = 0; i < leaderboard.size(); ++i) {
            sf::Text txtEntry;
            txtEntry.setFont(font); txtEntry.setCharacterSize(11);
            txtEntry.setFillColor(sf::Color(220, 220, 240));

            std::string entryLabel = "UNK";
            if (leaderboard[i].functionType == RASTRIGIN) entryLabel = "RAST";
            else if (leaderboard[i].functionType == BOOTH) entryLabel = "BOOT";
            else if (leaderboard[i].functionType == ACKLEY) entryLabel = "ACKL";
            else if (leaderboard[i].functionType == ROSENBROCK) entryLabel = "ROSB";

            std::stringstream ss;
            ss << (i + 1) << ". [" << entryLabel << "] " << std::scientific << std::setprecision(2) << leaderboard[i].absoluteScore;
            txtEntry.setString(ss.str());
            txtEntry.setPosition(WIDTH - 200, 295 + (i * 20));
            window.draw(txtEntry);
        }

        txtInstructions.setString("Controls:\n[Scroll Wheel] Zoom View\n[Tab] Cycle Landscapes\n[R] Toggle Auto-Reset\n[Shift+R] Force Reset\n[Arrows] Rotate View");
        txtInstructions.setPosition(WIDTH - 200, HEIGHT - 105);
        txtInstructions.setFillColor(sf::Color(160, 160, 170));
        txtInstructions.setCharacterSize(11);
        window.draw(txtInstructions);

        window.display();
    }

    return 0;
}
