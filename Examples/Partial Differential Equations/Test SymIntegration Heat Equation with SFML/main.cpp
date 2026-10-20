#include <SFML/Graphics.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "symintegrationc++.h"

// Window Configuration
const int WINDOW_WIDTH = 1200;
const int WINDOW_HEIGHT = 900;
float baseScale = 160.0f; 
float zoomFactor = 1.0f; 
const float HEIGHT_SCALE = 1.5f; 

float maxInitialTemp = 2.0f;
const int STEPS_PER_FRAME = 10;

// Legend Bar Layout Constants
const int LEGEND_X = 1100;
const int LEGEND_WIDTH = 25;
const int LEGEND_HEIGHT = 500;
const int LEGEND_Y = (WINDOW_HEIGHT - LEGEND_HEIGHT) / 2;

// States
Shape currentShape = CRESCENT;
bool isNeumann = false;
bool useImplicit = false;
bool isPaused = false;

float angleX = 110.0f * 3.14159265f / 180.0f;
float angleZ = -55.0f * 3.14159265f / 180.0f;

// 3D to 2D Projection
sf::Vector2f project(float x, float y, float z) {
    float currentScale = baseScale * zoomFactor;
    float x1 = x * std::cos(angleZ) - y * std::sin(angleZ);
    float y1 = x * std::sin(angleZ) + y * std::cos(angleZ);
    float x2 = x1;
    float y2 = y1 * std::cos(angleX) - (z * HEIGHT_SCALE) * std::sin(angleX);
    return sf::Vector2f((WINDOW_WIDTH / 2.0f - 100.0f) + x2 * currentScale, WINDOW_HEIGHT / 2.0f + y2 * currentScale);
}

// True Jet Colormap
sf::Color getTemperatureColor(float temp, float maxTemp) {
    if (maxTemp <= 0.001f) maxTemp = 0.001f;
    float v = temp / maxTemp;
    if (v < 0.0f) v = 0.0f;
    if (v > 1.0f) v = 1.0f;
    float r = 0.0f, g = 0.0f, b = 0.0f;
    if (v < 0.25f) { float t = v / 0.25f; r = 0.0f; g = t; b = 1.0f; }
    else if (v < 0.5f) { float t = (v - 0.25f) / 0.25f; r = 0.0f; g = 1.0f; b = 1.0f - t; }
    else if (v < 0.75f) { float t = (v - 0.5f) / 0.25f; r = t; g = 1.0f; b = 0.0f; }
    else { float t = (v - 0.75f) / 0.25f; r = 1.0f; g = 1.0f - t; b = 0.0f; }
    return sf::Color(static_cast<sf::Uint8>(r * 255), static_cast<sf::Uint8>(g * 255), static_cast<sf::Uint8>(b * 255));
}

int main() {
    // SFML 2.5.1 Anti-aliasing context setup
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Advanced 3D Heat Simulation Sandbox", sf::Style::Default, settings);
    window.setFramerateLimit(60);

    sf::Font font;
    bool fontLoaded = font.loadFromFile("oswald.ttf");
    if (!fontLoaded) {
        std::cerr << "Warning: Could not load oswald.ttf from folder. Text labels/HUD will not render.\n";
    }

    Grid U(N, std::vector<float>(N, 0.0f));
    Grid nextU(N, std::vector<float>(N, 0.0f));
    resetGrid(U, currentShape, maxInitialTemp);

    sf::VertexArray mesh(sf::Triangles, (N - 1) * (N - 1) * 6);
    sf::VertexArray legendMesh(sf::Quads);
    
    int slices = 100;
    float sliceHeight = static_cast<float>(LEGEND_HEIGHT) / slices;
    for (int s = 0; s < slices; ++s) {
        float norm1 = static_cast<float>(s) / slices;
        float norm2 = static_cast<float>(s + 1) / slices;
        sf::Color c1 = getTemperatureColor(norm1, 1.0f);
        sf::Color c2 = getTemperatureColor(norm2, 1.0f);
        float y1 = LEGEND_Y + LEGEND_HEIGHT - (s * sliceHeight);
        float y2 = LEGEND_Y + LEGEND_HEIGHT - ((s + 1) * sliceHeight);
        legendMesh.append(sf::Vertex(sf::Vector2f(static_cast<float>(LEGEND_X), y1), c1));
        legendMesh.append(sf::Vertex(sf::Vector2f(static_cast<float>(LEGEND_X + LEGEND_WIDTH), y1), c1));
        legendMesh.append(sf::Vertex(sf::Vector2f(static_cast<float>(LEGEND_X + LEGEND_WIDTH), y2), c2));
        legendMesh.append(sf::Vertex(sf::Vector2f(static_cast<float>(LEGEND_X), y2), c2));
    }

    const int numLabels = 5;
    std::vector<sf::Text> labels(numLabels);
    if (fontLoaded) {
        for (int i = 0; i < numLabels; ++i) {
            labels[i].setFont(font);
            labels[i].setCharacterSize(16);
            labels[i].setFillColor(sf::Color::White); // SFML 2.5.1 Compliant
        }
    }

    sf::Text hudText;
    if (fontLoaded) {
        hudText.setFont(font);
        hudText.setCharacterSize(16);
        hudText.setFillColor(sf::Color(230, 230, 255));
        hudText.setPosition(20.0f, 20.0f);
    }

    bool isMousePressed = false;
    sf::Vector2i lastMousePos;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();

            // Drag to rotate camera
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                isMousePressed = true;
                lastMousePos = sf::Mouse::getPosition(window);
            }
            if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
                isMousePressed = false;
            }
            if (event.type == sf::Event::MouseMoved && isMousePressed) {
                sf::Vector2i currentMousePos = sf::Mouse::getPosition(window);
                angleZ += (currentMousePos.x - lastMousePos.x) * 0.005f;
                angleX += (currentMousePos.y - lastMousePos.y) * 0.005f;
                lastMousePos = currentMousePos;
            }
            
            // SFML 2.5.1 Scroll Wheel syntax (MouseWheelScrolled)
            if (event.type == sf::Event::MouseWheelScrolled) {
                if (event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel) {
                    if (event.mouseWheelScroll.delta > 0) zoomFactor *= 1.1f;
                    else if (event.mouseWheelScroll.delta < 0) zoomFactor /= 1.1f;
                    zoomFactor = std::max(0.2f, std::min(zoomFactor, 5.0f));
                }
            }
            
            // Toggles
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Space) isPaused = !isPaused;
                if (event.key.code == sf::Keyboard::S) {
                    currentShape = static_cast<Shape>((currentShape + 1) % 3);
                    resetGrid(U, currentShape, maxInitialTemp);
                }
                if (event.key.code == sf::Keyboard::B) isNeumann = !isNeumann;
                if (event.key.code == sf::Keyboard::I) useImplicit = !useImplicit;
                if (event.key.code == sf::Keyboard::Up) maxInitialTemp += 0.5f;
                if (event.key.code == sf::Keyboard::Down) maxInitialTemp = std::max(0.5f, maxInitialTemp - 0.5f);
            }
        }

        // Right-click heat paint interaction
        if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
            sf::Vector2i mousePos = sf::Mouse::getPosition(window);
            float minDist = 1e9;
            int bestI = -1, bestJ = -1;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    sf::Vector2f sPos = project(-L + i * dx, -L + j * dx, U[i][j]);
                    float dist = std::hypot(sPos.x - mousePos.x, sPos.y - mousePos.y);
                    if (dist < minDist) { minDist = dist; bestI = i; bestJ = j; }
                }
            }
            if (minDist < 40.0f && bestI != -1 && bestJ != -1) {
                for (int di = -2; di <= 2; ++di) {
                    for (int dj = -2; dj <= 2; ++dj) {
                        int ni = bestI + di, nj = bestJ + dj;
                        if (ni >= 0 && ni < N && nj >= 0 && nj < N) U[ni][nj] = maxInitialTemp;
                    }
                }
            }
        }

        // Math Execution Pipeline Step
        if (!isPaused) {
            for (int step = 0; step < STEPS_PER_FRAME; ++step) {
                if (!useImplicit) {
                    stepExplicit(U, nextU);
                } else {
                    stepImplicit(U, nextU);
                }
                applyBoundaryConditions(U, isNeumann);
            }
        }

        // UI Label Logic
        if (fontLoaded) {
            for (int i = 0; i < numLabels; ++i) {
                float val = (maxInitialTemp / (numLabels - 1)) * i;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(1) << val;
                labels[i].setString(ss.str());
                float labelY = LEGEND_Y + LEGEND_HEIGHT - (static_cast<float>(i) / (numLabels - 1) * LEGEND_HEIGHT);
                sf::FloatRect textRect = labels[i].getLocalBounds();
                labels[i].setOrigin(0, textRect.height / 2.0f);
                labels[i].setPosition(static_cast<float>(LEGEND_X + LEGEND_WIDTH + 10), labelY);
            }

            std::stringstream hudStream;
            hudStream << "=========================================\n"
                      << "SYSTEM STATUS (HUD): " << (isPaused ? "[PAUSED]" : "[RUNNING]") << "\n"
                      << " [Space] Simulation State: " << (isPaused ? "Paused" : "Simulating") << "\n"
                      << " [S] Shape: " << (currentShape == CRESCENT ? "Crescent Moon" : currentShape == CYLINDER ? "Cylinder" : "Box") << "\n"
                      << " [B] Boundary: " << (isNeumann ? "Insulated (Neumann)" : "Freezing Cold (Dirichlet)") << "\n"
                      << " [I] Method: " << (useImplicit ? "Implicit (Crank-Nicolson)" : "Explicit (Forward Euler)") << "\n"
                      << " [Up/Down] Max Temp Limit: " << std::fixed << std::setprecision(1) << maxInitialTemp << "\n"
                      << "-----------------------------------------\n"
                      << " -> Scroll Mouse Wheel to Zoom In / Out\n"
                      << " -> Drag Left Mouse Click to Rotate Camera\n"
                      << " -> Drag Right Mouse Click to Inject Heat\n"
                      << "=========================================\n";
            hudText.setString(hudStream.str());
        }

        // Update 3D Vertex Map Arrays
        size_t vIdx = 0;
        for (int i = 0; i < N - 1; ++i) {
            for (int j = 0; j < N - 1; ++j) {
                sf::Vector2f p00 = project(-L + i*dx, -L + j*dx, U[i][j]);
                sf::Vector2f p10 = project(-L + (i+1)*dx, -L + j*dx, U[i+1][j]);
                sf::Vector2f p01 = project(-L + i*dx, -L + (j+1)*dx, U[i][j+1]);
                sf::Vector2f p11 = project(-L + (i+1)*dx, -L + (j+1)*dx, U[i+1][j+1]);

                mesh[vIdx++] = sf::Vertex(p00, getTemperatureColor(U[i][j], maxInitialTemp));
                mesh[vIdx++] = sf::Vertex(p10, getTemperatureColor(U[i+1][j], maxInitialTemp));
                mesh[vIdx++] = sf::Vertex(p01, getTemperatureColor(U[i][j+1], maxInitialTemp));

                mesh[vIdx++] = sf::Vertex(p10, getTemperatureColor(U[i+1][j], maxInitialTemp));
                mesh[vIdx++] = sf::Vertex(p11, getTemperatureColor(U[i+1][j+1], maxInitialTemp));
                mesh[vIdx++] = sf::Vertex(p01, getTemperatureColor(U[i][j+1], maxInitialTemp));
            }
        }

        // Draw Cycle
        window.clear(sf::Color(15, 15, 18));
        window.draw(mesh);
        window.draw(legendMesh);
        if (fontLoaded) {
            for (int i = 0; i < numLabels; ++i) window.draw(labels[i]);
            window.draw(hudText);
        }
        window.display();
    }

    return 0;
}
