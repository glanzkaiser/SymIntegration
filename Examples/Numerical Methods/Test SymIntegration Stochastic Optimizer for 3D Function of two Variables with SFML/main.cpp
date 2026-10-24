#include <SFML/Graphics.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "symintegrationc++.h"

// --- Configuration Constants ---
const int WINDOW_WIDTH = 1200;
const int WINDOW_HEIGHT = 800;


enum OptimizerType { DIRECT_MONTE_CARLO, STOCHASTIC_TUNNELING, PARALLEL_TEMPERING, NUM_OPTIMIZERS };

// --- 3D Grid Configuration ---
const int GRID_SIZE = 40; 
const float RANGE_X = 2.5f;
const float RANGE_Y = 2.5f;

// ============================================================================
// Local Objective Function Definition (Rendered explicitly as z = f(x, y))
// ============================================================================
double objective_function_3d(double x, double y, double s, double t) 
{
    double r2 = x * x + y * y;
    double term1 = std::exp(-r2 / (2.0 * s * s));
    double term2 = std::exp(-r2 / (2.0 * t * t));
    double penalty = 0.1 * r2;
    return term1 - term2 + penalty;
}

// Wrapper matching standard solver hook interface vectors
std::function<double(const std::vector<double>&)> make_solver_lambda(double s, double t) 
{
    return [s, t](const std::vector<double>& v) {
        return objective_function_3d(v[0], v[1], s, t);
    };
}

// ============================================================================
// Projection Math & Geometry Utilities
// ============================================================================
//struct Point3D { float x, y, z; };

sf::Vector2f project_isometric(Point3D p, float angleX, float angleZ, float scaleX, float scaleY, float scaleZ, sf::Vector2f center) 
{
    // Rotate around Z axis (Yaw)
    float x1 = p.x * std::cos(angleZ) - p.y * std::sin(angleZ);
    float y1 = p.x * std::sin(angleZ) + p.y * std::cos(angleZ);
    
    // Rotate around X axis (Pitch)
    float y2 = y1 * std::cos(angleX) - p.z * std::sin(angleX);
    float z2 = y1 * std::sin(angleX) + p.z * std::cos(angleX);

    // Apply viewport scale constants and relative center translation
    return sf::Vector2f(center.x + x1 * scaleX, center.y - y2 * scaleY + z2 * scaleZ);
}

// Simple Interactive Slider Class Component
struct Slider {
    sf::RectangleShape track;
    sf::CircleShape knob;
    float minValue, maxValue, currentValue;
    bool isDragging = false;
    sf::Text label;

    Slider(float x, float y, float w, float minV, float maxV, float startV, std::string title, const sf::Font& font) {
        minValue = minV; maxValue = maxV; currentValue = startV;
        track.setSize(sf::Vector2f(w, 4));
        track.setFillColor(sf::Color(180, 180, 180));
        track.setPosition(x, y);

        knob.setRadius(8);
        knob.setFillColor(sf::Color(40, 120, 220));
        knob.setOrigin(8, 8);
        updateKnobPosition();

        label.setFont(font);
        label.setCharacterSize(14);
        label.setFillColor(sf::Color(40, 40, 40));
        label.setPosition(x, y - 22);
    }

    void updateKnobPosition() {
        float ratio = (currentValue - minValue) / (maxValue - minValue);
        knob.setPosition(track.getPosition().x + ratio * track.getSize().x, track.getPosition().y + 2);
    }

    bool handleEvent(const sf::Event& event, const sf::RenderWindow& window) {
        if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
            sf::Vector2f mousePos(event.mouseButton.x, event.mouseButton.y);
            if (knob.getGlobalBounds().contains(mousePos)) {
                isDragging = true;
                return true; 
            }
        }
        if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
            if (isDragging) {
                isDragging = false;
                return true;
            }
        }
        if (isDragging && event.type == sf::Event::MouseMoved) {
            float mouseX = event.mouseMove.x;
            float trackX = track.getPosition().x;
            float trackW = track.getSize().x;
            float ratio = (mouseX - trackX) / trackW;
            ratio = std::max(0.0f, std::min(1.0f, ratio));
            currentValue = minValue + ratio * (maxValue - minValue);
            updateKnobPosition();
            return true;
        }
        return false;
    }

    void updateLabel(std::string prefix) {
        std::stringstream ss;
        ss << prefix << ": " << std::fixed << std::setprecision(3) << currentValue;
        label.setString(ss.str());
    }

    void draw(sf::RenderWindow& window) {
        window.draw(track);
        window.draw(knob);
        window.draw(label);
    }
};

// ============================================================================
// MAIN SYSTEM ORCHESTRATOR
// ============================================================================
int main() 
{
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "3D Stochastic Optimization", sf::Style::Default, settings);
    window.setFramerateLimit(60);

    sf::Font font;
    if (!font.loadFromFile("oswald.ttf")) {
        if (!font.loadFromFile("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")) {
            std::cerr << "Warning: System Font resources failed to load." << std::endl;
        }
    }

    // Initialize Viewport Position Mappings
    sf::Vector2f meshCenter(450.0f, 380.0f);
    float angleX = 55.0f * PI / 180.0f; 
    float angleZ = 45.0f * PI / 180.0f;
    
    // Global scaling zoom multiplier factor
    float zoomFactor = 1.0f;
    float baseScaleX = 140.0f, baseScaleY = 140.0f, baseScaleZ = 250.0f;

    OptimizerType currentOptimizer = DIRECT_MONTE_CARLO;
    Slider sliderS(920.0f, 80.0f, 240.0f, 0.2f, 1.5f, 0.75f, "Parameter s", font);
    Slider sliderT(920.0f, 140.0f, 240.0f, 0.5f, 2.5f, 1.0f, "Parameter t", font);

    std::vector<double> bestFoundX = {1.0, 1.0};
    double bestFoundEnergy = 0.0;
    bool triggerRerun = true;

    sf::Text headerText("Optimization Wireframe Engine", font, 22);
    headerText.setFillColor(sf::Color(20, 20, 20));
    headerText.setPosition(20, 15);

    sf::Text controlsText("[Tab]: Toggle Optimizer  |  [Left Click + Drag]: Rotate View  |  [Scroll Wheel]: Zoom", font, 13);
    controlsText.setFillColor(sf::Color(100, 100, 100));
    controlsText.setPosition(20, 48);

    sf::Text statsText("", font, 14);
    statsText.setFillColor(sf::Color(30, 30, 30));
    statsText.setPosition(920.0f, 200.0f);

    bool isRotatingView = false;
    sf::Vector2i lastMousePos;

    // ============================================================================
    // APPLICATION ENGINE RENDERING INTERACTION LOOP
    // ============================================================================
    while (window.isOpen()) 
    {
        sf::Event event;
        while (window.pollEvent(event)) 
        {
            if (event.type == sf::Event::Closed) window.close();

            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Tab) 
            {
                currentOptimizer = static_cast<OptimizerType>((currentOptimizer + 1) % NUM_OPTIMIZERS);
                triggerRerun = true;
            }

            // Handle UI Interaction Event Blocks (Sliders check consumption first)
            float oldS = sliderS.currentValue;
            float oldT = sliderT.currentValue;
            bool sliderActive = false;
            
            if (sliderS.handleEvent(event, window)) sliderActive = true;
            if (sliderT.handleEvent(event, window)) sliderActive = true;
            
            if (sliderS.currentValue != oldS || sliderT.currentValue != oldT) {
                triggerRerun = true; 
            }

            // --- 1. Zoom In/Out via Mouse Scroll Wheel Events ---
            if (event.type == sf::Event::MouseWheelScrolled) {
                if (event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel) {
                    zoomFactor += event.mouseWheelScroll.delta * 0.08f;
                    zoomFactor = std::max(0.2f, std::min(4.0f, zoomFactor));
                }
            }

            // --- 2. Camera Rotation via Left Click + Drag Interaction ---
            if (!sliderActive) {
                if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                    isRotatingView = true;
                    lastMousePos = sf::Mouse::getPosition(window);
                }
                if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
                    isRotatingView = false;
                }
                if (isRotatingView && event.type == sf::Event::MouseMoved) {
                    sf::Vector2i currentMousePos = sf::Mouse::getPosition(window);
                    angleZ += (currentMousePos.x - lastMousePos.x) * 0.005f;
                    angleX += (currentMousePos.y - lastMousePos.y) * 0.005f;
                    lastMousePos = currentMousePos;
                }
            }
        }

        float currentScaleX = baseScaleX * zoomFactor;
        float currentScaleY = baseScaleY * zoomFactor;
        float currentScaleZ = baseScaleZ * zoomFactor;

        // --- Execute Active Solver Runtime Context ---
        if (triggerRerun) 
        {
            double sVal = sliderS.currentValue;
            double tVal = sliderT.currentValue;
            auto target_hook = make_solver_lambda(sVal, tVal);

            if (currentOptimizer == DIRECT_MONTE_CARLO) {
                DMCMetrics res = run_direct_monte_carlo(15000, {-2.5, -2.5}, {2.5, 2.5}, target_hook);
                bestFoundX = res.best_x; bestFoundEnergy = res.best_energy;
            } 
            else if (currentOptimizer == STOCHASTIC_TUNNELING) {
                STUNMetrics res = run_stochastic_tunneling(8000, 0.06, 2.0, 0.04, {1.2, 1.2}, target_hook);
                bestFoundX = res.best_x; bestFoundEnergy = res.best_energy;
            } 
            else if (currentOptimizer == PARALLEL_TEMPERING) {
                PTMetrics res = run_parallel_tempering(300, 10, {1.2, 1.2}, {0.04, 0.2, 1.0, 4.0}, target_hook);
                bestFoundX = res.best_x; bestFoundEnergy = res.best_energy;
            }
            triggerRerun = false;
        }

        window.clear(sf::Color::White);

        // --- Render UI Text Elements ---
        sliderS.updateLabel("Gaussian Decay s");
        sliderT.updateLabel("Envelope Drag t");
        sliderS.draw(window);
        sliderT.draw(window);

        std::string optName[] = {"Direct Monte Carlo", "Stochastic Tunneling", "Parallel Tempering"};
        std::stringstream ss;
        ss << "Active Strategy: >> " << optName[currentOptimizer] << ""
           << "\nMin Energy Found:" << std::fixed << std::setprecision(8) << bestFoundEnergy << ""
           << "\nDiscovered Vector:" 
           << "\nX: " << std::setprecision(4) << bestFoundX[0] << "\nY: " << bestFoundX[1] << ""
           << "\nCurrent Zoom: " << std::setprecision(1) << zoomFactor * 100.0f << "%";
        statsText.setString(ss.str());
        window.draw(statsText);
        window.draw(headerText);
        window.draw(controlsText);

        // ============================================================================
        // 1. Draw 3D Spatial Axis Reference Markers At System Core Center
        // ============================================================================
        sf::VertexArray axisLines(sf::Lines, 6);
        float axisLength = 3.0f;
        sf::Vector2f originProjected = project_isometric({0,0,0}, angleX, angleZ, currentScaleX, currentScaleY, currentScaleZ, meshCenter);

        axisLines[0] = sf::Vertex(originProjected, sf::Color(230, 80, 80));
        axisLines[1] = sf::Vertex(project_isometric({axisLength, 0, 0}, angleX, angleZ, currentScaleX, currentScaleY, currentScaleZ, meshCenter), sf::Color(230, 80, 80));
        axisLines[2] = sf::Vertex(originProjected, sf::Color(80, 200, 80));
        axisLines[3] = sf::Vertex(project_isometric({0, axisLength, 0}, angleX, angleZ, currentScaleX, currentScaleY, currentScaleZ, meshCenter), sf::Color(80, 200, 80));
        axisLines[4] = sf::Vertex(originProjected, sf::Color(80, 80, 230));
        axisLines[5] = sf::Vertex(project_isometric({0, 0, axisLength}, angleX, angleZ, currentScaleX, currentScaleY, currentScaleZ, meshCenter), sf::Color(80, 80, 230));
        window.draw(axisLines);

        // ============================================================================
        // 2. Draw Mathematical Objective Function 3D Mesh Wireframe Context
        // ============================================================================
        std::vector<std::vector<sf::Vector2f>> projectedGrid(GRID_SIZE + 1, std::vector<sf::Vector2f>(GRID_SIZE + 1));
        double sVal = sliderS.currentValue;
        double tVal = sliderT.currentValue;

        for (int i = 0; i <= GRID_SIZE; ++i) {
            float gx = -RANGE_X + (2.0f * RANGE_X * i) / GRID_SIZE;
            for (int j = 0; j <= GRID_SIZE; ++j) {
                float gy = -RANGE_Y + (2.0f * RANGE_Y * j) / GRID_SIZE;
                float gz = static_cast<float>(objective_function_3d(gx, gy, sVal, tVal));
                projectedGrid[i][j] = project_isometric({gx, gy, gz}, angleX, angleZ, currentScaleX, currentScaleY, currentScaleZ, meshCenter);
            }
        }

        sf::VertexArray wireframe(sf::Lines, 0);
        sf::Color meshColor(190, 195, 205);
        for (int i = 0; i <= GRID_SIZE; ++i) {
            for (int j = 0; j <= GRID_SIZE; ++j) {
                if (i < GRID_SIZE) {
                    wireframe.append(sf::Vertex(projectedGrid[i][j], meshColor));
                    wireframe.append(sf::Vertex(projectedGrid[i+1][j], meshColor));
                }
                if (j < GRID_SIZE) {
                    wireframe.append(sf::Vertex(projectedGrid[i][j], meshColor));
                    wireframe.append(sf::Vertex(projectedGrid[i][j+1], meshColor));
                }
            }
        }
        window.draw(wireframe);

        // Overlay Best Solved Node Position Indicator Node inside the 3D Space
        if (bestFoundX.size() >= 2) {
            float targetZ = static_cast<float>(objective_function_3d(bestFoundX[0], bestFoundX[1], sVal, tVal));
            sf::Vector2f node2D = project_isometric({(float)bestFoundX[0], (float)bestFoundX[1], targetZ}, angleX, angleZ, currentScaleX, currentScaleY, currentScaleZ, meshCenter);
            sf::CircleShape trackerDot(6);
            trackerDot.setFillColor(sf::Color(255, 60, 60));
            trackerDot.setOutlineColor(sf::Color::White);
            trackerDot.setOutlineThickness(1.5f);
            trackerDot.setOrigin(6, 6);
            trackerDot.setPosition(node2D);
            window.draw(trackerDot);
        }

        // ============================================================================
        // 3. Draw Radial Cross-Section Mathematical Plot Engine Block Below
        // ============================================================================
        sf::Vector2f plotOrigin(120.0f, 720.0f);
        float plotWidth = 700.0f;
        float plotHeight = 160.0f;

        sf::RectangleShape plotBorder(sf::Vector2f(plotWidth, plotHeight));
        plotBorder.setPosition(plotOrigin.x, plotOrigin.y - plotHeight);
        plotBorder.setFillColor(sf::Color(250, 250, 250));
        plotBorder.setOutlineColor(sf::Color(210, 210, 210));
        plotBorder.setOutlineThickness(1.0f);
        window.draw(plotBorder);

        sf::VertexArray analyticalCurve(sf::LineStrip, 0);
        int curveSamples = 200;
        float maxRadius = 3.0f;
        
        for (int k = 0; k < curveSamples; ++k) {
            float r = (maxRadius * k) / (curveSamples - 1);
            float e = static_cast<float>(objective_function_3d(r, 0.0, sVal, tVal));
            
            float cx = plotOrigin.x + (r / maxRadius) * plotWidth;
            float cy = plotOrigin.y - (plotHeight / 2.0f) - (e * plotHeight * 1.5f);
            
            if (cy >= plotOrigin.y - plotHeight && cy <= plotOrigin.y) {
                analyticalCurve.append(sf::Vertex(sf::Vector2f(cx, cy), sf::Color(80, 80, 90)));
            }
        }
        window.draw(analyticalCurve);

        float bestR = std::sqrt(bestFoundX[0]*bestFoundX[0] + bestFoundX[1]*bestFoundX[1]);
        float pinX = plotOrigin.x + (bestR / maxRadius) * plotWidth;
        float pinY = plotOrigin.y - (plotHeight / 2.0f) - (static_cast<float>(bestFoundEnergy) * plotHeight * 1.5f);
        
        if (pinX >= plotOrigin.x && pinX <= plotOrigin.x + plotWidth) {
            sf::CircleShape graphPin(4);
            graphPin.setFillColor(sf::Color(240, 50, 50));
            graphPin.setOrigin(4, 4);
            graphPin.setPosition(pinX, pinY);
            window.draw(graphPin);
        }

        sf::Text graphLabel("Radial Cross-Section Landscape Curve vs Optimizer Discovered Position Pin", font, 11);
        graphLabel.setFillColor(sf::Color(90, 95, 105));
        graphLabel.setPosition(plotOrigin.x, plotOrigin.y - plotHeight - 18);
        window.draw(graphLabel);

        window.display();
    }

    return 0;
}