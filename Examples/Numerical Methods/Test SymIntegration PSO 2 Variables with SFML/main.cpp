#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "symintegrationc++.h"

// Window and Layout Constraints
const int WORKSPACE_WIDTH = 1100;
const int HEIGHT = 900;
const int UI_PANEL_WIDTH = 300;
const int TOTAL_WIDTH = WORKSPACE_WIDTH + UI_PANEL_WIDTH; // 1100 Total Width

const int NUM_PARTICLES = 40;
const float X_MIN = -5.0f;
const float X_MAX = 5.0f;
const float Y_MIN = -5.0f;
const float Y_MAX = 5.0f;

// Optimization Target Function
float multimodalLandscape(float x, float y) {
    return std::sin(x) * std::cos(y) + 0.1f * (x * x + y * y);
}

// Jet Color Map implementation for the landscape and colorbar matching MATLAB style
sf::Color jetColorMap(float val, float minV, float maxV) {
    float norm = (val - minV) / (maxV - minV);
    norm = std::max(0.0f, std::min(1.0f, norm));
    
    // Jet mapping math curves
    float r = std::max(0.0f, std::min(1.0f, 1.5f - std::abs(4.0f * norm - 3.0f)));
    float g = std::max(0.0f, std::min(1.0f, 1.5f - std::abs(4.0f * norm - 2.0f)));
    float b = std::max(0.0f, std::min(1.0f, 1.5f - std::abs(4.0f * norm - 1.0f)));
    
    return sf::Color(static_cast<sf::Uint8>(r * 255), 
                     static_cast<sf::Uint8>(g * 255), 
                     static_cast<sf::Uint8>(b * 255));
}

// Simple internal structure for handling native SFML basic slider interactions
struct Slider {
    sf::RectangleShape track;
    sf::CircleShape handle;
    float minVal;
    float maxVal;
    float currentVal;
    bool isDragging;
    std::string labelName;

    void init(float x, float y, float w, float minV, float maxV, float startV, std::string label) {
        minVal = minV;
        maxVal = maxV;
        currentVal = startV;
        isDragging = false;
        labelName = label;

        track.setSize(sf::Vector2f(w, 6.0f));
        track.setFillColor(sf::Color(70, 70, 70));
        track.setPosition(x, y);

        handle.setRadius(10.0f);
        handle.setOrigin(10.0f, 10.0f);
        handle.setFillColor(sf::Color(200, 200, 200));
        handle.setOutlineThickness(2.0f);
        handle.setOutlineColor(sf::Color(30, 30, 30));
        
        // Position handle relative to current initial variable mapping 
        float ratio = (currentVal - minVal) / (maxVal - minVal);
        handle.setPosition(x + ratio * w, y + 3.0f);
    }

    void handleEvent(const sf::Event& event, const sf::RenderWindow& window) {
        if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
            sf::Vector2f mousePos = window.mapPixelToCoords(sf::Vector2i(event.mouseButton.x, event.mouseButton.y));
            if (handle.getGlobalBounds().contains(mousePos.x, mousePos.y)) {
                isDragging = true;
            }
        }
        if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left) {
            isDragging = false;
        }
    }

    void update(const sf::RenderWindow& window) {
        if (isDragging) {
            sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
            sf::Vector2f mousePos = window.mapPixelToCoords(pixelPos);
            float trackLeft = track.getPosition().x;
            float trackWidth = track.getSize().x;
            
            float newX = std::max(trackLeft, std::min(trackLeft + trackWidth, mousePos.x));
            handle.setPosition(newX, handle.getPosition().y);
            
            float ratio = (newX - trackLeft) / trackWidth;
            currentVal = minVal + ratio * (maxVal - minVal);
        }
    }

    void draw(sf::RenderWindow& window, sf::Text& textTool) {
        window.draw(track);
        window.draw(handle);

        std::stringstream ss;
        ss << labelName << ": " << std::fixed << std::setprecision(2) << currentVal;
        textTool.setString(ss.str());
        textTool.setPosition(track.getPosition().x, track.getPosition().y - 25.0f);
        window.draw(textTool);
    }
};

int main() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    sf::RenderWindow window(sf::VideoMode(TOTAL_WIDTH, HEIGHT), "PSO 2D System Simulation - SFML 2.5.1", sf::Style::Titlebar | sf::Style::Close);
    window.setFramerateLimit(30);

    // Font Handling Check
    sf::Font font;
    if (!font.loadFromFile("oswald.ttf")) {
        std::cerr << "Error loading arial.ttf from runtime working path." << std::endl;
        return -1;
    }

    sf::Text uiText;
    uiText.setFont(font);
    uiText.setCharacterSize(14);
    uiText.setFillColor(sf::Color::White);

    ObjectiveFunction currentFunc = multimodalLandscape;

    // --- Generate Contour Jet Color Texture Background ---
    sf::Image bgImage;
    bgImage.create(WORKSPACE_WIDTH, HEIGHT);
    
    float minF = 1e9f;
    float maxF = -1e9f;
    
    for (int py = 0; py < HEIGHT; ++py) {
        for (int px = 0; px < WORKSPACE_WIDTH; ++px) {
            float x = X_MIN + (static_cast<float>(px) / WORKSPACE_WIDTH) * (X_MAX - X_MIN);
            float y = Y_MIN + (static_cast<float>(py) / HEIGHT) * (Y_MAX - Y_MIN);
            float f = currentFunc(x, y);
            if (f < minF) minF = f;
            if (f > maxF) maxF = f;
        }
    }

    for (int py = 0; py < HEIGHT; ++py) {
        for (int px = 0; px < WORKSPACE_WIDTH; ++px) {
            float x = X_MIN + (static_cast<float>(px) / WORKSPACE_WIDTH) * (X_MAX - X_MIN);
            float y = Y_MIN + (static_cast<float>(py) / HEIGHT) * (Y_MAX - Y_MIN);
            bgImage.setPixel(px, py, jetColorMap(currentFunc(x, y), minF, maxF));
        }
    }

    sf::Texture bgTexture;
    bgTexture.loadFromImage(bgImage);
    sf::Sprite bgSprite(bgTexture);

    // --- Generate Colorbar Legend Texture ---
    const int BAR_WIDTH = 30;
    const int BAR_HEIGHT = 300;
    sf::Image barImage;
    barImage.create(BAR_WIDTH, BAR_HEIGHT);
    for (int y = 0; y < BAR_HEIGHT; ++y) {
        // Invert y index logic so high values remain vertically near the top layout bounds
        float ratio = 1.0f - (static_cast<float>(y) / BAR_HEIGHT);
        float value = minF + ratio * (maxF - minF);
        sf::Color colorSample = jetColorMap(value, minF, maxF);
        for (int x = 0; x < BAR_WIDTH; ++x) {
            barImage.setPixel(x, y, colorSample);
        }
    }
    sf::Texture barTexture;
    barTexture.loadFromImage(barImage);
    sf::Sprite barSprite(barTexture);
    barSprite.setPosition(WORKSPACE_WIDTH + 40, 40);

    // --- Instantiate Interactive Component Elements ---
    SwarmState swarm;
    initSwarm(swarm, NUM_PARTICLES, X_MIN, X_MAX, Y_MIN, Y_MAX, currentFunc);

    Slider wSlider, c1Slider, c2Slider;
    wSlider.init(WORKSPACE_WIDTH + 40, 410, 200, 0.0f, 1.5f, 0.5f, "w (Inertia)");
    c1Slider.init(WORKSPACE_WIDTH + 40, 470, 200, 0.0f, 4.0f, 1.5f, "c1 (Cognitive)");
    c2Slider.init(WORKSPACE_WIDTH + 40, 530, 200, 0.0f, 4.0f, 1.5f, "c2 (Social)");

    // UI Separator Line 
    sf::Vertex uiDivider[] = {
        sf::Vertex(sf::Vector2f(WORKSPACE_WIDTH, 0), sf::Color(100, 100, 100)),
        sf::Vertex(sf::Vector2f(WORKSPACE_WIDTH, HEIGHT), sf::Color(100, 100, 100))
    };

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::R) {
                initSwarm(swarm, NUM_PARTICLES, X_MIN, X_MAX, Y_MIN, Y_MAX, currentFunc);
            }

            wSlider.handleEvent(event, window);
            c1Slider.handleEvent(event, window);
            c2Slider.handleEvent(event, window);
        }

        // Live calculation loops based on modern UI settings
        wSlider.update(window);
        c1Slider.update(window);
        c2Slider.update(window);

        PSOHyperparams activeParams = { wSlider.currentVal, c1Slider.currentVal, c2Slider.currentVal };
        updateSwarmMath(swarm, activeParams, X_MIN, X_MAX, Y_MIN, Y_MAX, currentFunc);

        // --- Render Target Phase ---
        window.clear(sf::Color(30, 30, 30));
        window.draw(bgSprite);

        // Draw Swarm Entities 
        for (int i = 0; i < NUM_PARTICLES; ++i) {
            float screenX = ((swarm.particles[i].x - X_MIN) / (X_MAX - X_MIN)) * WORKSPACE_WIDTH;
            float screenY = ((swarm.particles[i].y - Y_MIN) / (Y_MAX - Y_MIN)) * HEIGHT;

            sf::CircleShape particleDot(4.0f);
            particleDot.setOrigin(4.0f, 4.0f);
            particleDot.setPosition(screenX, screenY);
            particleDot.setFillColor(sf::Color::White);
            particleDot.setOutlineThickness(1.0f);
            particleDot.setOutlineColor(sf::Color::Black);
            window.draw(particleDot);

            float scale = 40.0f;
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(screenX, screenY), sf::Color::Black),
                sf::Vertex(sf::Vector2f(screenX + swarm.particles[i].vx * scale, screenY + swarm.particles[i].vy * scale), sf::Color::Black)
            };
            window.draw(line, 2, sf::Lines);
        }

        // --- Render UI Workspace Interface Panel ---
        window.draw(uiDivider, 2, sf::Lines);
        window.draw(barSprite);

        // Draw Colorbar Text Labels
        uiText.setCharacterSize(12);
        std::stringstream ssMax, ssMin;
        ssMax << std::fixed << std::setprecision(2) << maxF;
        uiText.setString(ssMax.str());
        uiText.setPosition(WORKSPACE_WIDTH + 80, 40);
        window.draw(uiText);

        ssMin << std::fixed << std::setprecision(2) << minF;
        uiText.setString(ssMin.str());
        uiText.setPosition(WORKSPACE_WIDTH + 80, 40 + BAR_HEIGHT - 15);
        window.draw(uiText);

        // Draw Live Global Best Metrics Output
        uiText.setCharacterSize(15);
        std::stringstream ssMetrics;
        ssMetrics << "--- Swarm Metrics ---"
                  << "Global Best Fitness:\n" << std::scientific << std::setprecision(5) << swarm.globalBestFitness << "\n\n"
                  << "Global Best X:\n" << std::fixed << std::setprecision(4) << swarm.globalBestX << "\n"
                  << "Global Best Y:\n" << std::fixed << std::setprecision(4) << swarm.globalBestY << "\n\n"
                  << "Press [R] to Reset Swarm";
        uiText.setString(ssMetrics.str());
        uiText.setPosition(WORKSPACE_WIDTH + 40, 600);
        window.draw(uiText);

        // Render sliders
        wSlider.draw(window, uiText);
        c1Slider.draw(window, uiText);
        c2Slider.draw(window, uiText);

        window.display();
    }

    return 0;
}
