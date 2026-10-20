#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

// Window Configuration
const int WINDOW_WIDTH = 1200;
const int WINDOW_HEIGHT = 900;
const int N = 80; // Grid size (N x N)
const float L = 2.5f; // Domain boundaries (-L to L)
const float SCALE = 160.0f; // Visual scale factor for the screen
const float HEIGHT_SCALE = 1.5f; // Vertical multiplier for temperature height

// Dynamic simulation scaling parameter
float maxInitialTemp = 2.0f;

// Legend Bar Layout Constants
const int LEGEND_X = 1100; // Placed far right, adjusted slightly for text space
const int LEGEND_WIDTH = 25;
const int LEGEND_HEIGHT = 500;
const int LEGEND_Y = (WINDOW_HEIGHT - LEGEND_HEIGHT) / 2; // Vertically centered

// Simulation physics constants
const float D = 1.0f; // Thermal diffusivity
const float dx = (2.0f * L) / (N - 1);
const float dt = 0.125f * (dx * dx) / D; // CFL explicit limit
const int STEPS_PER_FRAME = 10;

// Interaction States
enum Shape { CRESCENT, CYLINDER, BOX };
Shape currentShape = CRESCENT;
bool isNeumann = false; // False = Dirichlet (Freezing), True = Neumann (Insulated)
bool useImplicit = false; // False = Forward Euler, True = Crank-Nicolson

// Global viewing angles
float angleX = 110.0f * 3.14159265f / 180.0f;
float angleZ = -55.0f * 3.14159265f / 180.0f;

// 3D to 2D Projection Function
sf::Vector2f project(float x, float y, float z) {
    float x1 = x * std::cos(angleZ) - y * std::sin(angleZ);
    float y1 = x * std::sin(angleZ) + y * std::cos(angleZ);
    float x2 = x1;
    float y2 = y1 * std::cos(angleX) - (z * HEIGHT_SCALE) * std::sin(angleX);
    
    // Shift center to the left to fully isolate the legend bar
    return sf::Vector2f((WINDOW_WIDTH / 2.0f - 100.0f) + x2 * SCALE, WINDOW_HEIGHT / 2.0f + y2 * SCALE);
}

// True ParaView "Jet" Colormap
sf::Color getTemperatureColor(float temp, float maxTemp) {
    if (maxTemp <= 0.001f) maxTemp = 0.001f;
    float v = temp / maxTemp;
    if (v < 0.0f) v = 0.0f;
    if (v > 1.0f) v = 1.0f;
    
    float r = 0.0f, g = 0.0f, b = 0.0f;
    if (v < 0.25f) {
        float t = v / 0.25f;
        r = 0.0f; g = t; b = 1.0f;
    } else if (v < 0.5f) {
        float t = (v - 0.25f) / 0.25f;
        r = 0.0f; g = 1.0f; b = 1.0f - t;
    } else if (v < 0.75f) {
        float t = (v - 0.5f) / 0.25f;
        r = t; g = 1.0f; b = 0.0f;
    } else {
        float t = (v - 0.75f) / 0.25f;
        r = 1.0f; g = 1.0f - t; b = 0.0f;
    }
    return sf::Color(static_cast<sf::Uint8>(r * 255), static_cast<sf::Uint8>(g * 255), static_cast<sf::Uint8>(b * 255));
}

// Initialize Shapes
void resetGrid(std::vector<std::vector<float>>& U) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float x = -L + i * dx;
            float y = -L + j * dx;
            if (currentShape == CRESCENT) {
                float c = 2.0f, d = -1.0f, e = 1.0f, f = 0.5f, k = 1.2f, shift = 10.0f;
                float pattern = std::pow(c*c - std::pow(x/e - d, 2) - std::pow(y/f, 2), 2) + k * std::pow(c + d - x/e, 3) - shift;
                U[i][j] = (pattern < 0.0f) ? maxInitialTemp : 0.0f;
            } else if (currentShape == CYLINDER) {
                U[i][j] = (x*x + y*y < 1.0f) ? maxInitialTemp : 0.0f;
            } else if (currentShape == BOX) {
                U[i][j] = (std::abs(x) < 0.8f && std::abs(y) < 0.8f) ? maxInitialTemp : 0.0f;
            }
        }
    }
}

int main() {
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Advanced 3D Heat Simulation Sandbox", sf::Style::Default, settings);
    window.setFramerateLimit(60);

    // Load Font for Labels and HUD
    sf::Font font;
    bool fontLoaded = font.loadFromFile("oswald.ttf");
    if (!fontLoaded) {
        std::cerr << "Warning: Could not load arial.ttf from folder. Text labels/HUD will not render.\n";
    }

    std::vector<std::vector<float>> U(N, std::vector<float>(N, 0.0f));
    std::vector<std::vector<float>> nextU(N, std::vector<float>(N, 0.0f));
    resetGrid(U);

    sf::VertexArray mesh(sf::Triangles, (N - 1) * (N - 1) * 6);

    // Set up static geometry for the vertical color scale bar mesh
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

    // Initialize structures for Legend Text Labels (5 labels distributed from 0.0 to maxInitialTemp)
    const int numLabels = 5;
    std::vector<sf::Text> labels(numLabels);
    if (fontLoaded) {
        for (int i = 0; i < numLabels; ++i) {
            labels[i].setFont(font);
            labels[i].setCharacterSize(16);
            labels[i].setFillColor(sf::Color::White);
        }
    }

    // Configure HUD Text Overlay (Top-Left)
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
            if (event.type == sf::Event::Closed)
                window.close();

            // Camera control via mouse drag
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

            // Controls Key Toggles
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::S) {
                    currentShape = static_cast<Shape>((currentShape + 1) % 3);
                    resetGrid(U);
                }
                if (event.key.code == sf::Keyboard::B) {
                    isNeumann = !isNeumann;
                }
                if (event.key.code == sf::Keyboard::I) {
                    useImplicit = !useImplicit;
                }
                // Dynamically adjust maxInitialTemp with Arrow Keys on the fly
                if (event.key.code == sf::Keyboard::Up) {
                    maxInitialTemp += 0.5f;
                }
                if (event.key.code == sf::Keyboard::Down) {
                    maxInitialTemp = std::max(0.5f, maxInitialTemp - 0.5f);
                }
            }
        }

        // Dynamic Mouse Interaction (Right-click painting)
        if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
            sf::Vector2i mousePos = sf::Mouse::getPosition(window);
            float minDist = 1e9;
            int bestI = -1;
            int bestJ = -1;

            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    sf::Vector2f sPos = project(-L + i * dx, -L + j * dx, U[i][j]);
                    float dist = std::hypot(sPos.x - mousePos.x, sPos.y - mousePos.y);
                    if (dist < minDist) {
                        minDist = dist;
                        bestI = i;
                        bestJ = j;
                    }
                }
            }

            if (minDist < 40.0f && bestI != -1 && bestJ != -1) {
                for (int di = -2; di <= 2; ++di) {
                    for (int dj = -2; dj <= 2; ++dj) {
                        int ni = bestI + di, nj = bestJ + dj;
                        if (ni >= 0 && ni < N && nj >= 0 && nj < N) 
                            U[ni][nj] = maxInitialTemp;
                    }
                }
            }
        }

        // Processing Partial Differential Solver Steps
        for (int step = 0; step < STEPS_PER_FRAME; ++step) {
            if (!useImplicit) {
                for (int i = 1; i < N - 1; ++i) {
                    for (int j = 1; j < N - 1; ++j) {
                        float laplacian = (U[i+1][j] + U[i-1][j] + U[i][j+1] + U[i][j-1] - 4.0f * U[i][j]) / (dx * dx);
                        nextU[i][j] = U[i][j] + dt * D * laplacian;
                    }
                }
                U = nextU;
            } else {
                float current_dt = dt * 3.0f;
                float r = (D * current_dt) / (2.0f * dx * dx);
                std::vector<std::vector<float>> B(N, std::vector<float>(N, 0.0f));
                for (int i = 1; i < N - 1; ++i) {
                    for (int j = 1; j < N - 1; ++j) {
                        B[i][j] = U[i][j] * (1.0f - 4.0f * r) + r * (U[i+1][j] + U[i-1][j] + U[i][j+1] + U[i][j-1]);
                    }
                }
                nextU = U;
                for (int it = 0; it < 8; ++it) {
                    for (int i = 1; i < N - 1; ++i) {
                        for (int j = 1; j < N - 1; ++j) {
                            nextU[i][j] = (B[i][j] + r * (nextU[i+1][j] + nextU[i-1][j] + nextU[i][j+1] + nextU[i][j-1])) / (1.0f + 4.0f * r);
                        }
                    }
                }
                U = nextU;
            }

            // Compute Selected Boundary Profiles
            for (int k = 0; k < N; ++k) {
                if (isNeumann) {
                    // Left and Right bounds
                    U[k][0]   = U[k][1];
                    U[N-1][k] = U[N-2][k];
                    // Top and Bottom bounds
                    U[0][k]   = U[1][k];
                    U[k][N-1] = U[k][N-2];
                } else {
                    // Left and Right bounds
                    U[k][0]   = 0.0f;
                    U[N-1][k] = 0.0f;
                    // Top and Bottom bounds
                    U[0][k]   = 0.0f;
                    U[k][N-1] = 0.0f;
                }
            }
        }

        // Dynamically update Legend Text Labels based on current maxInitialTemp
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

            // Dynamically update On-Screen HUD Overlay
            std::stringstream hudStream;
            hudStream << "=========================================\n"
                      << "SYSTEM STATUS (HUD):\n"
                      << " [S] Shape: " << (currentShape == CRESCENT ? "Crescent Moon" : currentShape == CYLINDER ? "Cylinder" : "Box") << "\n"
                      << " [B] Boundary: " << (isNeumann ? "Insulated (Neumann)" : "Freezing Cold (Dirichlet)") << "\n"
                      << " [I] Method: " << (useImplicit ? "Implicit (Crank-Nicolson)" : "Explicit (Forward Euler)") << "\n"
                      << " [Up/Down] Max Temp Limit: " << std::fixed << std::setprecision(1) << maxInitialTemp << "\n"
                      << "-----------------------------------------\n"
                      << " 👉 Drag Left Mouse Click to Rotate Camera\n"
                      << " 👉 Drag Right Mouse Click to Inject Heat\n"
                      << "=========================================\n";
            hudText.setString(hudStream.str());
        }

        // Build 3D Mesh Output Vertex Data
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

        // Draw Sequence
        window.clear(sf::Color(15, 15, 18));
        window.draw(mesh);        // Render 3D Surface
        window.draw(legendMesh);  // Render Scale Bar

        // Render Text Scale Labels and HUD
        if (fontLoaded) {
            for (int i = 0; i < numLabels; ++i) {
                window.draw(labels[i]);
            }
            window.draw(hudText);
        }

        window.display();
    }

    return 0;
}
