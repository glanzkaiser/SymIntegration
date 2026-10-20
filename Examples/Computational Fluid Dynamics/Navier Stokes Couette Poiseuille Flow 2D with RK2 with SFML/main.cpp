#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

// Jet Colormap Interpolation Helper
sf::Color getJetColor(double v, double vmin, double vmax) 
{
	double dv = vmax - vmin;
	double x = 0.0;
	if (dv != 0.0) 
	{
		x = (v - vmin) / dv;
	}
	x = std::max(0.0, std::min(1.0, x));

	double r = std::max(0.0, std::min(1.0, double(4.0 * x - 1.5)));
	double g = std::max(0.0, std::min(1.0, double(4.0 * x - 0.5))) - std::max(0.0, std::min(1.0, double(4.0 * x - 3.5)));
	double b = std::max(0.0, std::min(1.0, double(1.5 - 4.0 * x))) + std::max(0.0, std::min(1.0, double(4.0 * x - 3.5))) * 0.5;

	return sf::Color(
	static_cast<sf::Uint8>(r * 255),
	static_cast<sf::Uint8>(g * 255),
	static_cast<sf::Uint8>(b * 255)
	);
}

// Computes spatial second derivative and includes the source term for pressure gradient
// Equation: du/dt = nu * d2u/dy2 - (1/rho)*(dp/dx)
std::vector<double> computeRHS(const std::vector<double>& u, double dy, double nu, double pressureGradient) 
{
	int N = u.size();
	std::vector<double> rhs(N, 0.0);
    
	// Central difference for interior points (2nd-order spatial accuracy)
	for (int i = 1; i < N - 1; ++i) 
	{
		double diffusion = (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (dy * dy);
		rhs[i] = nu * diffusion - pressureGradient; 
	}
	return rhs;
}

int main() {
	const unsigned int WINDOW_WIDTH = 1000;
	const unsigned int WINDOW_HEIGHT = 700;
	sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Optimized Flow Simulation (RK2)");
	window.setFramerateLimit(60);

	// Simulation Grid Parameters
	const int N = 60;                  // Increased grid density for better spatial accuracy
	double H = 400.0;                  
	double dy = H / (N - 1);           
    
	// Physical Parameters
	double nu = 2000.0;                
	double dt = 0.0005;                // Smaller time-step to preserve RK2 precision stability bounds
    
	// Simulation Modes
	enum Mode { MODE_VELOCITY, MODE_PRESSURE };
	Mode currentMode = MODE_VELOCITY;

	// Control Variables
	double U_target = 200.0;            // Controls top plate speed in Velocity Mode
	double P_gradient = 0.0;            // Controls pressure force source term in Pressure Mode

	// Velocity Arrays
	std::vector<double> u(N, 0.0);      
    
	// Visual layout metrics
	double paddingX = 120.0f;
	double plateTopY = 100.0f;
	double plateBottomY = plateTopY + static_cast<double>(H);

	sf::Font font;
	if (!font.loadFromFile("oswald.ttf")) 
	{
	// Fallback or system tracking
	}
    
	sf::Text hudText;
	hudText.setFont(font);
	hudText.setCharacterSize(16);
	hudText.setFillColor(sf::Color::White);
	hudText.setPosition(10.0f, 510.0f);

	bool spaceWasPressed = false;

	// Main Loop
	while (window.isOpen()) 
	{
		sf::Event event;
		while (window.pollEvent(event)) 
		{
			if (event.type == sf::Event::Closed)
			{
				window.close();
			}
		}
		// 1. Handle Keyboard Toggles & Parameters
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) 
		{
			if (!spaceWasPressed) 
			{ // Edge triggered toggle to prevent rapid fluttering
				currentMode = (currentMode == MODE_VELOCITY) ? MODE_PRESSURE : MODE_VELOCITY;
				spaceWasPressed = true;
			}
		} 
		else 
		{
			spaceWasPressed = false;
		}

		if (currentMode == MODE_VELOCITY) 
		{
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))    U_target += 4.0;
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))  U_target -= 4.0;
		} 
		else 
		{
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))    P_gradient += 0.5; 
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))  P_gradient -= 0.5;
		}

		// Boundary Conditions (Always enforced explicitly before stepping)
		u[0] = 0.0;          // Fixed bottom wall
		u[N - 1] = (currentMode == MODE_VELOCITY) ? U_target : 0.0; // In pure pressure mode, walls are static

		// 2. Physics Update: Runge-Kutta 2nd Order Solver (Heun's Method)
		// Sub-stepping inside the frame to allow fast physical propagation speeds safely
		for (int step = 0; step < 15; ++step) 
		{
			// Predictor Step (Euler Step to find intermediate value u_pred)
			std::vector<double> rhs1 = computeRHS(u, dy, nu, P_gradient);
			std::vector<double> u_pred(N, 0.0);
			for (int i = 1; i < N - 1; ++i) 
			{
				u_pred[i] = u[i] + dt * rhs1[i];
			}
			// Enforce BCs on predictor
			u_pred[0] = u[0];
			u_pred[N - 1] = u[N - 1];

			// Corrector Step (Averages current gradient with predicted boundary gradient)
			std::vector<double> rhs2 = computeRHS(u_pred, dy, nu, P_gradient);
			for (int i = 1; i < N - 1; ++i) 
			{
				u[i] = u[i] + 0.5 * dt * (rhs1[i] + rhs2[i]);
			}
		}

		// --- Drawing ---
		window.clear(sf::Color(12, 12, 16));

		// Draw Plates
		sf::RectangleShape topPlate(sf::Vector2f(500.0f, 6.0f));
		topPlate.setPosition(paddingX, plateTopY - 3.0f);
		topPlate.setFillColor(sf::Color(180, 180, 180));
		window.draw(topPlate);

		sf::RectangleShape bottomPlate(sf::Vector2f(500.0f, 6.0f));
		bottomPlate.setPosition(paddingX, plateBottomY - 3.0f);
		bottomPlate.setFillColor(sf::Color(180, 180, 180));
		window.draw(bottomPlate);

		// Calculate dynamic color mapping ranges to avoid cutoff boundaries
		double maxVelocityInField = 0.0;
		for (double val : u) 
		{
			if (std::abs(val) > maxVelocityInField) maxVelocityInField = std::abs(val);
		}
		double vMaxColor = std::max(maxVelocityInField, 150.0);
		double vMinColor = -vMaxColor;

		// If field flow is entirely positive, map colorbar zero to max for optimized contrast
		bool hasNegative = false;
		for (double val : u) 
		{ 
			if (val < -0.1) hasNegative = true; 
		}
		if (!hasNegative) vMinColor = 0.0;

		// Draw Vectors & Contour Profile
		sf::VertexArray profileCurve(sf::LineStrip, N);
		for (int i = 0; i < N; ++i) 
		{
			double yPos = plateBottomY - static_cast<double>(i * dy);
			double xVelocityOffset = static_cast<double>(u[i]);

			sf::Color nodeColor = getJetColor(u[i], vMinColor, vMaxColor);

			sf::Vertex line[] = {
				sf::Vertex(sf::Vector2f(paddingX, yPos), sf::Color(45, 45, 55, 120)),
				sf::Vertex(sf::Vector2f(paddingX + xVelocityOffset, yPos), nodeColor)
				};
			window.draw(line, 2, sf::Lines);

			profileCurve[i].position = sf::Vector2f(paddingX + xVelocityOffset, yPos);
			profileCurve[i].color = nodeColor;
		}
		window.draw(profileCurve);

		// Draw Vertical Colormap Scale Legend
		double legendX = 720.0f;
		double legendY = 150.0f;
		double legendWidth = 30.0f;
		double legendHeight = 300.0f;

		sf::VertexArray legendBar(sf::Quads, 4 * 100);
		double segmentHeight = legendHeight / 100.0f;
		for (int i = 0; i < 100; ++i) 
		{
			double yTop = legendY + legendHeight - (i + 1) * segmentHeight;
			double yBot = legendY + legendHeight - i * segmentHeight;
            
			double val = vMinColor + (i / 99.0) * (vMaxColor - vMinColor);
			sf::Color c = getJetColor(val, vMinColor, vMaxColor);

			legendBar[i * 4 + 0] = sf::Vertex(sf::Vector2f(legendX, yBot), c);
			legendBar[i * 4 + 1] = sf::Vertex(sf::Vector2f(legendX + legendWidth, yBot), c);
			legendBar[i * 4 + 2] = sf::Vertex(sf::Vector2f(legendX + legendWidth, yTop), c);
			legendBar[i * 4 + 3] = sf::Vertex(sf::Vector2f(legendX, yTop), c);
		}
		window.draw(legendBar);

		// Draw Text Labels for Legend
		sf::Text labelText;
		labelText.setFont(font);
		labelText.setCharacterSize(13);
		labelText.setFillColor(sf::Color::White);

		labelText.setString(std::to_string(static_cast<int>(vMaxColor)) + " px/s");
		labelText.setPosition(legendX + legendWidth + 10.0f, legendY);
		window.draw(labelText);

		labelText.setString(std::to_string(static_cast<int>(vMinColor)) + " px/s");
		labelText.setPosition(legendX + legendWidth + 10.0f, legendY + legendHeight - 15.0f);
		window.draw(labelText);

		labelText.setString("Velocity Scale:");
		labelText.setPosition(legendX, legendY - 22.0f);
		window.draw(labelText);

		// Draw Configuration HUD Text Info
		std::string modeString = (currentMode == MODE_VELOCITY) ? "VELOCITY DRIVEN (Couette Flow)" : "PRESSURE DRIVEN (Poiseuille Flow)";
		std::string controlString = (currentMode == MODE_VELOCITY) ? 
			"Top Plate Speed (U): " + std::to_string(static_cast<int>(U_target)) + " px/s" :
			"Pressure Gradient (dP/dx): " + std::to_string(P_gradient);

		hudText.setString("Couette/Poiseuille Fluid Solver [RK2 Accuracy]\n"
			"----------------------------------------------------\n"
			"Press [SPACE] to Toggle Control Mode\n"
			"Active Mode: " + modeString + "\n" +
			"Use UP / DOWN Arrows to Adjust Param\n" + controlString);
		window.draw(hudText);

		window.display();
	}

    return 0;
}