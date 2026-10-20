#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>


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

// Helper to draw an arrow head with a dynamically variable size and directional inversion
void drawArrowHead(sf::RenderWindow& window, sf::Vector2f tip, float length, float baseSize, sf::Color color) 
{
	if (std::abs(length) < 2.0f || baseSize < 1.0f) return;
	sf::VertexArray triangle(sf::Triangles, 3);
	float direction = (length > 0.0f) ? 1.0f : -1.0f;
	triangle[0].position = tip;
	triangle[1].position = sf::Vector2f(tip.x - direction * baseSize, tip.y - baseSize * 0.4f);
	triangle[2].position = sf::Vector2f(tip.x - direction * baseSize, tip.y + baseSize * 0.4f);
	for (int i = 0; i < 3; ++i) 
	{
		triangle[i].color = color;
	}
	window.draw(triangle);
}

// Thomas Algorithm (TDMA) Solver for Tridiagonal Matrix Systems
void solveTridiagonal(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& result) {
	int mN = d.size();
	std::vector<double> c_prime(mN, 0.0);
	std::vector<double> d_prime(mN, 0.0);
	c_prime[0] = c[0] / b[0];
	d_prime[0] = d[0] / b[0];
	for (int i = 1; i < mN; ++i) 
	{
		double denominator = b[i] - a[i] * c_prime[i - 1];
		if (i < mN - 1) 
		{
			c_prime[i] = c[i] / denominator;
		}
		d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denominator;
	}
	result[mN - 1] = d_prime[mN - 1];
	for (int i = mN - 2; i >= 0; i--) 
	{
		result[i] = d_prime[i] - c_prime[i] * result[i + 1];
	}
}

int main() {
	const unsigned int WINDOW_WIDTH = 1450; // Expanded width to fit all colorbars cleanly
	const unsigned int WINDOW_HEIGHT = 650;
	sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Thermo-Rheological Flow Solver (Multi-Colorbar Extension)");
	window.setFramerateLimit(60);

	// Simulation Grid Parameters
 	const int N = 60;
	double H = 400.0;
	double dy = H / (N - 1);
	double dt = 0.01;

	// Core Rheology Base Parameters
	double K0 = 3500.0;
	double n_index = 1.0;
	double min_visc = 50.0;
	double max_visc = 60000.0;

 	// Thermal & Dissipation Constants
	double alpha_thermal = 800.0;
	double T_wall_init = 20.0;
	double beta = 0.04;

	enum Mode { MODE_VELOCITY, MODE_PRESSURE };
	Mode currentMode = MODE_VELOCITY;
	bool isAdiabatic = false;
	double U_target = 200.0;
	double P_gradient = 0.0;

	// Coupled Storage Arrays
	std::vector<double> u(N, 0.0);
	std::vector<double> u_next(N, 0.0);
	std::vector<double> T(N, T_wall_init);
	std::vector<double> T_next(N, T_wall_init);
	std::vector<double> nu_eff(N, K0);

	// System Matrix Arrays
	std::vector<double> a(N, 0.0), b(N, 1.0), c(N, 0.0), d(N, 0.0);

	// Visual Spatial Metric Layouts
	float velOriginX = 60.0f;
	float viscOriginX = 500.0f;
	float tempOriginX = 940.0f;
	float plateTopY = 160.0f;
	float plateBottomY = plateTopY + static_cast<float>(H);
	float plotWidth = 240.0f;
	float cbWidth = 20.0f;
	float cbHeight = H;

	sf::Font font;
	if (!font.loadFromFile("oswald.ttf")) 
	{
        // Fallback or exit handling
	}

	sf::Text hudText;
	hudText.setFont(font);
	hudText.setCharacterSize(14);
	hudText.setFillColor(sf::Color::White);
	hudText.setPosition(10.0f, 10.0f);

	// Profile Header Labels (Updated colors to match Jet/White theme)
	sf::Text labelVel, labelVisc, labelTemp;
	labelVel.setFont(font);
	labelVel.setCharacterSize(16);
	labelVel.setFillColor(sf::Color::White);
	labelVel.setPosition(velOriginX, plateTopY - 30.0f);
	labelVel.setString("Velocity Profile (u)");
    
	labelVisc.setFont(font);
 	labelVisc.setCharacterSize(16);
	labelVisc.setFillColor(sf::Color::White);
	labelVisc.setPosition(viscOriginX, plateTopY - 30.0f);
	labelVisc.setString("Log Viscosity (log10 nu)");
    
	labelTemp.setFont(font);
	labelTemp.setCharacterSize(16);
	labelTemp.setFillColor(sf::Color::White);
	labelTemp.setPosition(tempOriginX, plateTopY - 30.0f);
	labelTemp.setString("Temperature Field (T)");

	bool spaceWasPressed = false;
	bool tWasPressed = false;

	// Main Execution Frame Loop
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

		// 1. Keyboard Inputs & State Toggles
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) 
		{
			if (!spaceWasPressed) 
			{
				currentMode = (currentMode == MODE_VELOCITY) ? MODE_PRESSURE : MODE_VELOCITY;
				std::fill(u.begin(), u.end(), 0.0);
				std::fill(T.begin(), T.end(), T_wall_init);
				spaceWasPressed = true;
			}
		} 
		else 
		{
			spaceWasPressed = false;
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::T)) 
		{
			if (!tWasPressed) 
			{
				isAdiabatic = !isAdiabatic;
				tWasPressed = true;
			}
		} 
		else 
		{
			tWasPressed = false;
		}

		if (currentMode == MODE_VELOCITY) 
		{
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) U_target += 4.0;
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) U_target -= 4.0;
		} 
		else 
		{
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) P_gradient += 3.0;
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) P_gradient -= 3.0;
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) 
		{
			n_index -= 0.01;
			if (n_index < 0.2) n_index = 0.2;
		}
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) 
		{
			n_index += 0.01;
			if (n_index > 2.0) n_index = 2.0;
		}

 		double u_bottom = 0.0;
		double u_top = (currentMode == MODE_VELOCITY) ? U_target : 0.0;

		// 2. Physics Simulation: Coupled Momentum & Energy Solver Loop
		for (int step = 0; step < 6; ++step) 
		{
			for (int i = 0; i < N; ++i) 
			{
				double shear_rate = 0.0;
				if (i == 0) 
				{
					shear_rate = (u[1] - u[0]) / dy;
				}
				else if (i == N - 1) 
				{
					shear_rate = (u[N - 1] - u[N - 2]) / dy;
				}
				else 
				{
					shear_rate = (u[i + 1] - u[i - 1]) / (2.0 * dy);
				}
				double abs_shear = std::abs(shear_rate);
				if (abs_shear < 1e-5) abs_shear = 1e-5;
		        
				double K_T = K0 * std::exp(-beta * (T[i] - T_wall_init));
				double computed_nu = K_T * std::pow(abs_shear, n_index - 1.0);
				nu_eff[i] = std::max(min_visc, std::min(max_visc, computed_nu));
			}

			for (int i = 1; i < N - 1; ++i) 
			{
				double nu_mid_right = 0.5 * (nu_eff[i] + nu_eff[i + 1]);
				double nu_mid_left = 0.5 * (nu_eff[i] + nu_eff[i - 1]);
				double alpha_right = (nu_mid_right * dt) / (2.0 * dy * dy);
				double alpha_left = (nu_mid_left * dt) / (2.0 * dy * dy);
                
				a[i] = -alpha_left;
				b[i] = 1.0 + alpha_left + alpha_right;
				c[i] = -alpha_right;
				d[i] = alpha_left * u[i - 1] + (1.0 - alpha_left - alpha_right) * u[i] + alpha_right * u[i + 1] - P_gradient * dt;
			}
			a[0] = 0.0; b[0] = 1.0; c[0] = 0.0; d[0] = u_bottom;
			a[N - 1] = 0.0; b[N - 1] = 1.0; c[N - 1] = 0.0; d[N - 1] = u_top;
            
			solveTridiagonal(a, b, c, d, u_next);
			u = u_next;

			double s_thermal = (alpha_thermal * dt) / (2.0 * dy * dy);
			for (int i = 1; i < N - 1; ++i) 
			{
				double shear_rate = (u[i + 1] - u[i - 1]) / (2.0 * dy);
				double viscous_dissipation_source = nu_eff[i] * (shear_rate * shear_rate) * 0.00004;
                
				a[i] = -s_thermal;
				b[i] = 1.0 + 2.0 * s_thermal;
  				c[i] = -s_thermal;
				d[i] = s_thermal * T[i - 1] + (1.0 - 2.0 * s_thermal) * T[i] + s_thermal * T[i + 1] + viscous_dissipation_source * dt;
			}

			if (isAdiabatic) 
			{
				a[0] = 0.0; b[0] = 1.0; c[0] = -1.0; d[0] = 0.0;
				a[N - 1] = -1.0; b[N - 1] = 1.0; c[N - 1] = 0.0; d[N - 1] = 0.0;
			} 
			else 
			{
				a[0] = 0.0; b[0] = 1.0; c[0] = 0.0; d[0] = T_wall_init;
				a[N - 1] = 0.0; b[N - 1] = 1.0; c[N - 1] = 0.0; d[N - 1] = T_wall_init;
			}
            
			solveTridiagonal(a, b, c, d, T_next);
			T = T_next;
		}

		// 3. Render Visualization Profiles
		window.clear(sf::Color(20, 20, 25));
		std::stringstream ss;
		ss << "=== Thermo-Rheological Solver State ===\n"
		<< "Mode: " << (currentMode == MODE_VELOCITY ? "Couette / Velocity-Driven" : "Poiseuille / Pressure-Driven") << " [Space to toggle]\n"
		<< "Thermal Bound: " << (isAdiabatic ? "Adiabatic (Insulated)" : "Isothermal (Fixed Temp)") << " [T to toggle]\n"
		<< "Flow Index (n): " << std::fixed << std::setprecision(2) << n_index << " [Left/Right arrows]\n"
		<< (currentMode == MODE_VELOCITY ? "Target Upper Velocity: " : "Pressure Gradient: ") << (currentMode == MODE_VELOCITY ? U_target : P_gradient) << " [Up/Down arrows]\n"
		<< "Max Temp: " << *std::max_element(T.begin(), T.end()) << " C";

		hudText.setString(ss.str());
		window.draw(hudText);

 		// Draw Channel Borders (Horizontal Bounds)
		sf::VertexArray plates(sf::Lines, 6);
		plates[0] = sf::Vertex(sf::Vector2f(velOriginX, plateTopY), sf::Color(100, 100, 100));
		plates[1] = sf::Vertex(sf::Vector2f(velOriginX + plotWidth, plateTopY), sf::Color(100, 100, 100));
		plates[2] = sf::Vertex(sf::Vector2f(viscOriginX, plateTopY), sf::Color(100, 100, 100));
		plates[3] = sf::Vertex(sf::Vector2f(viscOriginX + plotWidth, plateTopY), sf::Color(100, 100, 100));
		plates[4] = sf::Vertex(sf::Vector2f(tempOriginX, plateTopY), sf::Color(100, 100, 100));
		plates[5] = sf::Vertex(sf::Vector2f(tempOriginX + plotWidth, plateTopY), sf::Color(100, 100, 100));
		window.draw(plates);

		sf::VertexArray bottomPlates(sf::Lines, 6);
		bottomPlates[0] = sf::Vertex(sf::Vector2f(velOriginX, plateBottomY), sf::Color(100, 100, 100));
		bottomPlates[1] = sf::Vertex(sf::Vector2f(velOriginX + plotWidth, plateBottomY), sf::Color(100, 100, 100));
		bottomPlates[2] = sf::Vertex(sf::Vector2f(viscOriginX, plateBottomY), sf::Color(100, 100, 100));
		bottomPlates[3] = sf::Vertex(sf::Vector2f(viscOriginX + plotWidth, plateBottomY), sf::Color(100, 100, 100));
		bottomPlates[4] = sf::Vertex(sf::Vector2f(tempOriginX, plateBottomY), sf::Color(100, 100, 100));
		bottomPlates[5] = sf::Vertex(sf::Vector2f(tempOriginX + plotWidth, plateBottomY), sf::Color(100, 100, 100));
		window.draw(bottomPlates);

		sf::VertexArray leftAxes(sf::Lines, 6);
		leftAxes[0] = sf::Vertex(sf::Vector2f(velOriginX, plateTopY), sf::Color(100, 100, 100));
		leftAxes[1] = sf::Vertex(sf::Vector2f(velOriginX, plateBottomY), sf::Color(100, 100, 100));
		leftAxes[2] = sf::Vertex(sf::Vector2f(viscOriginX, plateTopY), sf::Color(100, 100, 100));
		leftAxes[3] = sf::Vertex(sf::Vector2f(viscOriginX, plateBottomY), sf::Color(100, 100, 100));
		leftAxes[4] = sf::Vertex(sf::Vector2f(tempOriginX, plateTopY), sf::Color(100, 100, 100));
		leftAxes[5] = sf::Vertex(sf::Vector2f(tempOriginX, plateBottomY), sf::Color(100, 100, 100));
		window.draw(leftAxes);

		sf::VertexArray velLines(sf::LineStrip, N);
		sf::VertexArray viscLines(sf::LineStrip, N);
		sf::VertexArray tempLines(sf::LineStrip, N);

		double max_u_val = std::max(10.0, *std::max_element(u.begin(), u.end()));
		double min_u_val = std::min(0.0, *std::min_element(u.begin(), u.end()));
		if (min_u_val > 0.0) min_u_val = 0.0;
		double u_range = (max_u_val - min_u_val == 0.0) ? 1.0 : (max_u_val - min_u_val);

		double logMin = std::log10(min_visc);
 		double logMax = std::log10(max_visc);
		double max_T = std::max(T_wall_init + 5.0, *std::max_element(T.begin(), T.end()));

		sf::VertexArray linesToDraw(sf::Lines);

		for (int i = 0; i < N; ++i) 
		{
			float yPos = plateBottomY - static_cast<float>(i) * (static_cast<float>(H) / (N - 1));

 			// Velocity Node (Updated to JetColor mapping)
			float xVel = velOriginX + static_cast<float>((u[i] - min_u_val) / u_range) * plotWidth;
			velLines[i].position = sf::Vector2f(xVel, yPos);
			sf::Color vColor = getJetColor(u[i], min_u_val, max_u_val);
			velLines[i].color = vColor;
            
			float xVelZero = velOriginX + static_cast<float>((0.0 - min_u_val) / u_range) * plotWidth;
			sf::Color fadedVColor = vColor; fadedVColor.a = 75;
			linesToDraw.append(sf::Vertex(sf::Vector2f(xVelZero, yPos), fadedVColor));
			linesToDraw.append(sf::Vertex(sf::Vector2f(xVel, yPos), fadedVColor));

			// Effective Viscosity Node (Updated to JetColor mapping)
			double logVisc = std::log10(nu_eff[i]);
			float xVisc = viscOriginX + static_cast<float>((logVisc - logMin) / (logMax - logMin)) * plotWidth;
			viscLines[i].position = sf::Vector2f(xVisc, yPos);
			sf::Color viscColor = getJetColor(logVisc, logMin, logMax);
			viscLines[i].color = viscColor;
            
			sf::Color fadedViscColor = viscColor; fadedViscColor.a = 75;
			linesToDraw.append(sf::Vertex(sf::Vector2f(viscOriginX, yPos), fadedViscColor));
			linesToDraw.append(sf::Vertex(sf::Vector2f(xVisc, yPos), fadedViscColor));

			// Temperature Node
			float xTemp = tempOriginX + static_cast<float>((T[i] - T_wall_init) / (max_T - T_wall_init)) * plotWidth;
			tempLines[i].position = sf::Vector2f(xTemp, yPos);
			sf::Color tColor = getJetColor(T[i], T_wall_init, max_T);
			tempLines[i].color = tColor;
            
			sf::Color fadedTColor = tColor; fadedTColor.a = 75;
			linesToDraw.append(sf::Vertex(sf::Vector2f(tempOriginX, yPos), fadedTColor));
			linesToDraw.append(sf::Vertex(sf::Vector2f(xTemp, yPos), fadedTColor));
		}

		window.draw(linesToDraw);
		window.draw(velLines);
		window.draw(viscLines);
		window.draw(tempLines);

		// Vector arrow render pass
		for (int i = 0; i < N; i += 2) 
		{
 			float yPos = plateBottomY - static_cast<float>(i) * (static_cast<float>(H) / (N - 1));
			const float MIN_ARROW_SIZE = 3.0f;
			const float MAX_ARROW_SIZE = 12.0f;

 			// Velocity Vector Color Mapping
			float xVel = velOriginX + static_cast<float>((u[i] - min_u_val) / u_range) * plotWidth;
 			float xVelZero = velOriginX + static_cast<float>((0.0 - min_u_val) / u_range) * plotWidth;
			float velSizeFactor = MIN_ARROW_SIZE;
			if (currentMode == MODE_VELOCITY) 
			{
				double max_possible_u = std::max(1.0, std::abs(U_target));
				double velocity_ratio = std::abs(u[i]) / max_possible_u;
				velSizeFactor = MIN_ARROW_SIZE + static_cast<float>(velocity_ratio) * (MAX_ARROW_SIZE - MIN_ARROW_SIZE);
			} 
			else 
			{
				double p_ratio = std::min(1.0, std::abs(P_gradient) / 150.0);
				velSizeFactor = MIN_ARROW_SIZE + static_cast<float>(p_ratio) * (MAX_ARROW_SIZE - MIN_ARROW_SIZE);
			}
			velSizeFactor = std::max(MIN_ARROW_SIZE, std::min(MAX_ARROW_SIZE, velSizeFactor));
			sf::Color vArrowColor = getJetColor(u[i], min_u_val, max_u_val);
			vArrowColor.a = 200;
			drawArrowHead(window, sf::Vector2f(xVel, yPos), xVel - xVelZero, velSizeFactor, vArrowColor);

			// Viscosity Vector Color Mapping
			double logVisc = std::log10(nu_eff[i]);
			float xVisc = viscOriginX + static_cast<float>((logVisc - logMin) / (logMax - logMin)) * plotWidth;
			float viscRatio = static_cast<float>((logVisc - logMin) / (logMax - logMin));
			float viscSizeFactor = MIN_ARROW_SIZE + viscRatio * (MAX_ARROW_SIZE - MIN_ARROW_SIZE);
			sf::Color viscArrowColor = getJetColor(logVisc, logMin, logMax);
			viscArrowColor.a = 200;
			drawArrowHead(window, sf::Vector2f(xVisc, yPos), xVisc - viscOriginX, viscSizeFactor, viscArrowColor);

			// Temperature Vector Color Mapping
			float xTemp = tempOriginX + static_cast<float>((T[i] - T_wall_init) / (max_T - T_wall_init)) * plotWidth;
			sf::Color tArrowColor = getJetColor(T[i], T_wall_init, max_T);
			tArrowColor.a = 200;
			float tempRatio = static_cast<float>((T[i] - T_wall_init) / (max_T - T_wall_init));
			float tempSizeFactor = MIN_ARROW_SIZE + tempRatio * (MAX_ARROW_SIZE - MIN_ARROW_SIZE);
			drawArrowHead(window, sf::Vector2f(xTemp, yPos), xTemp - tempOriginX, tempSizeFactor, tArrowColor);
		}

		window.draw(labelVel);
		window.draw(labelVisc);
		window.draw(labelTemp);

		// =================================================================
		// COLORBAR LEGEND RENDERING STEP (ALL NOW USING JETCOLOR)
		// =================================================================
		int colorbarSteps = 100;
		float spacingOffset = 30.0f;

		// --- 1. Velocity Colorbar (JetColor) ---
		float cbVelX = velOriginX + plotWidth + spacingOffset;
		sf::VertexArray cbVel(sf::TrianglesStrip, colorbarSteps * 2);
		for (int i = 0; i < colorbarSteps; ++i) 
		{
			float interFactor = static_cast<float>(i) / (colorbarSteps - 1);
			float yBarPos = plateTopY + cbHeight - (interFactor * cbHeight);
			double uValue = min_u_val + interFactor * (max_u_val - min_u_val);
			sf::Color cColor = getJetColor(uValue, min_u_val, max_u_val);
			cbVel[i * 2].position = sf::Vector2f(cbVelX, yBarPos);
			cbVel[i * 2].color = cColor;
			cbVel[i * 2 + 1].position = sf::Vector2f(cbVelX + cbWidth, yBarPos);
			cbVel[i * 2 + 1].color = cColor;
		}
		window.draw(cbVel);
        
		sf::RectangleShape borderVel(sf::Vector2f(cbWidth, cbHeight));
		borderVel.setPosition(cbVelX, plateTopY);
		borderVel.setFillColor(sf::Color::Transparent);
		borderVel.setOutlineColor(sf::Color(150, 150, 150));
		borderVel.setOutlineThickness(1.0f);
		window.draw(borderVel);

		// --- 2. Viscosity Colorbar (JetColor) ---
	 	float cbViscX = viscOriginX + plotWidth + spacingOffset;
		sf::VertexArray cbVisc(sf::TrianglesStrip, colorbarSteps * 2);
		for (int i = 0; i < colorbarSteps; ++i) 
		{
			float interFactor = static_cast<float>(i) / (colorbarSteps - 1);
			float yBarPos = plateTopY + cbHeight - (interFactor * cbHeight);
			double vLogVal = logMin + interFactor * (logMax - logMin);
			sf::Color cColor = getJetColor(vLogVal, logMin, logMax);
			cbVisc[i * 2].position = sf::Vector2f(cbViscX, yBarPos);
			cbVisc[i * 2].color = cColor;
			cbVisc[i * 2 + 1].position = sf::Vector2f(cbViscX + cbWidth, yBarPos);
			cbVisc[i * 2 + 1].color = cColor;
		}
		window.draw(cbVisc);
        
		sf::RectangleShape borderVisc(sf::Vector2f(cbWidth, cbHeight));
		borderVisc.setPosition(cbViscX, plateTopY);
		borderVisc.setFillColor(sf::Color::Transparent);
		borderVisc.setOutlineColor(sf::Color(150, 150, 150));
		borderVisc.setOutlineThickness(1.0f);
		window.draw(borderVisc);

		// --- 3. Temperature Colorbar (JetColor) ---
		float cbTempX = tempOriginX + plotWidth + spacingOffset;
		sf::VertexArray cbTemp(sf::TrianglesStrip, colorbarSteps * 2);
		for (int i = 0; i < colorbarSteps; ++i) 
		{
			float interFactor = static_cast<float>(i) / (colorbarSteps - 1);
			float yBarPos = plateTopY + cbHeight - (interFactor * cbHeight);
			double tValue = T_wall_init + interFactor * (max_T - T_wall_init);
			sf::Color cColor = getJetColor(tValue, T_wall_init, max_T);
			cbTemp[i * 2].position = sf::Vector2f(cbTempX, yBarPos);
			cbTemp[i * 2].color = cColor;
			cbTemp[i * 2 + 1].position = sf::Vector2f(cbTempX + cbWidth, yBarPos);
			cbTemp[i * 2 + 1].color = cColor;
		}
		window.draw(cbTemp);
        
		sf::RectangleShape borderTemp(sf::Vector2f(cbWidth, cbHeight));
		borderTemp.setPosition(cbTempX, plateTopY);
		borderTemp.setFillColor(sf::Color::Transparent);
		borderTemp.setOutlineColor(sf::Color(150, 150, 150));
		borderTemp.setOutlineThickness(1.0f);
		window.draw(borderTemp);

		// --- Labels Construction and Rendering Loop ---
		sf::Text scaleMax, scaleMid, scaleMin;
		scaleMax.setFont(font); scaleMax.setCharacterSize(11); scaleMax.setFillColor(sf::Color::White);
		scaleMid.setFont(font); scaleMid.setCharacterSize(11); scaleMid.setFillColor(sf::Color::White);
		scaleMin.setFont(font); scaleMin.setCharacterSize(11); scaleMin.setFillColor(sf::Color::White);

		// Velocity scale text labels
		std::stringstream sMaxV, sMidV, sMinV;
		sMaxV << std::fixed << std::setprecision(1) << max_u_val << " m/s";
		sMidV << std::fixed << std::setprecision(1) << (min_u_val + (max_u_val - min_u_val) * 0.5) << " m/s";
		sMinV << std::fixed << std::setprecision(1) << min_u_val << " m/s";
		scaleMax.setString(sMaxV.str()); scaleMax.setPosition(cbVelX + cbWidth + 6.0f, plateTopY - 6.0f);
		scaleMid.setString(sMidV.str()); scaleMid.setPosition(cbVelX + cbWidth + 6.0f, plateTopY + (cbHeight * 0.5f) - 6.0f);
		scaleMin.setString(sMinV.str()); scaleMin.setPosition(cbVelX + cbWidth + 6.0f, plateTopY + cbHeight - 6.0f);
		window.draw(scaleMax); window.draw(scaleMid); window.draw(scaleMin);

		// Viscosity scale text labels
		std::stringstream sMaxVi, sMidVi, sMinVi;
		sMaxVi << std::fixed << std::setprecision(2) << logMax;
		sMidVi << std::fixed << std::setprecision(2) << (logMin + (logMax - logMin) * 0.5);
		sMinVi << std::fixed << std::setprecision(2) << logMin;
		scaleMax.setString(sMaxVi.str()); scaleMax.setPosition(cbViscX + cbWidth + 6.0f, plateTopY - 6.0f);
		scaleMid.setString(sMidVi.str()); scaleMid.setPosition(cbViscX + cbWidth + 6.0f, plateTopY + (cbHeight * 0.5f) - 6.0f);
		scaleMin.setString(sMinVi.str()); scaleMin.setPosition(cbViscX + cbWidth + 6.0f, plateTopY + cbHeight - 6.0f);
		window.draw(scaleMax); window.draw(scaleMid); window.draw(scaleMin);

		// Temperature scale text labels
		std::stringstream sMaxT, sMidT, sMinT;
		sMaxT << std::fixed << std::setprecision(1) << max_T << " C";
	 	sMidT << std::fixed << std::setprecision(1) << (T_wall_init + (max_T - T_wall_init) * 0.5) << " C";
		sMinT << std::fixed << std::setprecision(1) << T_wall_init << " C";
		scaleMax.setString(sMaxT.str()); scaleMax.setPosition(cbTempX + cbWidth + 6.0f, plateTopY - 6.0f);
		scaleMid.setString(sMidT.str()); scaleMid.setPosition(cbTempX + cbWidth + 6.0f, plateTopY + (cbHeight * 0.5f) - 6.0f);
		scaleMin.setString(sMinT.str()); scaleMin.setPosition(cbTempX + cbWidth + 6.0f, plateTopY + cbHeight - 6.0f);
		window.draw(scaleMax); window.draw(scaleMid); window.draw(scaleMin);

		window.display();
	}

    return 0;
}
