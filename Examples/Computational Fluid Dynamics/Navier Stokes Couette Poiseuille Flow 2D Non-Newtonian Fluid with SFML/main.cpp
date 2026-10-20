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

// Thomas Algorithm (TDMA) Solver for Tridiagonal Matrix Systems
void solveTridiagonal(const std::vector<double>& a, 
                      const std::vector<double>& b, 
                      const std::vector<double>& c, 
                      std::vector<double>& d, 
                      std::vector<double>& u_next) 
{
	int N = d.size();
	std::vector<double> c_prime(N, 0.0);
	std::vector<double> d_prime(N, 0.0);

	for (int i = 1; i < N; ++i) 
	{
		c_prime[i] = c[i] / b[i];
		d_prime[i] = d[i] / b[i];
	}

	for (int i = 1; i < N; ++i) 
	{
		double denominator = b[i] - a[i] * c_prime[i - 1];
		if (i < N - 1) 
		{
			c_prime[i] = c[i] / denominator;
		}
		d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denominator;
	}

	u_next[N - 1] = d_prime[N - 1];
	for (int i = N - 2; i >= 0; --i) 
	{
		u_next[i] = d_prime[i] - c_prime[i] * u_next[i + 1];
	}
}

int main() 
{
	const unsigned int WINDOW_WIDTH = 1100;
	const unsigned int WINDOW_HEIGHT = 800;
	sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Couette Flow Non-Newtonian Rheology Dashboard");
	window.setFramerateLimit(60);

	const int N = 60;                  
	double H = 400.0;                  
	double dy = H / (N - 1);           
	double dt = 0.01;                  
    
	double K = 2000.0;                 
	double n_index = 1.0;              
	double min_visc = 100.0;           
	double max_visc = 50000.0;         

	enum Mode { MODE_VELOCITY, MODE_PRESSURE };
	Mode currentMode = MODE_VELOCITY;

	double U_target = 200.0;            
	double P_gradient = 0.0;            

	std::vector<double> u(N, 0.0);      
	std::vector<double> u_next(N, 0.0); 
	std::vector<double> nu_eff(N, K);   

	std::vector<double> a(N, 0.0); 
	std::vector<double> b(N, 1.0); 
	std::vector<double> c(N, 0.0); 
	std::vector<double> d(N, 0.0); 

	double velOriginX  = 120.0f;  
	double viscOriginX = 620.0f;  
	double plateTopY   = 150.0f;
	double plateBottomY = plateTopY + static_cast<double>(H);
	double plotWidth    = 300.0f; 

	sf::Font font;
	if (!font.loadFromFile("oswald.ttf")) 
	{

	}
    
	sf::Text hudText;
	hudText.setFont(font);
	hudText.setCharacterSize(15);
	hudText.setFillColor(sf::Color::White);
	hudText.setPosition(10.0f, 600.0f);

	bool spaceWasPressed = false;

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

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) 
		{
			if (!spaceWasPressed) 
			{
				currentMode = (currentMode == MODE_VELOCITY) ? MODE_PRESSURE : MODE_VELOCITY;
				std::fill(u.begin(), u.end(), 0.0);
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
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))    P_gradient += 2.0; 
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))  P_gradient -= 2.0;
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) 
		{
			n_index -= 0.01; if (n_index < 0.2) n_index = 0.2;
		}
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) 
		{
			n_index += 0.01; if (n_index > 2.0) n_index = 2.0;
		}

		double u_bottom = 0.0;
		double u_top = (currentMode == MODE_VELOCITY) ? U_target : 0.0;

		for (int step = 0; step < 5; ++step) 
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

			double computed_nu = K * std::pow(abs_shear, n_index - 1.0);
			nu_eff[i] = std::max(min_visc, std::min(max_visc, computed_nu));
			}

			for (int i = 1; i < N - 1; ++i) 
			{
				double nu_mid_right = 0.5 * (nu_eff[i] + nu_eff[i + 1]);
				double nu_mid_left  = 0.5 * (nu_eff[i] + nu_eff[i - 1]);

				double alpha_right = (nu_mid_right * dt) / (2.0 * dy * dy);
				double alpha_left  = (nu_mid_left * dt) / (2.0 * dy * dy);

				a[i] = -alpha_left;
				b[i] = 1.0 + alpha_left + alpha_right;
				c[i] = -alpha_right;
    				d[i] = alpha_left * u[i - 1] + (1.0 - alpha_left - alpha_right) * u[i] + alpha_right * u[i + 1] - P_gradient * dt;
			}

			b[0] = 1.0; c[0] = 0.0; d[0] = u_bottom;
 			a[N - 1] = 0.0; b[N - 1] = 1.0; d[N - 1] = u_top;

			solveTridiagonal(a, b, c, d, u_next);
			u = u_next;
		}

		window.clear(sf::Color(12, 14, 18));

		auto drawPlates = [&](double originX) 
		{
			sf::RectangleShape topPlate(sf::Vector2f(plotWidth, 4.0f));
			topPlate.setPosition(originX, plateTopY - 2.0f);
			topPlate.setFillColor(sf::Color(140, 140, 140));
			window.draw(topPlate);

			sf::RectangleShape bottomPlate(sf::Vector2f(plotWidth, 4.0f));
			bottomPlate.setPosition(originX, plateBottomY - 2.0f);
			bottomPlate.setFillColor(sf::Color(140, 140, 140));
			window.draw(bottomPlate);

			sf::RectangleShape zeroLine(sf::Vector2f(1.0f, H));
			zeroLine.setPosition(originX, plateTopY);
			zeroLine.setFillColor(sf::Color(60, 60, 70, 120));
			window.draw(zeroLine);
		};
		drawPlates(velOriginX);
		drawPlates(viscOriginX);

		double maxVelocityInField = 0.0;
		for (double val : u) 
		{ 
			if (std::abs(val) > maxVelocityInField) 
			{
				maxVelocityInField = std::abs(val); 
			}
		}
		double vMaxColor = std::max(maxVelocityInField, 150.0);
		double vMinColor = -vMaxColor;

		bool hasNegative = false;
		for (double val : u) 
		{ 
			if (val < -0.1) 
			{
				hasNegative = true; 
			}
		}
		if (!hasNegative) vMinColor = 0.0;

		double maxViscInField = *std::max_element(nu_eff.begin(), nu_eff.end());
		double minViscInField = *std::min_element(nu_eff.begin(), nu_eff.end());

		sf::VertexArray velCurve(sf::LineStrip, N);
		for (int i = 0; i < N; ++i) 
		{
			double yPos = plateBottomY - static_cast<double>(i * dy);
			double xVelocityOffset = static_cast<double>(u[i]);

			sf::Color vColor = getJetColor(u[i], vMinColor, vMaxColor);

			sf::Vertex line[] = {
				sf::Vertex(sf::Vector2f(velOriginX, yPos), sf::Color(45, 45, 55, 100)),
				sf::Vertex(sf::Vector2f(velOriginX + xVelocityOffset, yPos), vColor)
				};
			window.draw(line, 2, sf::Lines);

			velCurve[i].position = sf::Vector2f(velOriginX + xVelocityOffset, yPos);
			velCurve[i].color = vColor;
		}
		window.draw(velCurve);

		sf::VertexArray viscCurve(sf::LineStrip, N);
		for (int i = 0; i < N; ++i) 
		{
			double yPos = plateBottomY - static_cast<double>(i * dy);
            
			double normVisc = 0.0f;
  			if (maxViscInField != minViscInField) 
			{
				normVisc = static_cast<double>((nu_eff[i] - minViscInField) / (maxViscInField - minViscInField));
			}
			double xViscOffset = normVisc * plotWidth;

			sf::Color viscColor = getJetColor(nu_eff[i], minViscInField, maxViscInField);

			sf::Vertex line[] = {
				sf::Vertex(sf::Vector2f(viscOriginX, yPos), sf::Color(45, 45, 55, 100)),
				sf::Vertex(sf::Vector2f(viscOriginX + xViscOffset, yPos), viscColor)
 			};
			window.draw(line, 2, sf::Lines);

			viscCurve[i].position = sf::Vector2f(viscOriginX + xViscOffset, yPos);
			viscCurve[i].color = viscColor;
		}
		window.draw(viscCurve);

		auto drawColorbarLegend = [&](double legendX, double minVal, double maxVal, std::string title, std::string units) 
		{
			double legendY = plateTopY;
			double legendW = 20.0f;
			double legendH = H;

			sf::VertexArray bar(sf::Quads, 4 * 100);
			double segH = legendH / 100.0f;
			for (int i = 0; i < 100; ++i) 
			{
				double yTop = legendY + legendH - (i + 1) * segH;
				double yBot = legendY + legendH - i * segH;
				double val = minVal + (i / 99.0) * (maxVal - minVal);
				sf::Color c = getJetColor(val, minVal, maxVal);

				bar[i * 4 + 0] = sf::Vertex(sf::Vector2f(legendX, yBot), c);
				bar[i * 4 + 1] = sf::Vertex(sf::Vector2f(legendX + legendW, yBot), c);
				bar[i * 4 + 2] = sf::Vertex(sf::Vector2f(legendX + legendW, yTop), c);
    				bar[i * 4 + 3] = sf::Vertex(sf::Vector2f(legendX, yTop), c);
  			}
			window.draw(bar);

  			sf::Text label; label.setFont(font); label.setCharacterSize(12); label.setFillColor(sf::Color::White);
            
 			label.setString(std::to_string(static_cast<int>(maxVal)) + " " + units);
			label.setPosition(legendX + legendW + 6.0f, legendY); window.draw(label);

			label.setString(std::to_string(static_cast<int>(minVal)) + " " + units);
			label.setPosition(legendX + legendW + 6.0f, legendY + legendH - 14.0f); window.draw(label);

			label.setString(title); label.setPosition(legendX, legendY - 20.0f); window.draw(label);
		};

		drawColorbarLegend(velOriginX + plotWidth + 20.0f, vMinColor, vMaxColor, "Velocity profile:", "px/s");
		drawColorbarLegend(viscOriginX + plotWidth + 20.0f, minViscInField, maxViscInField, "Viscosity profile (nu):", "m2/s");

		std::string fluidType = "Newtonian Fluid";
		if (n_index < 0.99) fluidType = "Pseudoplastic (Shear-Thinning)";
		else if (n_index > 1.01) fluidType = "Dilatant (Shear-Thickening)";

		std::string modeString = (currentMode == MODE_VELOCITY) ? "VELOCITY DRIVEN (Couette Flow)" : "PRESSURE DRIVEN (Poiseuille Flow)";
		std::string controlString = (currentMode == MODE_VELOCITY) ? 
		"Top Plate Speed (U): " + std::to_string(static_cast<int>(U_target)) + " px/s" :
		"Pressure Gradient (dP/dx): " + std::to_string(static_cast<int>(P_gradient));

		std::stringstream ss; ss << std::fixed << std::setprecision(2) << n_index;

		hudText.setString("Dual-Profile Interactive Rheology Dashboard [Implicit CN Solver]\n"
			"----------------------------------------------------------------------------------------\n"
			"\nPress [SPACE] to Toggle Control Mode | Active Mode: " + modeString + " " +
			"\nUse UP / DOWN Arrows to Adjust: " + controlString + " " +
			"\nUse LEFT / RIGHT Arrows to Tune Fluid Behavior Index (n): " + ss.str() + " [" + fluidType + "]");
		window.draw(hudText);

		window.display();
	}

    return 0;
}
