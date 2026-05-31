#include "symintegral/symintegrationc++.h"

using namespace std;

int main() 
{
	cout << "\033[2J\033[1;1H"; // Clear screen
	print_centered(BOLD BLUE "Iris Species Predictor" RESET, '=', 60);
	cout << endl;

	FNN_NoHiddenLayer_Iris_ConjugateGradient nn;
	nn.load_model("model.txt");

	vector<string> species_names = {"Iris-setosa", "Iris-versicolor", "Iris-virginica"};

	while (true) 
	{
		cout << BLUE << "\nEnter four input features (or type 'exit' to quit):" << RESET << endl;
		cout << BLUE << "1. Sepal length (cm): " << RESET;
		std::string input;
		std::getline(std::cin, input);
        
		if (input == "exit") 
		{
			break;
		}

		double x1 = std::stod(input);
		double x2 = get_double_input(BLUE "2. Sepal width (cm):  " RESET);
		double x3 = get_double_input(BLUE "3. Petal length (cm): " RESET);
		double x4 = get_double_input(BLUE "4. Petal width (cm):  " RESET);

		vector<double> inputs = {x1, x2, x3, x4};

		cout << "\n" << BLUE << "Predicting..." << RESET << endl;
		for (int i = 0; i <= 100; ++i) 
		{
			print_progress_bar(i / 100.0);
			delay(); // Use our simple delay function instead of sleep
		}
		cout << endl;

		vector<double> prediction = nn.predict(inputs);
		int predicted_class = std::distance(prediction.begin(), std::max_element(prediction.begin(), prediction.end()));
		
		cout << "\n" << GREEN << "Prediction Results:" << RESET << endl;
		cout << GREEN << "-------------------" << RESET << endl;
		cout << BLUE << "Predicted Species: " << RESET << BOLD << species_names[predicted_class] << RESET << endl;
		cout << BLUE << "\nProbabilities:" << RESET << endl;
		
		for (size_t i = 0; i < species_names.size(); ++i) 
		{
			//cout << " " << std::fixed << std::setprecision(2) << prediction[i] * 100 << "%   " << species_names[i] << std::left << ": " ;

			//cout << std::setw(46) << std::left << std::setw(16)  << ": ";
			//print_progress_bar(prediction[i], 20);
			cout << " " << std::fixed << std::setprecision(2) << prediction[i] * 100 << "%  \t"   << species_names[i] << endl;
		}
		cout << endl;
 	}

	print_centered(BLUE "Thank you for using the Iris Species Predictor. Goodbye!" RESET, '-', 70);

	return 0;
}