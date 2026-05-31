
#include "symintegrationc++.h"


using namespace std;

int main() 
{
	try 
	{
	std::cout << BLUE << BOLD << "\n======== Iris Classification with Neural Network ========" << RESET << std::endl;
	std::cout << BLUE << "Loading dataset..." << RESET << std::endl;

	std::vector<int> labels;
	std::vector<std::string> label_names;
	vector<vector<double>> dataset = read_dataset_iris("Iris.csv", labels, label_names);

	DEBUG("Dataset loaded. Size: " << dataset.size() << " samples, " << dataset[0].size() << " features");
	DEBUG("Number of labels: " << labels.size());
	DEBUG("Label names: ");
	for (const auto& name : label_names) 
	{
		DEBUG(" - " << name);
	}

	std::vector<size_t> indices(dataset.size());
	std::iota(indices.begin(), indices.end(), 0);
	std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

	//  The dataset is split into 80% training and 20% testing sets
	size_t train_size = dataset.size() * 0.8;
	std::vector<std::vector<double>> X_train, X_test, X_all;
	std::vector<int> y_train, y_test, y_all;
	// before shuffled by std::shuffle, labels[0] to labels[49] has value of 0, representing Setosa
	// labels[indices[i]] has value between 0, 1, and 2 to represent the species
	for (size_t i = 0; i < dataset.size(); ++i) 
	{
		X_all.push_back(dataset[indices[i]]);
		y_all.push_back(labels[indices[i]]);
		if (i < train_size) 
		{
			X_train.push_back(dataset[indices[i]]);
			y_train.push_back(labels[indices[i]]);
		} 
		else 
		{
			X_test.push_back(dataset[indices[i]]);
			y_test.push_back(labels[indices[i]]);
		}
	}
	
	DEBUG("Train set size: " << X_train.size());
	DEBUG("Test set size: " << X_all.size());

	std::cout << BLUE << BOLD << "\nInitializing Neural Network..." << RESET << std::endl;
	FNN_NoHiddenLayer_Iris nn;

	std::cout << BLUE << BOLD << "\nStarting Training Process" << RESET << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	nn.train(X_train, y_train, 1000);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	std::cout << GREEN << BOLD << "\nTraining completed in " << elapsed.count() << " seconds." << RESET << std::endl;

	evaluate_model(nn, X_all, y_all, label_names); // evaluate model on all the data

	nn.save_model("model.txt");
	std::cout << GREEN << BOLD << "\nModel trained and saved successfully!" << RESET << std::endl;

	} catch (const std::exception& e) {
	std::cerr << RED << BOLD << "An error occurred: " << e.what() << RESET << std::endl;
	return 1;
	} // end try

	return 0;
}