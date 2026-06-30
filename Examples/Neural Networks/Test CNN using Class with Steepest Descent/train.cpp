
#include "symintegrationc++.h"


using namespace std;

int main() 
{
	try 
	{
	std::cout << BLUE << BOLD << "\n======== Convolutional Neural Network ========" << RESET << std::endl;
	std::cout << BLUE << "Loading dataset..." << RESET << std::endl;

	vector<vector<vector<int>>> images = loadIDX3("train-images-idx3-ubyte");	
	vector<int> labels = loadIDX1("train-labels-idx1-ubyte");
        
	cout << "Successfully loaded " << images.size() << " images!" << endl;
	cout << "Dimension of Image : " << images[0].size() << "x" << images[0][0].size() << endl;

	DEBUG("Dataset loaded. Size: " << images.size() << " images.");
	DEBUG("Number of labels: " << labels.size());

	vector<int> indices(images.size());
	std::iota(indices.begin(), indices.end(), 0);
	std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

	//  The dataset is split into 80% training and 20% testing sets
	int train_size = images.size() * 0.8;
	vector<vector<vector<int>>> X_train, X_test, X_all;
	vector<int> y_train, y_test, y_all;
	
	// labels[indices[i]] has value between 0, 1, and 2 to represent the species
	for (int i = 0; i < images.size(); ++i) 
	{
		X_all.push_back(images[indices[i]]);
		y_all.push_back(labels[indices[i]]);
		if (i < train_size) 
		{
			X_train.push_back(images[indices[i]]);
			y_train.push_back(labels[indices[i]]);
		} 
		else 
		{
			X_test.push_back(images[indices[i]]);
			y_test.push_back(labels[indices[i]]);
		}
	}
	
	DEBUG("Train set size: " << X_train.size());
	DEBUG("Test set size: " << X_test.size());
	save3DMatrixint(X_train,"shuffledXtrain.txt");
	save3DMatrixint(X_test,"shuffledXtest.txt");
	saveVectorint(y_train,"shuffledytrain.txt");
	saveVectorint(y_test,"shuffledytest.txt");
	
	std::cout << BLUE << BOLD << "\nInitializing Convolutional Neural Network..." << RESET << std::endl;
	CNN_LeNet5 cnn;
	vector<double> prediction = cnn.predict(X_train[0]);


	} catch (const std::exception& e) {
	std::cerr << RED << BOLD << "An error occurred: " << e.what() << RESET << endl;
	return 1;
	} // end try

	return 0;
}