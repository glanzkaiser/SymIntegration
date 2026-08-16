// http://yang.amp.i.kyoto-u.ac.jp/~yyama/computer/FAQ/eigen/eigensystem.html
// https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

// g++ -o result main.cpp 

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/QR>

using namespace Eigen;

int main() {
    // Define a matrix A (example: 3x3 matrix)
    MatrixXd A(3, 3);
    A << 12, -51, 4,
         6, 167, -68,
        -4, 24, -41;

    // Perform the QR decomposition using ColPivHouseholderQR
    // ColPivHouseholderQR provides column pivoting, which is more robust
    ColPivHouseholderQR<MatrixXd> qr(A);

    // Get the Q and R matrices
    MatrixXd Q = qr.householderQ();
    MatrixXd R = qr.matrixR().triangularView<Upper>();

    std::cout << "Original Matrix A:\n" << A << std::endl;
    std::cout << "\nMatrix Q:\n" << Q << std::endl;
    std::cout << "\nMatrix R (Upper Triangular):\n" << R << std::endl;

    // Verify the decomposition (A should be close to Q * R)
    MatrixXd A_reconstructed = Q * R;
    std::cout << "\nReconstructed A (Q * R):\n" << A_reconstructed << std::endl;

    return 0;
}

