#include <iostream>
#include <vector>
#include "Matrix.h"
#include "MatrixOperations.h"

int main() {
    std::vector<std::vector<long double>> matrix(7, std::vector<long double>(7));

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            matrix[i][j] = rand() % 9 + 1; 
        }
    }
    
    Matrix matr(matrix);
    std::cout << matr.determinant() << std::endl;
    std::cout << matr.trace() << std::endl;
    std::cout << matr.rank() << std::endl;


    MatrixOperations a;
    std::vector<std::vector<long double>> test1 = matr.hausholderAlgo().first;
    std::vector<std::vector<long double>> test2 = matr.hausholderAlgo().second;
    Matrix A(test1);
    Matrix Q(test2);
    A.printMatrix();
    std::cout << std::endl;

    Q.printMatrix();
    std::cout << std::endl;

    Matrix result(a.multiplication(Q, A));

    result.printMatrix();
    std::cout << std::endl;

    Matrix isOrthogonal(a.multiplication(Q, Q.transposeMatrix()));
    isOrthogonal.printMatrix();


}