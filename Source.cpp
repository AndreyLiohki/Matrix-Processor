#include <iostream>
#include <vector>
#include "Matrix.h"
#include "MatrixOperations.h"

int main() {
    std::vector<std::vector<long double>> matrix = {
        {3, 2, 4, 4, 5},
        {2, 4, 3, 1, 2},
        {1, 1, 5, 6, 3},
        {7, 4, 9, 8, 9},
        {3, 5, 4, 8, 7}
    };
    
    Matrix matr(matrix);
    std::pair<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>> par = matr.luDecomposition();

    MatrixOperations a;
    Matrix matr1(a.multiplication(par.first, par.second));
    matr1.printMatrix();
}