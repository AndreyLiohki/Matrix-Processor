#include <iostream>
#include <vector>
#include "Matrix.h"
#include "MatrixOperations.h"

int main() {
    std::vector<std::vector<long double>> matrix = {
        {1,2,3},
        {4,5,6},
        {7,8,10},

    };
    
    Matrix matr(matrix);
    std::cout << matr.determinant() << std::endl;
    std::cout << matr.trace() << std::endl;
    std::cout << matr.rank() << std::endl;


    MatrixOperations a;
    std::vector<std::vector<long double>> test = matr.hausholderAlgo();
    Matrix create(test);
    create.printMatrix();
    std::cout << std::endl;
    Matrix obj1 = a.multiplication(create.transposeMatrix(), create);
    Matrix obj2 = a.multiplication(create, create.transposeMatrix());

    obj1.printMatrix();
    std::cout << std::endl;
    obj2.printMatrix();

}