#pragma once
#include <vector>
#include <iostream>
#include <iomanip>

class Matrix {
public:

	Matrix(size_t rows, size_t columns) : rows(rows), columns(columns),
		matrix(rows, std::vector<long double>(columns, 0.0)) {
		if (rows == 0 || columns == 0) {
			throw std::invalid_argument("Rows and columns must be greater than 0.");
		}
	}

	Matrix(const std::vector<std::vector<long double>>& inputMatrix) {
		if (inputMatrix.empty()) {
			throw std::invalid_argument("Input data cannot be empty.");
		}
		rows = inputMatrix.size();
		columns = inputMatrix[0].size();

		for (const auto& row : inputMatrix) {
			if (row.size() != columns) {
				throw std::invalid_argument("All rows must have same size.");
			}
		}

		matrix = inputMatrix;
	}

	double get(size_t row, size_t col) const{
		if (row >= rows || col >= columns) {
			throw std::out_of_range("Index out of range.");
		}
		return matrix[row][col];
	}

	void set(size_t row, size_t col, const double value) {
		if (row >= rows || col >= columns) {
			throw std::out_of_range("Index out of range.");
		}
		matrix[row][col] = value;
	}

	size_t getNumberOfRows() const  {
		return rows;
	}
	
	size_t getNumberOfColumns() const {
		return columns;
	}

	void printMatrix() const {
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < columns; ++j) {
				std::cout << std::setw(15) << matrix[i][j];
			}
			std::cout << std::endl;
		}
	}

	std::vector<std::vector<long double>> transposeMatrix() const {
			std::vector<std::vector<long double>> transposedMatrix(columns, std::vector<long double>(rows));
		
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < columns; ++j) {
				transposedMatrix[j][i] = matrix[i][j];
			}
		}

		return transposedMatrix;
	}

	std::pair<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>> luDecomposition() const {
		if (rows != columns) {
			throw std::invalid_argument("LU-decomposition doesn't exist for non-square matrix");
		}
		
		std::vector<std::vector<long double>> L(rows, std::vector<long double>(columns));
		std::vector<std::vector<long double>> U(rows, std::vector<long double>(columns));

		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < columns; ++j) {
				L[i][j] = 0;
				U[i][j] = 0;
			}
			L[i][i] = 1;
		}

		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < columns; ++j) {

				if (i <= j) {
					long double summ = 0;
					for (size_t k = 0; k < i; ++k) {
						summ += L[i][k] * U[k][j];
					}
					U[i][j] = matrix[i][j] - summ;
				}

				if (i > j) {
					long double summ = 0;
					for (size_t k = 0; k < j; ++k) {
						summ += L[i][k] * U[k][j];
					}
					if (abs(U[j][j]) > EPSILON)
						L[i][j] = (matrix[i][j] - summ) / U[j][j];
					else
						throw std::invalid_argument("Matrix is singular. LU-decomposition doesn't exist.");
				}
			}
		}

		return { L, U };
	}

	long double determinant() const {
		if (rows != columns) {
			throw std::invalid_argument("LU-decomposition doesn't exist for non-square matrix");
		}

		std::pair<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>> decomposition = luDecomposition();

		long double det = 1;
		for (size_t i = 0; i < rows; ++i) {
			det *= decomposition.first[i][i] * decomposition.second[i][i];
		}
		return det;
	}

	long rank() const{
		long rank = 0;

		std::vector<std::vector<long double>> computeRank = matrix;

		for (size_t row = 0; row < rows; ++row) {
			if (computeRank[row][rank] != 0) {
				for (size_t i = row + 1; i < rows; ++i) {
					long double ratio = computeRank[i][rank] / computeRank[row][rank];
					for (size_t j = rank; j < columns; ++j) {
						computeRank[i][j] -= ratio * computeRank[row][j];
					}
				}
				rank++;
			}
			else {
				bool nonZeroFound = false;
				for (size_t i = row + 1; i < rows; ++i) {
					if (computeRank[i][rank] != 0) {
						std::swap(computeRank[row], computeRank[i]);
						nonZeroFound = true;
						break;
					}
				}
				if (nonZeroFound) {
					for (size_t i = row + 1; i < rows; ++i) {
						long double ratio = computeRank[i][rank] / computeRank[row][rank];
						for (size_t j = rank; j < columns; ++j) {
							computeRank[i][j] -= ratio * computeRank[row][j];
						}
					}
					rank++;
				}
			}
		}
		return rank;
	}

	long double trace() const {
		if (rows != columns) {
			throw std::invalid_argument("LU-decomposition doesn't exist for non-square matrix");
		}

		long double tr = 0;

		for (size_t i = 0; i < rows; ++i) {
			tr += matrix[i][i];
		}

		return tr;
	}

	//std::vector<std::vector<long double>> orthogonalMatrix() const {
	//	std::vector<std::vector<long double>> result(rows, std::vector<long double>(columns));

	//	for (size_t i = 0; i < columns; ++i) {
	//		std::vector<long double> v (rows, 0);
	//		for (size_t j = 0; j < rows; ++j) {
	//			v[j] = matrix[j][i];
	//		}
	//		for (size_t j = 0; j < i; ++j) {
	//			std::vector<long double> proj = vectorProjection(result[j], v);
	//			v = vectorSubstraction(v, proj);
	//		}
	//		result[i] = normalization(v);
	//	}
	//	return result;

	//}

private:
	const long double EPSILON = 1e-12;

	std::vector<std::vector<long double>> matrix;
	size_t rows;
	size_t columns;

	long double scalarProduct(const std::vector<long double>& a, const std::vector<long double>& b) const{
		double value = 0;
		for (size_t i = 0; i < a.size(); ++i) {
			value += a[i] * b[i];
		}
		return value;
	}

	std::vector<long double> vectorProjection(const std::vector<long double>& a, const std::vector<long double>& b) const{
		long double coeff = scalarProduct(a, b) / scalarProduct(b, b);
		std::vector<long double> projection(a.size());
		for (size_t i = 0; i < a.size(); ++i) {
			projection[i] = coeff * b[i];
		}
		return projection;
	}

	std::vector<long double> vectorSubstraction(const std::vector<long double>& a, const std::vector<long double>& b) const {
		std::vector<long double> result(a.size());
		for (size_t i = 0; i < a.size(); ++i) {
			result[i] = a[i] - b[i];
		}
		return result;
	}

	std::vector<long double> normalization(const std::vector<long double>& a) const  {
		long double norm = std::sqrt(scalarProduct(a, a));
		std::vector<long double> result(a.size(), 0);
		for (size_t i = 0; i < a.size(); ++i) {
			result[i] = a[i] / norm;
		}
		return result;
	}

};
