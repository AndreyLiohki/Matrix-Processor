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
	
	std::pair<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>> hausholderAlgo() {
		std::vector<std::vector<long double>> A = matrix;
		std::vector<std::vector<long double>> Q(A.size(), std::vector<long double>(A.size(), 0));
		size_t m = A.size();
		size_t n = A[0].size();

		for (size_t i = 0; i < m; ++i) {
			Q[i][i] = 1;
		}

		for (size_t k = 0; k < m-1; ++k) {
			std::vector<long double> x(m - k);
			for (size_t i = k; i < m; ++i) {
				x[i - k] = A[i][k];
			}

			std::vector<long double> v = hausholderVector(x);

			std::vector<std::vector<long double>> Hk(m, std::vector<long double>(m, 0));
			for (size_t i = 0; i < m; ++i) {
				Hk[i][i] = 1;
			}

			for (size_t i = 0; i < m - k; ++i) {
				for (size_t j = 0; j < m - k; ++j) {
					Hk[k + i][k + j] -= 2.0 * v[i] * v[j];
				}
			}

			Q = stlMatrixMultiplication(Q, Hk);

			A = stlMatrixMultiplication(Hk, A);
		}

		for (size_t i = 0; i < Q.size(); ++i) {
			for (size_t j = 0; j < Q[0].size(); ++j) {
				if (std::abs(Q[i][j]) < EPSILON) {
					Q[i][j] = 0;
				}
			}
		}

		for (size_t i = 0; i < A.size(); ++i) {
			for (size_t j = 0; j < A[0].size(); ++j) {
				if (std::abs(A[i][j]) < EPSILON) {
					A[i][j] = 0;
				}
			}
		}

		return {A, Q};
	}

	std::pair<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>> hausholderAlgo(std::vector<std::vector<long double>> Ak) {
		std::vector<std::vector<long double>> A = Ak;
		std::vector<std::vector<long double>> Q(Ak.size(), std::vector<long double>(Ak.size(), 0));
		size_t m = Ak.size();
		size_t n = Ak[0].size();

		for (size_t i = 0; i < m; ++i) {
			Q[i][i] = 1;
		}

		for (size_t k = 0; k < m-1; ++k) {
			std::vector<long double> x(m - k);
			for (size_t i = k; i < m; ++i) {
				x[i - k] = Ak[i][k];
			}

			std::vector<long double> v = hausholderVector(x);

			std::vector<std::vector<long double>> Hk(m, std::vector<long double>(m, 0));
			for (size_t i = 0; i < m; ++i) {
				Hk[i][i] = 1;
			}

			for (size_t i = 0; i < m - k; ++i) {
				for (size_t j = 0; j < m - k; ++j) {
					Hk[k + i][k + j] -= 2.0 * v[i] * v[j];
				}
			}

			Q = stlMatrixMultiplication(Q, Hk);

			A = stlMatrixMultiplication(Hk, Ak);
		}

		for (size_t i = 0; i < Q.size(); ++i) {
			for (size_t j = 0; j < Q[0].size(); ++j) {
				if (std::abs(Q[i][j]) < EPSILON) {
					Q[i][j] = 0;
				}
			}
		}

		for (size_t i = 0; i < A.size(); ++i) {
			for (size_t j = 0; j < A[0].size(); ++j) {
				if (std::abs(A[i][j]) < EPSILON) {
					A[i][j] = 0;
				}
			}
		}

		return { A, Q };
	}

	std::vector<std::vector<long double>> getShurMatrix() {
		std::pair<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>> kIteration = hausholderAlgo(matrix);
		std::vector<std::vector<long double>> Ak = stlMatrixMultiplication(kIteration.first, kIteration.second);
		std::vector<long double> currentEigenvalues(Ak.size());
		std::vector<long double> previousEigenvalues(Ak.size());
		for (size_t i = 0; i < Ak.size(); ++i) {
			currentEigenvalues[i] = Ak[i][i];
		}
		for (size_t i = 0; i < 200; ++i) {
			kIteration = hausholderAlgo(Ak);
			Ak = stlMatrixMultiplication(kIteration.first, kIteration.second);
			previousEigenvalues = currentEigenvalues;

			
			for (size_t i = 0; i < Ak.size(); ++i) {
				currentEigenvalues[i] = Ak[i][i];
			}
			bool converged = true;
			for (size_t i = 0; i < currentEigenvalues.size(); ++i) {
				if (std::abs(currentEigenvalues[i] - previousEigenvalues[i]) > EPSILON) {
					converged = false;
					break; 
				}
			}

			if (converged) {
				std::cout << "Converged on " << i + 1 << " iteration." << std::endl;
				break;
			}
				
		}
		return Ak;
	}
private:

	const long double EPSILON = 1e-15;

	std::vector<std::vector<long double>> matrix;
	size_t rows;
	size_t columns;

	long double sign(long double x) {
		if (x > 0) {
			return 1;
		}
		if (x < 0) {
			return -1;
		}
		return 0;
	}

	long double scalarProduct(const std::vector<long double>& a, const std::vector<long double>& b) const{
		long double value = 0;
		for (size_t i = 0; i < a.size(); ++i) {
			value += a[i] * b[i];
		}
		return value;
	}

	std::vector<long double> vectorProjection(const std::vector<long double>& a, const std::vector<long double>& b) const {
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
		if (norm < 1e-5) {
			throw std::runtime_error("dependent vectors");
		}
		std::vector<long double> result(a.size(), 0);
		for (size_t i = 0; i < a.size(); ++i) {
			result[i] = a[i] / norm;
		}
		return result;
	}

	std::vector<long double> hausholderVector(const std::vector<long double>& x) {
		long double euclidianNorm = std::sqrt(scalarProduct(x, x));
		std::vector<long double> result(x.size());

		for (size_t i = 0; i < x.size(); ++i) {
			if (i == 0)
				result[i] = x[i] - sign(x[0]) * euclidianNorm;
			else
				result[i] = x[i];
		}

		long double norm = std::sqrt(scalarProduct(result, result));

		if (norm < EPSILON) {

			throw std::runtime_error("Degenerate vector : Householder transformation not possible.");
		}
		for (size_t i = 0; i < x.size(); ++i) {
			result[i] /= norm;
		}
		return result;
	}

	std::vector<std::vector<long double>> hausholderMatrix(const std::vector<long double>& v, size_t size) {
		std::vector<std::vector<long double>> result(size, std::vector<long double>(size, 0));
		
		for (size_t i = 0; i < size; ++i) {
			result[i][i] = 1;
		}

		for (size_t i = 0; i < v.size(); ++i) {
			for (size_t j = 0; j < v.size(); ++j) {
				result[i][j] -= 2.0 * v[i] * v[j];
			}
		}

		return result;
	}

	std::vector<std::vector<long double>> stlMatrixMultiplication(const std::vector<std::vector<long double>>& a, const std::vector<std::vector<long double>>& b) {
		if (a[0].size() != b.size()) {
			throw std::invalid_argument("wrong matrices sizes");
		}
		std::vector<std::vector<long double>> answer(a.size(), std::vector<long double>(b[0].size()));

		for (size_t i = 0; i < a.size(); ++i) {
			for (size_t j = 0; j < b.size(); ++j) {
				long double element = multiplicationElement(a, b, i, j);
				if (std::abs(element - 0) > EPSILON)
					answer[i][j] = element;
				else
					answer[i][j] = 0;

			}
		}
		return answer;
	}

	long double multiplicationElement(const std::vector<std::vector<long double>>& a, const std::vector<std::vector<long double>>& b, size_t row, size_t col) {
		long double result = 0;
		for (size_t i = 0; i < a[0].size(); ++i) {
			result += a[row][i] * b[i][col];
		}
		return result;
	}
};
