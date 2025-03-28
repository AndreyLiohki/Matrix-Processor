#pragma once
#include <vector>
#include <iostream>
#include "Matrix.h"

class MatrixOperations {
public:
	static Matrix add(const Matrix& a, const Matrix& b) {
		if (a.getNumberOfRows() != b.getNumberOfRows() || a.getNumberOfColumns() != b.getNumberOfColumns()) {
			throw std::invalid_argument("Both matrices must have same size.");
		}

		Matrix additionResult(a.getNumberOfRows(), b.getNumberOfColumns());

		for (size_t i = 0; i < a.getNumberOfRows(); ++i) {
			for (size_t j = 0; j < a.getNumberOfColumns(); ++j) {
				additionResult.set(i, j, a.get(i, j) + b.get(i, j));
			}
		}
		return additionResult;

	}

	static Matrix substraction(const Matrix& a, const Matrix& b) {
		if (a.getNumberOfRows() != b.getNumberOfRows() || a.getNumberOfColumns() != b.getNumberOfColumns()) {
			throw std::invalid_argument("Both matrices must have same size.");
		}

		Matrix substractionResult(a.getNumberOfRows(), b.getNumberOfColumns());

		for (size_t i = 0; i < a.getNumberOfRows(); ++i) {
			for (size_t j = 0; j < a.getNumberOfColumns(); ++j) {
				substractionResult.set(i, j, a.get(i, j) - b.get(i, j));
			}
		}
		return substractionResult;

	}

	static Matrix multiplication(const Matrix& a, const Matrix& b) {
		if (a.getNumberOfColumns() != b.getNumberOfRows()) {
			throw std::invalid_argument("wrong matrices sizes");
		}
		Matrix multiplicationResult(a.getNumberOfRows(), b.getNumberOfColumns());

		for (size_t i = 0; i < a.getNumberOfRows(); ++i) {
			for (size_t j = 0; j < b.getNumberOfColumns(); ++j) {
				long double element = multiplicationElement(a, b, i, j);
				multiplicationResult.set(i, j, element);

			}
		}
		return multiplicationResult;
	}
private:
	static long double multiplicationElement(const Matrix& a, const Matrix& b, const size_t row, const size_t col) {
		long double result = 0;
		for (size_t i = 0; i < a.getNumberOfColumns(); ++i) {
			result += a.get(row, i) * b.get(i, col);
		}
		return result;
	}
};
