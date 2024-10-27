#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstring> 
#include "SparseMatrixCSR.h"

namespace matrix {

	// ����������� �� ���������
	SparseMatrixCSR::SparseMatrixCSR() : row_number_(0), column_number_(0), non_zero_elements_(0), values_(nullptr), col_indices_(nullptr), row_ptr_(nullptr) {}

	// ����������� ��� ������ ������� �� ����� � ������� MTX
	SparseMatrixCSR::SparseMatrixCSR(const std::filesystem::path& matrix_path) {
		std::ifstream fin(matrix_path);
		if (!fin) {
			throw std::runtime_error("Error: Unable to open file");
		}

		// ������ ����� �����, �������� � ��������� ���������
		fin >> row_number_ >> column_number_ >> non_zero_elements_;

		// �������� ������ ��� ��������, �������� �������� � ���������� �� ������
		values_ = new double[non_zero_elements_];
		col_indices_ = new size_t[non_zero_elements_];
		row_ptr_ = new size_t[row_number_ + 1]{ 0 }; // �������������� ������ row_ptr ������

		// ��������� ������ ��� �������� ��������� ��������� � ������ ������
		size_t* row_counts = new size_t[row_number_]{ 0 };

		size_t row, col;
		double value;
		for (size_t i = 0; i < non_zero_elements_; ++i) {
			fin >> row >> col >> value;
			--row; --col;
			values_[i] = value;
			col_indices_[i] = col;
			row_counts[row]++;
		}

		// ��������� row_ptr (������� ������ ������ ������ � values_)
		row_ptr_[0] = 0;
		for (size_t i = 1; i <= row_number_; ++i) {
			row_ptr_[i] = row_ptr_[i - 1] + row_counts[i - 1];
		}

		delete[] row_counts; // ����������� ��������� ������
		fin.close();
	}

	// ����� ����������� � �����������
	SparseMatrixCSR::SparseMatrixCSR(size_t rows, size_t cols, size_t non_zeros, double* values, size_t* col_indices, size_t* row_ptr)
		: row_number_(rows), column_number_(cols), non_zero_elements_(non_zeros) {
		values_ = new double[non_zeros];
		col_indices_ = new size_t[non_zeros];
		row_ptr_ = new size_t[rows + 1];

		std::memcpy(values_, values, non_zeros * sizeof(double));
		std::memcpy(col_indices_, col_indices, non_zeros * sizeof(size_t));
		std::memcpy(row_ptr_, row_ptr, (rows + 1) * sizeof(size_t));
	}

	// ����������� �����������
	SparseMatrixCSR::SparseMatrixCSR(const SparseMatrixCSR& other) :
		row_number_(other.row_number_),
		column_number_(other.column_number_),
		non_zero_elements_(other.non_zero_elements_) {

		values_ = new double[non_zero_elements_];
		col_indices_ = new size_t[non_zero_elements_];
		row_ptr_ = new size_t[row_number_ + 1];

		std::memcpy(values_, other.values_, non_zero_elements_ * sizeof(double));
		std::memcpy(col_indices_, other.col_indices_, non_zero_elements_ * sizeof(size_t));
		std::memcpy(row_ptr_, other.row_ptr_, (row_number_ + 1) * sizeof(size_t));
	}

	// �������� ������������ ������������
	SparseMatrixCSR& SparseMatrixCSR::operator=(const SparseMatrixCSR& other) {
		if (this != &other) {
			// ����������� ������ ������
			delete[] values_;
			delete[] col_indices_;
			delete[] row_ptr_;

			// �������� ����� ��������
			row_number_ = other.row_number_;
			column_number_ = other.column_number_;
			non_zero_elements_ = other.non_zero_elements_;

			values_ = new double[non_zero_elements_];
			col_indices_ = new size_t[non_zero_elements_];
			row_ptr_ = new size_t[row_number_ + 1];

			std::memcpy(values_, other.values_, non_zero_elements_ * sizeof(double));
			std::memcpy(col_indices_, other.col_indices_, non_zero_elements_ * sizeof(size_t));
			std::memcpy(row_ptr_, other.row_ptr_, (row_number_ + 1) * sizeof(size_t));
		}
		return *this;
	}

	// ����������� �����������
	SparseMatrixCSR::SparseMatrixCSR(SparseMatrixCSR&& other) noexcept :
		row_number_(other.row_number_),
		column_number_(other.column_number_),
		non_zero_elements_(other.non_zero_elements_),
		values_(other.values_),
		col_indices_(other.col_indices_),
		row_ptr_(other.row_ptr_) {
		// �������� ������ � ���������
		other.values_ = nullptr;
		other.col_indices_ = nullptr;
		other.row_ptr_ = nullptr;
	}

	// �������� ������������ ������������
	SparseMatrixCSR& SparseMatrixCSR::operator=(SparseMatrixCSR&& other) noexcept {
		if (this != &other) {
			// ����������� ������ ������
			delete[] values_;
			delete[] col_indices_;
			delete[] row_ptr_;

			// ���������� ����� ������
			row_number_ = other.row_number_;
			column_number_ = other.column_number_;
			non_zero_elements_ = other.non_zero_elements_;
			values_ = other.values_;
			col_indices_ = other.col_indices_;
			row_ptr_ = other.row_ptr_;

			// �������� ������ � ���������
			other.values_ = nullptr;
			other.col_indices_ = nullptr;
			other.row_ptr_ = nullptr;
		}
		return *this;
	}

	// ����������
	SparseMatrixCSR::~SparseMatrixCSR() {
		delete[] values_;
		delete[] col_indices_;
		delete[] row_ptr_;
	}

	// ����� ��� ������ ������� � ���� � ������� MTX
	void SparseMatrixCSR::WriteToFile(const std::filesystem::path& file_path) const {
		std::ofstream fout(file_path);
		if (!fout) {
			throw std::runtime_error("Error: Unable to open file");
		}

		fout << row_number_ << " " << column_number_ << " " << non_zero_elements_ << "\n";

		for (size_t i = 0; i < row_number_; ++i) {
			for (size_t j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
				fout << i + 1 << " " << col_indices_[j] + 1 << " " << values_[j] << "\n"; 
			}
		}

		fout.close();
	}

	// ����� ��� ���������� ����������
	SparseMatrixCSR SparseMatrixCSR::Slice(size_t row_start, size_t row_end, size_t col_start, size_t col_end) const {
		// �������� ������
		if (row_start >= row_end || col_start >= col_end || row_start >= row_number_ || row_end > row_number_ || col_start >= column_number_ || col_end > column_number_) {
			throw std::out_of_range("Slice indices are out of range");
		}

		size_t slice_non_zero_elements = 0;
		// ������� ��������� ��������� � ����������
		for (size_t i = row_start; i < row_end; ++i) {
			for (size_t j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
				if (col_indices_[j] >= col_start && col_indices_[j] < col_end) {
					slice_non_zero_elements++;
				}
			}
		}

		// �������� �������� ��� ����������
		double* slice_values = new double[slice_non_zero_elements];
		size_t* slice_col_indices = new size_t[slice_non_zero_elements];
		size_t* slice_row_ptr = new size_t[row_end - row_start + 1]{ 0 };

		size_t slice_index = 0;
		for (size_t i = row_start; i < row_end; ++i) {
			for (size_t j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
				if (col_indices_[j] >= col_start && col_indices_[j] < col_end) {
					slice_values[slice_index] = values_[j];
					slice_col_indices[slice_index] = col_indices_[j] - col_start; // ������������ ������ �������
					slice_row_ptr[i - row_start + 1]++;
					slice_index++;
				}
			}
		}

		// ��������� ���������� ������� � slice_row_ptr
		for (size_t i = 1; i <= row_end - row_start; ++i) {
			slice_row_ptr[i] += slice_row_ptr[i - 1];
		}

		return SparseMatrixCSR(row_end - row_start, col_end - col_start, slice_non_zero_elements, slice_values, slice_col_indices, slice_row_ptr);
	}

	// ��������� �������� �� �������
	double SparseMatrixCSR::get_element(size_t row, size_t col) const {
		for (size_t j = row_ptr_[row]; j < row_ptr_[row + 1]; ++j) {
			if (col_indices_[j] == col) {
				return values_[j];
			}
		}
		return 0.0; 
	}

		

	SparseMatrixCSR SparseMatrixCSR::operator+(const SparseMatrixCSR& other) const {
		if (row_number_ != other.row_number_ || column_number_ != other.column_number_) {
			throw std::invalid_argument("Matrix dimensions must match for addition");
		}

		size_t new_non_zero_elements = 0;
		size_t* new_col_indices = new size_t[non_zero_elements_ + other.non_zero_elements_];
		double* new_values = new double[non_zero_elements_ + other.non_zero_elements_];
		size_t* new_row_ptr = new size_t[row_number_ + 1]{ 0 };

		for (size_t i = 0; i < row_number_; ++i) {
			size_t this_ptr = row_ptr_[i];
			size_t other_ptr = other.row_ptr_[i];
			new_row_ptr[i] = new_non_zero_elements;

			while (this_ptr < row_ptr_[i + 1] || other_ptr < other.row_ptr_[i + 1]) {
				if (this_ptr < row_ptr_[i + 1] && (other_ptr >= other.row_ptr_[i + 1] || col_indices_[this_ptr] < other.col_indices_[other_ptr])) {
					new_values[new_non_zero_elements] = values_[this_ptr];
					new_col_indices[new_non_zero_elements] = col_indices_[this_ptr];
					++new_non_zero_elements;
					++this_ptr;
				}
				else if (other_ptr < other.row_ptr_[i + 1] && (this_ptr >= row_ptr_[i + 1] || other.col_indices_[other_ptr] < col_indices_[this_ptr])) {
					new_values[new_non_zero_elements] = other.values_[other_ptr];
					new_col_indices[new_non_zero_elements] = other.col_indices_[other_ptr];
					++new_non_zero_elements;
					++other_ptr;
				}
				else {
					double sum = values_[this_ptr] + other.values_[other_ptr];
					if (sum != 0) { // �������� �� ��������� ��������
						new_values[new_non_zero_elements] = sum;
						new_col_indices[new_non_zero_elements] = col_indices_[this_ptr];
						++new_non_zero_elements;
					}
					++this_ptr;
					++other_ptr;
				}
			}
		}
		new_row_ptr[row_number_] = new_non_zero_elements;

		return SparseMatrixCSR(row_number_, column_number_, new_non_zero_elements, new_values, new_col_indices, new_row_ptr);
	}


	// ���������� ��������� ��������� ��� ����������� ������
	SparseMatrixCSR SparseMatrixCSR::operator-(const SparseMatrixCSR& other) const {
		if (row_number_ != other.row_number_ || column_number_ != other.column_number_) {
			throw std::invalid_argument("Matrix dimensions must match for subtraction");
		}

		size_t new_non_zero_elements = 0;
		size_t* new_col_indices = new size_t[non_zero_elements_ + other.non_zero_elements_];
		double* new_values = new double[non_zero_elements_ + other.non_zero_elements_];
		size_t* new_row_ptr = new size_t[row_number_ + 1]{ 0 };

		for (size_t i = 0; i < row_number_; ++i) {
			size_t this_ptr = row_ptr_[i];
			size_t other_ptr = other.row_ptr_[i];
			new_row_ptr[i] = new_non_zero_elements;

			while (this_ptr < row_ptr_[i + 1] || other_ptr < other.row_ptr_[i + 1]) {
				if (this_ptr < row_ptr_[i + 1] && (other_ptr >= other.row_ptr_[i + 1] || col_indices_[this_ptr] < other.col_indices_[other_ptr])) {
					new_values[new_non_zero_elements] = values_[this_ptr];
					new_col_indices[new_non_zero_elements] = col_indices_[this_ptr];
					++new_non_zero_elements;
					++this_ptr;
				}
				else if (other_ptr < other.row_ptr_[i + 1] && (this_ptr >= row_ptr_[i + 1] || other.col_indices_[other_ptr] < col_indices_[this_ptr])) {
					new_values[new_non_zero_elements] = -other.values_[other_ptr];
					new_col_indices[new_non_zero_elements] = other.col_indices_[other_ptr];
					++new_non_zero_elements;
					++other_ptr;
				}
				else {
					double difference = values_[this_ptr] - other.values_[other_ptr];
					if (difference != 0) { // �������� �� ��������� ��������
						new_values[new_non_zero_elements] = difference;
						new_col_indices[new_non_zero_elements] = col_indices_[this_ptr];
						++new_non_zero_elements;
					}
					++this_ptr;
					++other_ptr;
				}
			}
		}
		new_row_ptr[row_number_] = new_non_zero_elements;

		return SparseMatrixCSR(row_number_, column_number_, new_non_zero_elements, new_values, new_col_indices, new_row_ptr);
	}

		// ������� �����
		SparseMatrixCSR SparseMatrixCSR::operator-() const {
			double* new_values = new double[non_zero_elements_];
			for (size_t i = 0; i < non_zero_elements_; ++i) {
				new_values[i] = -values_[i];
			}
			return SparseMatrixCSR(row_number_, column_number_, non_zero_elements_, new_values, col_indices_, row_ptr_);
		}

		// ��������� �� ������� ������
		double* SparseMatrixCSR::operator*(const double* vector) const {
			double* result = new double[row_number_] { 0.0 };
			for (size_t i = 0; i < row_number_; ++i) {
				for (size_t j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j) {
					result[i] += values_[j] * vector[col_indices_[j]];
				}
			}
			return result;
		}

		// ��������� �� ������
		SparseMatrixCSR SparseMatrixCSR::operator*(double scalar) const {
			double* new_values = new double[non_zero_elements_];
			for (size_t i = 0; i < non_zero_elements_; ++i) {
				new_values[i] = values_[i] * scalar;
			}
			return SparseMatrixCSR(row_number_, column_number_, non_zero_elements_, new_values, col_indices_, row_ptr_);
		}

		
		// ���������� ���������� ����������
		double SparseMatrixCSR::operator[](std::tuple<size_t, size_t> ind) const {
			size_t row = std::get<0>(ind);
			size_t col = std::get<1>(ind);

			// �������� ������
			if (row >= row_number_ || col >= column_number_) {
				throw std::out_of_range("Index is out of range");
			}

			return get_element(row, col);
		}

		double& SparseMatrixCSR::operator[](std::tuple<size_t, size_t> ind) {
			size_t row = std::get<0>(ind);
			size_t col = std::get<1>(ind);

			// �������� ������
			if (row >= row_number_ || col >= column_number_) {
				throw std::out_of_range("Index is out of range");
			}

			// ������� ������ �������� � CSR
			for (size_t j = row_ptr_[row]; j < row_ptr_[row + 1]; ++j) {
				if (col_indices_[j] == col) {
					return values_[j]; // ���������� ������ �� ������������ �������
				}
			}

			// ����������, ���� ���������� � �������� ��������, �. � � ���������������
			throw std::runtime_error("Attempt to modify a zero element. Use the addition method instead.");
		}
}

