#pragma once
#include <filesystem>

namespace matrix {

	class SparseMatrixCSR {
		size_t row_number_;           // ����� �����
		size_t column_number_;        // ����� ��������
		size_t non_zero_elements_;    // ���������� ��������� ���������

		double* values_;              // �������� ��������� ���������
		size_t* col_indices_;         // ������� �������� ��� ������� ��������
		size_t* row_ptr_;             // ��������� �� ������ ����� � values_

	public:
		SparseMatrixCSR(); // ����������� �� ���������
		SparseMatrixCSR(const std::filesystem::path& matrix_path); // ����������� �� �����
		SparseMatrixCSR(size_t rows, size_t cols, size_t non_zeros, double* values, size_t* col_indices, size_t* row_ptr); // ����������� � �����������
		SparseMatrixCSR(const SparseMatrixCSR& other);
		SparseMatrixCSR& operator=(const SparseMatrixCSR& other);
		SparseMatrixCSR(SparseMatrixCSR&& other) noexcept;
		SparseMatrixCSR& operator=(SparseMatrixCSR&& other) noexcept;
		~SparseMatrixCSR(); // ����������

		void WriteToFile(const std::filesystem::path& file_path) const; // ����� ������� � ����
		SparseMatrixCSR Slice(size_t row_start, size_t row_end, size_t col_start, size_t col_end) const; // ���������� ����������

		double get_element(size_t row, size_t col) const; // ��������� �������� �� �������

		// ��������� ���������� ����������
		SparseMatrixCSR operator+(const SparseMatrixCSR& other) const; // �������� �������� ������
		SparseMatrixCSR operator-(const SparseMatrixCSR& other) const; // �������� ��������� ������
		SparseMatrixCSR operator-() const;                             // ������� �����
		double* operator*(const double* vector) const;                 // ��������� �� ������� ������
		SparseMatrixCSR operator*(double scalar) const;                // ��������� �� ������
		double operator[](std::tuple<size_t, size_t> ind) const;          // ���������� �� �������������
		double& operator[](std::tuple<size_t, size_t> ind);               // ���������� �� �������������
	};
}
