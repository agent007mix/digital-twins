#pragma once
#include <filesystem>

namespace matrix {

	class SparseMatrixCSR {
		size_t row_number_;           // Число строк
		size_t column_number_;        // Число столбцов
		size_t non_zero_elements_;    // Количество ненулевых элементов

		double* values_;              // Значения ненулевых элементов
		size_t* col_indices_;         // Индексы столбцов для каждого значения
		size_t* row_ptr_;             // Указатели на начало строк в values_

	public:
		SparseMatrixCSR(); // Конструктор по умолчанию
		SparseMatrixCSR(const std::filesystem::path& matrix_path); // Конструктор из файла
		SparseMatrixCSR(size_t rows, size_t cols, size_t non_zeros, double* values, size_t* col_indices, size_t* row_ptr); // Конструктор с параметрами
		SparseMatrixCSR(const SparseMatrixCSR& other);
		SparseMatrixCSR& operator=(const SparseMatrixCSR& other);
		SparseMatrixCSR(SparseMatrixCSR&& other) noexcept;
		SparseMatrixCSR& operator=(SparseMatrixCSR&& other) noexcept;
		~SparseMatrixCSR(); // Деструктор

		void WriteToFile(const std::filesystem::path& file_path) const; // Вывод матрицы в файл
		SparseMatrixCSR Slice(size_t row_start, size_t row_end, size_t col_start, size_t col_end) const; // Извлечение подматрицы

		double get_element(size_t row, size_t col) const; // Получение элемента по индексу

		// Добавляем объявления операторов
		SparseMatrixCSR operator+(const SparseMatrixCSR& other) const; // Бинарное сложение матриц
		SparseMatrixCSR operator-(const SparseMatrixCSR& other) const; // Бинарное вычитание матриц
		SparseMatrixCSR operator-() const;                             // Унарный минус
		double* operator*(const double* vector) const;                 // Умножение на плотный вектор
		SparseMatrixCSR operator*(double scalar) const;                // Умножение на скаляр
		double operator[](std::tuple<size_t, size_t> ind) const;          // Индексация по мультииндексу
		double& operator[](std::tuple<size_t, size_t> ind);               // Индексация по мультииндексу
	};
}
