#include <iostream>
#include <tuple>
#include <fstream>
#include "SparseMatrixCSR.h"

int main() {
    std::filesystem::path input_file = "input.txt";
    std::filesystem::path input_file2 = "input2.txt";
    std::filesystem::path output_file = "output_matrix.txt";

    try {
        // Создание экземпляра разреженной матрицы из файла
        matrix::SparseMatrixCSR sparse_matrix(input_file);
        matrix::SparseMatrixCSR sparse(input_file2);
        // Вывод исходной матрицы в файл
        sparse_matrix.WriteToFile(output_file);
        std::cout << "Sparse matrix successfully read from " << input_file << " and written to " << output_file << ".\n";

        // Тест сложения матриц
        matrix::SparseMatrixCSR sum_matrix = sparse_matrix+sparse;
        std::filesystem::path sum_output_file = "sum_matrix.txt";
        sum_matrix.WriteToFile(sum_output_file);
        std::cout << "Matrix addition successful. Result written to " << sum_output_file << ".\n";

        // Тест вычитания матриц
        matrix::SparseMatrixCSR diff_matrix = sparse_matrix - sparse;
        std::filesystem::path diff_output_file = "diff_matrix.txt";
        diff_matrix.WriteToFile(diff_output_file);
        std::cout << "Matrix subtraction successful. Result written to " << diff_output_file << ".\n";

        // Тест унарного минуса
        matrix::SparseMatrixCSR neg_matrix = -sparse_matrix;
        std::filesystem::path neg_output_file = "neg_matrix.txt";
        neg_matrix.WriteToFile(neg_output_file);
        std::cout << "Unary minus successful. Result written to " << neg_output_file << ".\n";

        // Тест умножения на плотный вектор
        size_t vector_size = 5;
        double* dense_vector = new double[vector_size];
        for (size_t i = 0; i < vector_size; ++i) {
            dense_vector[i] = 2.0;  // Пример вектора
        }
        double* result_vector = sparse * dense_vector;

        std::filesystem::path vector_output_file = "vector_result.txt";
        std::ofstream vector_out(vector_output_file);
        if (!vector_out) {
            throw std::runtime_error("Error: Unable to open file for vector result");
        }
        for (size_t i = 0; i < vector_size; ++i) {
            vector_out << result_vector[i] << "\n";
        }
        vector_out.close();
        std::cout << "Matrix-vector multiplication successful. Result written to " << vector_output_file << ".\n";
        delete[] dense_vector;
        delete[] result_vector;

        // Тест умножения на скаляр
        
        double scalar = 2.0;
        matrix::SparseMatrixCSR scaled_matrix = sparse_matrix * scalar;
        std::filesystem::path scaled_output_file = "scaled_matrix.txt";
        scaled_matrix.WriteToFile(scaled_output_file);
        std::cout << "Matrix-scalar multiplication successful. Result written to " << scaled_output_file << ".\n";

        // Тест индексации по мультииндексу
        
        double value = sparse_matrix[{4, 4}];
        sparse_matrix[{4, 4}] = 5.5;
        std::cout << "Element : " << value  << "\n";
        value = sparse_matrix[{4, 4}];
        std::cout << "Element new: " << value << "\n";
        
        // Извлечение подматрицы для проверки Slice
        size_t row_start = 3, row_end = 5;
        size_t col_start = 3, col_end = 5;
        matrix::SparseMatrixCSR submatrix = sparse_matrix.Slice(row_start, row_end, col_start, col_end);

        std::filesystem::path submatrix_output_file = "submatrix_output.txt";
        submatrix.WriteToFile(submatrix_output_file);
        std::cout << "Submatrix successfully written to " << submatrix_output_file << ".\n";

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
