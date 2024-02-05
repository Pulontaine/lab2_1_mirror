#include <iostream>
#include <initializer_list> 
#include <cmath>
#include <iomanip>
#pragma once
using namespace std;

namespace linalg
{
    class Matrix{
    private:
    
        double *m_ptr;
        int m_rows;
        int m_columns;

    public:

        int rows() const noexcept;
        int columns() const noexcept;
        bool empty() const noexcept;
        Matrix& reshape(int rows, int cols);

        Matrix() noexcept;
        Matrix(int rows);
        Matrix(int rows, int cols);
        Matrix(const Matrix &other) noexcept; //копирование
        Matrix(Matrix&& moved) noexcept; //перемещение
        Matrix(initializer_list<double> init_list) noexcept; //initializer_list только с числами внутри
        Matrix(initializer_list<initializer_list<double>> init_list); //initializer_list с initializer_list внутри
        ~Matrix();

        Matrix& operator = (const Matrix& copy) noexcept; //оператор копирования
        Matrix& operator = (Matrix&& moved) noexcept; //оператор перемещения

        double& operator () (int i, int j); //оператор вызова функции 
        double operator () (int i, int j) const; //оператор вызова функции (но матрица константна)
        friend ostream& operator << (ostream& os, const Matrix& out) noexcept; //перегрузка оператора вывода

        Matrix operator + (const Matrix& right) const; // Matrix + Matrix
        Matrix& operator += (const Matrix& right); // Matrix += Matrix
        Matrix operator - (const Matrix& right) const; // Matrix - Matrix
        Matrix& operator -= (const Matrix& right); // Matrix -= Matrix
        Matrix operator * (const Matrix& right) const; // Matrix * Matrix
        friend Matrix operator * (double num, const Matrix& rval); // double * Matrix
        Matrix operator * (const double& num) const; // Matrix * double
        Matrix& operator *= (const Matrix& right); // Matrix *= Matrix
        Matrix& operator *= (double num); // Matrix *= double
        friend bool operator == (const Matrix& left, const Matrix& right) noexcept; // Matrix == Matrix
        friend bool operator != (const Matrix& left, const Matrix& right) noexcept; // Matrix != Matrix

        double norm() const noexcept; // Норма Фробениуса
        double trace() const; // След матрицы
        friend Matrix transpose(const Matrix& m) noexcept; //функция транспонирования 
        friend Matrix concatenate(const Matrix& m1, const Matrix& m2); //функция конкатенирования
        double det() const; //Определитель матрицы
        friend Matrix power(const Matrix &m, int degree); //возведение матрицы в степень

        Matrix& gauss_forward() noexcept; //прямой ход Гаусса
        Matrix& gauss_backward(); //обратный ход Гаусса
        friend Matrix invert(const Matrix &mat); //обратная матрица
        friend Matrix solve(const Matrix &A, const Matrix &f); //решение СЛУ
        int rank() const noexcept; //ранг матрицы

        Matrix operator-() const noexcept; //унарный минус
        double minor(int num_row, int num_col); //минор матрицы
        friend Matrix vertical_concatenate(const Matrix &m1, const Matrix &m2); //вертикальная конкатенация
    };


}; // namespace linalg