#include <iostream>
#include <initializer_list> 
#include <cmath>
#include <iomanip>
#pragma once
using namespace std;

namespace linalg
{
    template <class T>
    class Matrix{
    private:
    
        T *m_ptr;
        int m_rows;
        int m_columns;
        int m_capacity;

    public:

        int rows() const noexcept;
        int columns() const noexcept;
        bool empty() const noexcept;
        Matrix<T> reshape(int rows, int cols);

        Matrix() noexcept;
        Matrix(int rows);
        Matrix(int rows, int cols);

        template <typename C>
        Matrix(const Matrix<C> &other) noexcept; //копирование

        Matrix(const Matrix<T> &other) noexcept; //копирование

        Matrix(Matrix<T> &&moved) noexcept; //перемещение
        ~Matrix();

        template <typename C>
        Matrix(initializer_list<C> init_list) noexcept; //initializer_list только с числами внутри

        template <typename C>
        Matrix(initializer_list<initializer_list<C>> init_list); //initializer_list с initializer_list внутри

        template <typename C>
        Matrix<T> operator = (const Matrix<C>& copy) noexcept; //оператор копирования

        Matrix<T> operator = (const Matrix<T>& copy) noexcept; //оператор копирования

        Matrix<T> operator = (Matrix<T>&& moved) noexcept; //оператор перемещения

        T& operator () (int i, int j); //оператор вызова функции 
        T operator () (int i, int j) const; //оператор вызова функции (но матрица константна)

        template <typename C>
        friend ostream& operator << (ostream& os, const Matrix<C>& out) noexcept; //перегрузка оператора вывода

        template <typename C>
        Matrix<decltype(T{} + C{})> operator + (const Matrix<C>& right) const; // Matrix + Matrix

        template <typename C>
        Matrix<T>& operator += (const Matrix<C>& right); // Matrix += Matrix

        template <typename C>
        Matrix<decltype(C{} - T{})> operator - (const Matrix<C>& right) const; // Matrix - Matrix

        template <typename C>
        Matrix<T>& operator -= (const Matrix<C>& right); // Matrix -= Matrix

        template <typename C>
        Matrix<decltype(C{} * T{})> operator * (const Matrix<C>& right) const; // Matrix * Matrix

        template <typename C>
        friend Matrix<decltype(double{} * C{})> operator * (double num, const Matrix<C>& rval); // double * Matrix

        template <typename C>
        Matrix<decltype(C{} * double{})> operator * (const double& num) const; // Matrix * double

        template <typename C>
        Matrix<T>& operator *= (const Matrix<C>& right); // Matrix *= Matrix

        Matrix<T>& operator *= (double num); // Matrix *= double

        template <typename C>
        bool operator == (const Matrix<C>& right) noexcept; // Matrix == Matrix

        template <typename C>
        bool operator != (const Matrix<C>& right) noexcept; // Matrix != Matrix
        
        double norm() const noexcept; // Норма Фробениуса
        double trace() const; // След матрицы

        template <typename C>
        friend Matrix<C> transpose(const Matrix<C>& m) noexcept; //функция транспонирования 

        template <typename C>
        friend Matrix<C> concatenate(const Matrix<C>& m1, const Matrix<C>& m2); //функция конкатенирования

        T det() const; //Определитель матрицы

        template <typename C>
        friend Matrix<C> power(const Matrix<C> &m, int degree); //возведение матрицы в степень

        Matrix<T>& gauss_forward() noexcept; //прямой ход Гаусса
        Matrix<T>& gauss_backward(); //обратный ход Гаусса

        template <typename C>
        friend Matrix<C> invert(const Matrix<C> &mat); //обратная матрица

        template <typename C>
        friend Matrix<C> solve(const Matrix<C> &A, const Matrix<C> &f); //решение СЛУ
        
        int rank() const noexcept; //ранг матрицы

        template <typename C>
        Matrix<C> operator-() const noexcept; //унарный минус

        T minor(int num_row, int num_col); //минор матрицы

        template <typename C>
        friend Matrix<C> vertical_concatenate(const Matrix<C> &m1, const Matrix<C> &m2); //вертикальная конкатенация

        int capacity() const;

        void reserve(int n) const;

        void shrink_to_fit() const;

        void clear() const;
    };


}; // namespace linalg