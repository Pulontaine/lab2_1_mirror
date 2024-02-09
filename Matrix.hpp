#include "Matrix.h"
#include "test.h"

namespace linalg{

template <typename T>
int Matrix<T>::rows() const noexcept{
        return m_rows;
}

template <typename T>   
int Matrix<T>::columns() const noexcept{
        return m_columns;
}

template <typename T>   
int Matrix<T>::capacity() const noexcept{
        return m_capacity;
}

template <typename T>
T Matrix<T>::ptr(int i) const noexcept{
        return this->m_ptr[i];
}

template <typename T>
bool Matrix<T>::empty() const noexcept{

        return (m_rows == 0 || m_ptr == nullptr || m_columns == 0 || m_capacity == 0);
}

template <typename T>
Matrix<T> Matrix<T>::reshape(int rows, int cols){
        if((*this).empty()){
                cout << "\nEmpty matrix reshape (no shanges)";
                return *this;
        }
        else if(rows <= 0 || cols <= 0 || rows*cols != m_columns*m_rows){
                throw runtime_error("\nERROR: you can't reshape matrix this way.\n");
        }
        else{
                m_columns = cols;
                m_rows = rows;
                m_capacity = cols * rows;
        }
        return *this;
}

template <typename T>
Matrix<T>::Matrix() noexcept{                                       //дефолтный конструктор
        m_columns = 0;
        m_rows = 0;
        m_ptr = nullptr;
        m_capacity = 0;
}

template <typename T>
Matrix<T>::Matrix(int rows): m_rows(rows), m_columns(1), m_capacity(rows){
        if (rows <= 0)
        {
                throw runtime_error("\nError: Size should contain only positive numbers.\n");
        }
        m_ptr = new T[rows];

        // T ptr_in_memory[rows];
        // m_ptr = new(ptr_in_memory) T;
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols) : m_rows(rows), m_columns(cols), m_capacity(rows * cols)  //конструктор с заданным кол-вом строк и/или столбцов
{
        if (rows <= 0 || cols <= 0)
        {
                throw runtime_error("\nError: Size should contain only positive numbers.\n");
        }
        m_ptr = new T[rows * cols];
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &other) noexcept{                                    // конструктор копирования
        m_columns = other.m_columns;
        m_rows = other.m_rows;
        m_capacity = other.capacity();
        if(!other.empty()){
                m_ptr = new T[other.m_rows * other.m_columns];
                for (int i = 0; i < other.m_columns * other.m_rows; i++)
                {
                        m_ptr[i] = (T)other.m_ptr[i];
                }
        }
        else m_ptr = nullptr;
}

template <typename T>
template <typename C>
Matrix<T>::Matrix(const Matrix<C> &other) noexcept{                                    // конструктор копирования
        m_columns = other.columns();
        m_rows = other.rows();
        m_capacity = other.capacity();
        if(!other.empty()){
                m_ptr = new T[other.rows() * other.columns()];
                for (int i = 0; i < other.columns() * other.rows(); i++)
                {
                        m_ptr[i] = (T)other.ptr(i);
                }
        }
        else m_ptr = nullptr;
}



template <typename T>
Matrix<T>::Matrix(Matrix &&moved) noexcept                                  // конструктор перемещения
{ 
        m_rows = moved.rows();
        m_columns = moved.columns();
        m_ptr = moved.m_ptr;
        m_capacity = moved.capacity();
        moved.m_rows = 0;
        moved.m_columns = 0;
        moved.m_ptr = nullptr;
        moved.m_capacity = 0;
}

template <typename T>
template <typename C>
Matrix<T>::Matrix(initializer_list<C> init_list) noexcept{             //конструктор для initializer list
        this->m_rows = init_list.size();
        this->m_columns = 1;
        this->m_capacity = m_rows;
        this->m_ptr = new T[m_rows];

        try{
                m_ptr[0] = static_cast<T>(*(init_list.begin() + 0));
        }
        catch (runtime_error exc) {
		cerr << exc.what();
	}

        for (int i = 0; i < m_rows; i++){
                m_ptr[i] = *(init_list.begin() + i);
        }
}

    template <typename T>
    template <typename U>
    Matrix<T>::Matrix(initializer_list<initializer_list<U>> init_list)
    {
        int tmp = 0;
        auto check = *init_list.begin();
        this->m_rows = init_list.size();
        for (auto i : init_list)
        {
            if (i.size() != check.size())
            {
                throw runtime_error("\nError! Wrong sizes of matrix.\n");
            }
            this->m_columns = i.size();
        }
        this->m_ptr = new T[this->m_rows * this->m_columns];
        auto it = init_list.begin();

        for (int i = 0; i < m_rows * m_columns; i++)
        {
            m_ptr[i] = static_cast<T>(*(it->begin() + i - tmp));
            if ((i + 1) % m_columns == 0)
            {
                tmp += m_columns;
                it++;
            }
        }
        this->m_capacity = this->m_rows * this->m_columns;
    }


template <typename T>
Matrix<T>::~Matrix(){
        delete[] m_ptr;
}

template <typename T>
template <typename C>
Matrix<T> Matrix<T>::operator=(const Matrix<C> &copy) noexcept{                           //копирующий оператор 
        m_rows = copy.rows();
        m_columns = copy.columns();
        m_capacity = copy.capacity();

        if(!this->empty()){
                delete[] m_ptr;
        }
        if(copy.empty()){
                m_ptr = nullptr;
        }
        else{
                m_ptr = new T[copy.rows() * copy.columns()];
                for (int i = 0; i < copy.columns() * copy.rows(); i++)
                {
                        m_ptr[i] = static_cast<T>(copy.ptr(i));
                }
        }
        return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator=(const Matrix<T> &copy) noexcept{                           //копирующий оператор 
        m_rows = copy.m_rows;
        m_columns = copy.m_columns;
        m_capacity = copy.capacity();
        if(!this->empty()){
                delete[] m_ptr;
        }
        if(copy.empty()){
                m_ptr = nullptr;
        }
        else{
                m_ptr = new T[copy.m_rows * copy.m_columns];
                for (int i = 0; i < copy.m_columns * copy.m_rows; i++)
                {
                        m_ptr[i] = static_cast<T>(copy.m_ptr[i]);
                }
        }
        return *this;
}
template <typename T>
Matrix<T> Matrix<T>::operator=(Matrix<T> &&moved) noexcept{                               //перемещающий оператор
        if (&moved != this)
        {
                m_ptr = static_cast<T>(moved.m_ptr);
                m_columns = moved.m_columns;
                m_rows = moved.m_rows;
                m_capacity = moved.m_capacity;
                moved.m_ptr = nullptr;
                moved.m_columns = 0;
                moved.m_rows = 0;
                moved.m_capacity = 0;
        }
        return *this;
}

template <typename T>
T &Matrix<T>::operator () (int i, int j){
        if(this->empty()){
                throw runtime_error("Error: matrix is empty.\n");
        }
        else if(i < 0 || j < 0 || i >= m_rows || j >= m_columns){
                throw runtime_error("Error: index out of range.\n");
        }
        return *(this->m_ptr + i * this->m_columns + j);
}

template <typename T>
T Matrix<T>::operator () (int i, int j) const{
        if(this->empty()){
                throw runtime_error("Error: matrix is empty.\n");
        }
        else if(i < 0 || j < 0 || i >= m_rows || j >= m_columns){
                throw runtime_error("Error: index out of range.\n");
        }
        return *(this->m_ptr + i * this->m_columns + j);
}

template <typename T>
template <typename C>
Matrix<decltype(T{} + C{})> Matrix<T>::operator + (const Matrix<C>& right) const{
        Matrix sum(*this);
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.rows() || m_columns != right.columns()){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }
        for(int i = 0; i < m_columns*m_rows; i++){
                sum.m_ptr[i] += right.ptr(i);
        }
        return sum;
}

template <typename T>
template <typename C>
Matrix<T>& Matrix<T>::operator += (const Matrix<C>& right){
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.rows() || m_columns != right.columns()){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }

        try{
                static_cast<T>(right.ptr(0));
        }
        catch (runtime_error exc) {
		cerr << exc.what();
	}

        for(int i = 0; i < m_columns*m_rows; i++){
                m_ptr[i] += right.ptr(i);
        }
        return *this;
}

template <typename T>
template <typename C>
Matrix<decltype(C{} - T{})> Matrix<T>::operator - (const Matrix<C>& right) const{
        Matrix dif(*this);
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.rows() || m_columns != right.columns()){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }
        for(int i = 0; i < m_columns*m_rows; i++){
                dif.m_ptr[i] -= right.ptr(i);
        }
        return dif;
}

template <typename T>
template <typename C>
Matrix<T>& Matrix<T>::operator -= (const Matrix<C>& right){
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.rows() || m_columns != right.columns()){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }

        try{
                static_cast<T>(right.ptr(0));
        }
        catch (runtime_error exc) {
		cerr << exc.what();
	}

        for(int i = 0; i < m_columns*m_rows; i++){
                m_ptr[i] -= right.ptr(i);
        }
        return *this;
}

template <typename T>
template <typename C>
Matrix<decltype(C{} * T{})> Matrix<T>::operator * (const Matrix<C>& other) const{ //Matrix * Matrix
        if ((this->empty()) || (other.empty())){
		throw std::runtime_error("\nError: empty matrix multiplication.\n");
        }
	else if (m_columns != other.rows()){
		throw std::runtime_error("\nError: can't multiply matrixes.\n");
        }
        Matrix tmp(this->m_rows, other.columns());
        int counter = -1;
        for (int i = 0; i < tmp.m_columns * tmp.m_rows; i++)
        {
            tmp.m_ptr[i] = 0;
            if (!(i % tmp.m_columns))
                counter++;
            for (int j = 0; j < this->m_columns; j++)
            {
                tmp.m_ptr[i] += this->m_ptr[counter * this->m_columns + j] * other.ptr(i % other.columns() + j * other.columns());
            }
        }
        return tmp;

}

template <typename T>
template <typename C>
Matrix<decltype(C{} * double{})> Matrix<T>::operator * (const double& num) const{ // Matrix * double
        if(this->empty()){
                throw runtime_error("\nError: empty matrix.\n");
        }
        Matrix mult(*this);
        for(int i = 0; i < mult.m_columns*mult.m_rows; i++){
                mult.m_ptr[i] = num * mult.m_ptr[i];
        }
        return mult;
}

template <typename T>
template <typename C>
Matrix<T>& Matrix<T>::operator *= (const Matrix<C>& other){ // Matrix *= Matrix
        if(this->empty() || other.empty()){
                throw runtime_error("\nError: matrix is empty.\n");
        }
        else if(this->m_columns != other.m_rows){
                throw runtime_error("\nError: matrix 1 and 2 have different shapes.\n");
        }

        try{
                static_cast<T>(other.ptr(0));
        }
        catch (runtime_error exc) {
		cerr << exc.what();
	}

        Matrix tmp(this->m_rows, other.m_columns);
        int counter = -1;
        for (int i = 0; i < tmp.m_columns * tmp.m_rows; i++){
            tmp.m_ptr[i] = 0;
            if (!(i % tmp.m_columns)) counter++;
            for (int j = 0; j < this->m_columns; j++){
                tmp.m_ptr[i] += this->m_ptr[counter * this->m_columns + j] * other.m_ptr[i % other.m_columns + j * other.m_columns];
            }
        }
        *this = tmp;
        return *this;

}

template <typename T>
Matrix<decltype(double{} * T{})> operator * (double num, const Matrix<T> &rval){ // double * Matrix
        if(rval.empty()){
                throw runtime_error("\nError: matrix is empty.\n");
        }
        Matrix mult(rval);
        for(int i = 0; i < rval.m_columns*rval.m_rows; i++){
                mult.m_ptr[i] = num * rval.m_ptr[i];
        }
        return mult;
}

template <typename T>
Matrix<T>& Matrix<T>::operator *= (double num){ // Matrix *= double
        if(this->empty()){
                throw runtime_error("\nError: empty matrix.\n");
        }
        for(int i = 0; i < this->m_columns*this->m_rows; i++){
                this->m_ptr[i] = this->m_ptr[i] * num;
        }
        return *this;
}

template <typename T>
template <typename C>
bool Matrix<T>::operator == (const Matrix<C>& right) noexcept{ // Matrix == Matrix
        if(this->m_columns != right.columns() || this->m_rows != right.rows()){
                return false;
        }
        
        if(!this->empty() && !right.empty()){
                for(int i = 0; i < this->m_columns * this->m_rows; i++){
                        if(this->m_ptr[i] != right.ptr(i)) return false;
                }
        }
        return true;
}

template <typename T>
template <typename C>
bool Matrix<T>::operator != (const Matrix<C>& right) noexcept{ // Matrix != Matrix
        if(this->m_columns != right.columns() || this->m_rows != right.rows()){
                return true;
        }
        if(!this->empty() && !right.empty()){
                for(int i = 0; i < this->m_columns * this->m_rows; i++){
                        if(this->m_ptr[i] != right.ptr(i)) return true;
                }
        }
        return false;
}

template <typename T>
double Matrix<T>::norm() const noexcept{ //Норма Фробениуса
        double norma = 0;
        for(int i = 0; i < this->m_columns*this->m_rows; i++){
                norma += (this->m_ptr[i]) * (this->m_ptr[i]);
        }
        norma = pow(norma, 0.5);
        return norma;
}

template <typename T>
double Matrix<T>::trace() const{ //След матрицы
        if(this->m_columns != this->m_rows){
                throw runtime_error("\nError: Matrix is not square.\n");
        }
        double traceM = 0;
        for(int i = 0; i < this->m_rows; i++){
                traceM += this->m_ptr[i*this->m_columns + i];
        }
        return traceM;
}

template <typename T>
int Matrix<T>::width(T number) { //Вспомогательная функция для нахождения числа символов
        int width = 0;
        ostringstream stringNumber; //воспользуемся классом ostringstream
        stringNumber << setprecision(10) << number; //в strs теперь double не более 10 знаков
        string str = stringNumber.str(); //double теперь в виде string
        width = (int)str.size(); //получаем кол-во символов в строке
        return width;
}

// template <typename T>
// ostream& operator << (ostream& os, const Matrix<T>& M) noexcept { //перегрузка оператора вывода
//         if (!M.empty()) {
//                 int* max_width;
//                 max_width = new int[M.m_columns]; //в массиве будут самые длинные числа в столбцах
//                 for (int j = 0; j < M.m_columns; j++) { 
//                         max_width[j] = width(M.ptr(j));
//                         for (int i = 1; i < M.m_rows; i++) { //проходим по каждому ряду в столбце
//                                 if (max_width[j] < width(M.ptr(i * M.m_columns + j))) //если меньше макс, макс = width
//                                         max_width[j] = width(M.ptr(i * M.m_columns + j)); //определяем самое "длинное" число в каждом столбце
//                 } }

//                 for (int i = 0; i < M.m_rows; i++) {
//                         os << "|";
//                         for (int j = 0; j < M.m_columns; j++) {
//                                 os << setw(max_width[j] + 1) << setprecision(10) << right << M.m_ptr[i * M.m_columns + j];
//                         }
//                         os << " |\n";
//                 }
//                 delete[] max_width;
//         }
//         else os << "( )\n";
//         return os << endl;
// };

template <typename T>
ostream& operator << (ostream& os, const Matrix<T>& M) noexcept{
        for (int i = 0; i < M.m_rows; i++) {
                        os << "| ";
                        for (int j = 0; j < M.m_columns; j++) {
                                os << M.m_ptr[i * M.m_columns + j];
                                os << " ";
                        }
                        os << "|\n";
                }
        return os << endl;
}

template <typename T>
Matrix<T> transpose(const Matrix<T>& m) noexcept{ //функция транспонирования матрицы
        if(m.empty()) {
                cout << "\nYou transpose empty matrix(no changes).\n";
                return Matrix<T>();
        }
        Matrix<T> N(m);
        for(int i = 0; i < m.m_columns; i++){
                for(int j = 0; j < m.m_rows; j++){
                        N.m_ptr[j*m.m_columns+i] = m.m_ptr[i*m.m_rows+j];
                }
        }
        return N;
}

template <typename T>
Matrix<T> concatenate(const Matrix<T>& m1, const Matrix<T>& m2){ //функция конкатенации
        if(m1.empty() || m2.empty()){
                throw runtime_error("\nError: can't concatenate empty matrix.\n");
        }
        if(m1.m_rows != m2.m_rows){
                throw runtime_error("\nError: can't concatenate matrixes with different quantity of rows.\n");
        }
        int newcolumns = m1.m_columns + m2.m_columns;
        Matrix con(m1.m_rows, newcolumns);
        for(int i = 0; i < m1.m_rows; i++){
                for(int j = 0; j < m1.m_columns; j++){
                        con.m_ptr[i*newcolumns + j] = m1.m_ptr[i*m1.m_columns + j];
                }
                for(int q = 0; q < m2.m_columns; q++){
                        con.m_ptr[i*newcolumns + m1.m_columns + q] = m2.m_ptr[i*m2.m_columns + q];
                }
        }
        return con;
}

template <typename T>
Matrix<T> power(const Matrix<T> &m, int degree){ //Возведение матрицы в степень
        if(m.m_columns != m.m_rows) throw runtime_error("\nError: matrix is not square shape.\n");
        if(m.m_columns == 0 || m.m_rows == 0 || m.m_capacity == 0) throw runtime_error("\nError: matrix is empty.\n");
        if(degree == 0){
                Matrix E(m.m_columns, m.m_columns);
                for(int i = 0; i < m.m_columns; i++){
                        E.m_ptr[i*m.m_columns+i] = 1;
                        for(int j = 0; j < m.m_rows; j++){
                                if (i != j) E.m_ptr[i*m.m_columns+j] = 0;
                        }
                }
                return E;
        }
        Matrix m2(m);
        if(degree < 0) {
                m2 = invert(m2);       
        }
        Matrix positive(m2);
        cout << m2;
        cout << positive;
        for(int i = 1; i < abs(degree); i++){
                m2 = m2 * positive;
        }
        return m2;
}

inline static bool is_equal(double first, double second) { //из-за типа double при маленькой ошибке считаем два числа равными
        return fabs(first - second) <= numeric_limits<double>::epsilon() * 100;
    }

template <typename T>
Matrix<T> &Matrix<T>::gauss_forward() noexcept { //Прямой ход метода Гаусса
        int shift = 0;
        for (int i = 0; i < m_rows - 1; i++) { // m_rows - 1, т. к. не проверяем послелнюю строку
            if (shift + i >= m_columns) break; //shift - смещение для треугольного вида
            if (is_equal(m_ptr[i * m_columns + i + shift], 0)) { //маленькая ошибка из-за типа double 
                for (int p = i + 1; p < m_rows; p++){     //идем по всем следующим за i строкам
                    if (!is_equal(m_ptr[p * m_columns + i + shift], 0)) {
                        for (int g = i + shift; g < m_columns; g++) { //если первый элемент в следующей строке не 0
                            T element = m_ptr[i * m_columns + g]; //то меняем строки местами
                            m_ptr[i * m_columns + g] = m_ptr[p * m_columns + g];
                            m_ptr[p * m_columns + g] = element;
                        }
                        break;
                    }
                }
            }
            if (is_equal(m_ptr[i * m_columns + i + shift], 0)) { // если два нуля в двух смежных строках
                shift++; // то остаемся в той же строке, но смещаемся на один элемент вправо
                i--;
                continue; // не можем выбрать ведущий элемент 0, поэтому переходим в цикле принудительно
            }
            for (int j = i + 1; j < m_rows; j++){
                T coefficient = m_ptr[j * m_columns + i + shift] / m_ptr[i * m_columns + i + shift];
                if (coefficient == 0) continue; //если элемент в строке уже ноль
                for (int k = i; k < m_columns; k++) {
                    m_ptr[j * m_columns + k + shift] -= m_ptr[i * m_columns + k + shift] * coefficient; //хотя бы 1 элемент в строке теперь 0
                }
            }
        }
        return *this;
    }

template <typename T>
Matrix<T> &Matrix<T>::gauss_backward() {
        if ((m_columns - 1 != m_rows) && (this->m_columns / 2) != this->m_rows) 
                throw runtime_error("\nError! Matrix size is wrong or empty\n");
        int shift, real_shift = 0;
        Matrix tmp;
        if (this->m_columns - 1 == this->m_rows){
            shift = 1;
            tmp.m_columns = this->m_columns - 1;
            tmp.m_rows = this->m_rows;
        }
        else{
            shift = this->m_columns / 2;
            tmp.m_columns = this->m_columns / 2;
            tmp.m_rows = this->m_rows;
        }
        tmp.m_ptr = new T[tmp.m_columns * tmp.m_rows];
        for (int i = 0; i < tmp.m_columns * tmp.m_rows; i++){
            if (i % tmp.m_columns == 0 && i != 0)
                real_shift += shift;
            tmp.m_ptr[i] = this->m_ptr[i + real_shift];
        }
        if (!tmp.det()){
            throw runtime_error("\nError! Determinant = 0.\n");
        }

        for (int i = m_rows - 1; i >= 0; i--){
            T norm_coef = this->m_ptr[i * this->m_columns + i];
            this->m_ptr[i * this->m_columns + i] /= norm_coef;
            for (int j = this->m_rows; j < this->m_columns; j++){
                this->m_ptr[i * this->m_columns + j] /= norm_coef;
            }
            if (i != 0){
                for (int j = i - 1; j >= 0; j--){
                    T coef = this->m_ptr[j * this->m_columns + i];
                    this->m_ptr[j * this->m_columns + i] -= this->m_ptr[i * this->m_columns + i] * coef;
                    for (int k = this->m_rows; k < this->m_columns; k++){
                        this->m_ptr[j * this->m_columns + k] -= this->m_ptr[i * this->m_columns + k] * coef;
                    }
                }
            }
        }
        return *this;
    }

template <typename T>
T Matrix<T>::det() const{
        int degree = 0;
        if (this->m_columns != this->m_rows) throw runtime_error("\nError! Matrix is not square.\n");
        if (m_columns == 0 || m_rows == 0) throw runtime_error("\nError! Matrix is empty.\n");
        Matrix tmp(*this);
        for (int i = 0; i < this->m_rows - 1; i++){
            if (is_equal(tmp.m_ptr[i * this->m_columns + i], 0)){
                degree++;
                for (int p = i + 1; p < this->m_rows; p++){
                    if (!is_equal(tmp.m_ptr[p * this->m_columns + i], 0)){
                        for (int g = i; g < tmp.m_columns; g++){
                            T val = tmp.m_ptr[i * this->m_columns + g];
                            tmp.m_ptr[i * this->m_columns + g] = tmp.m_ptr[p * this->m_columns + g];
                            tmp.m_ptr[p * this->m_columns + g] = val;
                        }
                        break;
                    }
                }
            }
            if (is_equal(tmp.m_ptr[i * this->m_columns + i], 0)){
                return 0;
            }
            for (int j = i + 1; j < this->m_rows; j++){
                T coef = tmp.m_ptr[j * this->m_columns + i] / tmp.m_ptr[i * this->m_columns + i];
                for (int k = i; k < m_columns; k++){
                    tmp.m_ptr[j * this->m_columns + k] -= tmp.m_ptr[i * this->m_columns + k] * coef;
                }
            }
        }
        T res = tmp.m_ptr[0];
        for (int i = tmp.m_columns; i < this->m_rows * this->m_columns; i++) {
            if (i % tmp.m_columns == i / tmp.m_columns)
                res *= tmp.m_ptr[i];
        }
        return res * pow((-1), degree);
}

template <typename T>
Matrix<T> invert(const Matrix<T> &mat){
        if (!mat.det()) {
            throw runtime_error("\nError: determinant = 0 (division by zero).\n");
        }
        int shift = mat.m_columns;
        int real_shift = 0;
        Matrix E(mat.m_rows, mat.m_columns);
        for (int i = 0; i < E.m_rows * E.m_columns; i++){
            E.m_ptr[i] = 0;
            if (i % E.m_columns == i / E.m_columns)
                E.m_ptr[i] = 1;
        }
        Matrix tmp = concatenate(mat, E);
        tmp.gauss_forward();
        tmp.gauss_backward();
        Matrix res(tmp.m_rows, mat.m_columns);
        for (int i = 0; i < res.m_columns * res.m_rows; i++){
            if (i % res.m_columns == 0)
                real_shift += shift;
            res.m_ptr[i] = tmp.m_ptr[i + real_shift];
        }

        return res;
}

template <typename T>
Matrix<T> solve(const Matrix<T> &A, const Matrix<T> &f) {
        if (!A.det())
            throw runtime_error("\nError! Determinant of system matrix is 0.\n");
        Matrix m = concatenate(A, f);
        Matrix sol(m.m_rows);
        m.gauss_forward();
        m.gauss_backward();
        for (int i = 0; i < sol.m_rows; i++){
            sol.m_ptr[i] = m.m_ptr[A.m_columns + i * m.m_columns];
        }
        return sol;
}

template <typename T>
int Matrix<T>::rank() const noexcept{ //ранг матрицы
        Matrix m(*this);
        int s = 0;
        int flag = 0;
        m.gauss_forward();
        for(int i = 0; i < m.m_rows; i++){
                flag = 0;
                for(int q = 0; q < m.m_columns; q++){
                        if(m.m_ptr[i*m.m_columns + q] == 0){
                                flag++;
                        }
                }
                if(flag == m.m_columns){
                        s++;
                }
        }
        return m.m_rows - s;
}


//========================= Дополнительные задания

template <typename T>
template <typename C>
Matrix<C> Matrix<T>::operator - () const noexcept{ //унарный минус  
    Matrix result(m_rows, m_columns);
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_columns; j++) {
            result(i, j) = -(*this)(i, j);
        }
    }
    return result;
}

template <typename T>
T Matrix<T>::minor(int num_row, int num_col){ //минор матрицы (индексация с 0)
        if (num_row < 0 || num_col < 0)
                throw runtime_error("\nError! Can't use negative numbers.\n");
        if (num_row > m_rows - 1 || num_col > m_columns - 1)
                throw runtime_error("\nError! Out of matrix.\n");
        Matrix tmp(m_rows - 1, m_columns - 1);
        int s = 0;
        for(int i = 0; i < m_columns*m_rows; i++){
                if(i/m_columns == num_row) continue;
                if(i/m_rows == num_col) continue;
                tmp.m_ptr[s] = this->m_ptr[i];
                s++;
        }
        return tmp.det();
    }

template <typename T>
Matrix<T> vertical_concatenate(const Matrix<T> &m1, const Matrix<T> &m2){ //вертикальная конкатенация
        if(m1.empty() || m2.empty()){
                throw runtime_error("\nError: can't concatenate empty matrix.\n");
        }
        if(m1.m_rows != m2.m_rows){
                throw runtime_error("\nError: can't concatenate matrixes with different quantity of columns.\n");
        }
        Matrix tmp(m1.m_rows + m2.m_rows, m1.m_columns);
        for (int i = 0; i < tmp.m_rows * tmp.m_columns; i++){
            if ((i < (m1.m_columns * m1.m_rows))) tmp.m_ptr[i] = m1.m_ptr[i];
            else tmp.m_ptr[i] = m2.m_ptr[i - m1.m_columns * m1.m_rows];
        }
        return tmp;
    }

template <typename T>
void Matrix<T>::reserve(int n) const{
        this->m_capacity = n;
        this->m_ptr = (T*)realloc(this->m_ptr, n);
}

template <typename T>
void Matrix<T>::shrink_to_fit() const{
        this->m_ptr = (T*)realloc(this->m_ptr, this->m_rows * this-> m_columns);
        this->m_capacity = this->m_rows * this-> m_columns;
}

template <typename T>
void Matrix<T>::clear() const{
        this->m_columns = 0;
        this->m_rows = 0;
}


}