#include "Matrix.h"
#include "test.h"

namespace linalg{

int Matrix::rows() const noexcept{
        return m_rows;
}

int Matrix::columns() const noexcept{
        return m_columns;
}

bool Matrix::empty() const noexcept{

        return (m_rows == 0 || m_ptr == nullptr || m_columns == 0);
}

Matrix& Matrix::reshape(int rows, int cols){
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
        }
        return *this;
}

Matrix::Matrix() noexcept{                                       //дефолтный конструктор
        m_columns = 0;
        m_rows = 0;
        m_ptr = nullptr;
}

Matrix::Matrix(int rows): m_rows(rows), m_columns(1){
        if (rows <= 0)
        {
                throw runtime_error("\nError: Size should contain only positive numbers.\n");
        }
        m_ptr = new double[rows];
}

Matrix::Matrix(int rows, int cols) : m_rows(rows), m_columns(cols)  //конструктор с заданным кол-вом строк и/или столбцов
{
        if (rows <= 0 || cols <= 0)
        {
                throw runtime_error("\nError: Size should contain only positive numbers.\n");
        }
        m_ptr = new double[rows * cols];
}

Matrix::Matrix(const Matrix &other) noexcept{                                    // конструктор копирования
        m_columns = other.m_columns;
        m_rows = other.m_rows;
        if(!other.empty()){
                m_ptr = new double[other.m_rows * other.m_columns];
                for (int i = 0; i < other.m_columns * other.m_rows; i++)
                {
                        m_ptr[i] = other.m_ptr[i];
                }
        }
        else m_ptr = nullptr;
}

Matrix::Matrix(Matrix &&moved) noexcept                                  // конструктор перемещения
{ 
        m_rows = moved.rows();
        m_columns = moved.columns();
        m_ptr = moved.m_ptr;
        moved.m_rows = 0;
        moved.m_columns = 0;
        moved.m_ptr = nullptr;
}

Matrix::Matrix(initializer_list<double> init_list) noexcept{             //конструктор для initializer list
        this->m_rows = init_list.size();
        this->m_columns = 1;
        this->m_ptr = new double[m_rows];
        for (int i = 0; i < m_rows; i++){
                m_ptr[i] = *(init_list.begin() + i);
        }
}

Matrix::Matrix(initializer_list<initializer_list<double>> init_list){    //конструктор для initializer list с initializer list внутри
        int tmp = 0;
        auto check = *init_list.begin();
        this->m_rows = init_list.size();
        for (auto i : init_list){
                if (i.size() != check.size()){
                        throw runtime_error("\nError! Wrong sizes of matrix.\n");
                }
                this->m_columns = i.size();
        }
        this->m_ptr = new double[this->m_rows * this->m_columns];
        auto it = init_list.begin();

        for (int i = 0; i < m_rows * m_columns; i++){
                m_ptr[i] = *(it->begin() + i - tmp);
                if ((i + 1) % m_columns == 0){
                        tmp += m_columns;
                        it++;
                }
        }
}

Matrix::~Matrix(){
        delete[] m_ptr;
}

Matrix &Matrix::operator=(const Matrix &copy) noexcept{                           //копирующий оператор 
        m_rows = copy.m_rows;
        m_columns = copy.m_columns;
        if(!this->empty()){
                delete[] m_ptr;
        }
        if(copy.empty()){
                m_ptr = nullptr;
        }
        else{
                m_ptr = new double[copy.m_rows * copy.m_columns];
                for (int i = 0; i < copy.m_columns * copy.m_rows; i++)
                {
                        m_ptr[i] = copy.m_ptr[i];
                }
        }
        return *this;
}

Matrix &Matrix::operator=(Matrix &&moved) noexcept{                               //перемещающий оператор
        if (&moved != this)
        {
                m_ptr = moved.m_ptr;
                m_columns = moved.m_columns;
                m_rows = moved.m_rows;
                moved.m_ptr = nullptr;
                moved.m_columns = 0;
                moved.m_rows = 0;
        }
        return *this;
}

double &Matrix::operator () (int i, int j){
        if(this->empty()){
                throw runtime_error("Error: matrix is empty.\n");
        }
        else if(i < 0 || j < 0 || i >= m_rows || j >= m_columns){
                throw runtime_error("Error: index out of range.\n");
        }
        return *(this->m_ptr + i * this->m_columns + j);
}

double Matrix::operator () (int i, int j) const{
        if(this->empty()){
                throw runtime_error("Error: matrix is empty.\n");
        }
        else if(i < 0 || j < 0 || i >= m_rows || j >= m_columns){
                throw runtime_error("Error: index out of range.\n");
        }
        return *(this->m_ptr + i * this->m_columns + j);
}

Matrix Matrix::operator + (const Matrix& right) const{
        Matrix sum(*this);
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.m_rows || m_columns != right.m_columns){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }
        for(int i = 0; i < m_columns*m_rows; i++){
                sum.m_ptr[i] += right.m_ptr[i];
        }
        return sum;
}

Matrix& Matrix::operator += (const Matrix& right){
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.m_rows || m_columns != right.m_columns){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }
        for(int i = 0; i < m_columns*m_rows; i++){
                m_ptr[i] += right.m_ptr[i];
        }
        return *this;
}

Matrix Matrix::operator - (const Matrix& right) const{
        Matrix dif(*this);
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.m_rows || m_columns != right.m_columns){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }
        for(int i = 0; i < m_columns*m_rows; i++){
                dif.m_ptr[i] -= right.m_ptr[i];
        }
        return dif;
}

Matrix& Matrix::operator -= (const Matrix& right){
        if(this->empty() || right.empty()){
                throw runtime_error("Error: empty matrix.\n");
        }
        else if(m_rows != right.m_rows || m_columns != right.m_columns){
                throw runtime_error("Error: matrix 1 and 2 have different shapes.\n");
        }
        for(int i = 0; i < m_columns*m_rows; i++){
                m_ptr[i] -= right.m_ptr[i];
        }
        return *this;
}

Matrix Matrix::operator * (const Matrix& other) const{ //Matrix * Matrix
        if ((this->empty()) || (other.empty())){
		throw std::runtime_error("\nError: empty matrix multiplication.\n");
        }
	else if (m_columns != other.m_rows){
		throw std::runtime_error("\nError: can't multiply matrixes.\n");
        }
        Matrix tmp(this->m_rows, other.m_columns);
        int counter = -1;
        for (int i = 0; i < tmp.m_columns * tmp.m_rows; i++)
        {
            tmp.m_ptr[i] = 0;
            if (!(i % tmp.m_columns))
                counter++;
            for (int j = 0; j < this->m_columns; j++)
            {
                tmp.m_ptr[i] += this->m_ptr[counter * this->m_columns + j] * other.m_ptr[i % other.m_columns + j * other.m_columns];
            }
        }
        return tmp;

}

Matrix Matrix::operator * (const double& num) const{ // Matrix * double
        if(this->empty()){
                throw runtime_error("\nError: empty matrix.\n");
        }
        Matrix mult(*this);
        for(int i = 0; i < mult.m_columns*mult.m_rows; i++){
                mult.m_ptr[i] = num * mult.m_ptr[i];
        }
        return mult;
}

Matrix& Matrix::operator *= (const Matrix& other){ // Matrix *= Matrix
        if(this->empty() || other.empty()){
                throw runtime_error("\nError: matrix is empty.\n");
        }
        else if(this->m_columns != other.m_rows){
                throw runtime_error("\nError: matrix 1 and 2 have different shapes.\n");
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


Matrix operator * (double num, const Matrix &rval){ // double * Matrix
        if(rval.empty()){
                throw runtime_error("\nError: matrix is empty.\n");
        }
        Matrix mult(rval);
        for(int i = 0; i < rval.m_columns*rval.m_rows; i++){
                mult.m_ptr[i] = num * rval.m_ptr[i];
        }
        return mult;
}

Matrix& Matrix::operator *= (double num){ // Matrix *= double
        if(this->empty()){
                throw runtime_error("\nError: empty matrix.\n");
        }
        for(int i = 0; i < this->m_columns*this->m_rows; i++){
                this->m_ptr[i] = this->m_ptr[i] * num;
        }
        return *this;
}

bool operator == (const Matrix& left, const Matrix& right) noexcept{ // Matrix == Matrix
        if(left.m_columns != right.m_columns || left.m_rows != right.m_rows){
                return false;
        }
        if(!left.empty() && !right.empty()){
                for(int i = 0; i < left.m_columns*left.m_rows; i++){
                        if(left.m_ptr[i] != right.m_ptr[i]) return false;
                }
        }
        return true;
}

bool operator != (const Matrix& left, const Matrix& right) noexcept{ // Matrix != Matrix
        if(left.m_columns != right.m_columns || left.m_rows != right.m_rows){
                return true;
        }
        if(!left.empty() && !right.empty()){
                for(int i = 0; i < left.m_columns*left.m_rows; i++){
                        if(left.m_ptr[i] != right.m_ptr[i]) return true;
                }
        }
        return false;
}

double Matrix::norm() const noexcept{ //Норма Фробениуса
        double norma = 0;
        for(int i = 0; i < this->m_columns*this->m_rows; i++){
                norma += (this->m_ptr[i]) * (this->m_ptr[i]);
        }
        norma = pow(norma, 0.5);
        return norma;
}

double Matrix::trace() const{ //След матрицы
        if(this->m_columns != this->m_rows){
                throw runtime_error("\nError: Matrix is not square.\n");
        }
        double traceM = 0;
        for(int i = 0; i < this->m_rows; i++){
                traceM += this->m_ptr[i*this->m_columns + i];
        }
        return traceM;
}


int width(double number) { //Вспомогательная функция для нахождения числа символов
        int width = 0;
        ostringstream stringNumber; //воспользуемся классом ostringstream
        stringNumber << setprecision(10) << number; //в strs теперь double не более 10 знаков
        string str = stringNumber.str(); //double теперь в виде string
        width = (int)str.size(); //получаем кол-во символов в строке
        return width;
}

ostream& operator << (ostream& os, const Matrix& M) noexcept { //перегрузка оператора вывода
        if (!M.empty()) {
                int* max_width;
                max_width = new int[M.m_columns]; //в массиве будут самые длинные числа в столбцах
                for (int j = 0; j < M.m_columns; j++) { 
                        max_width[j] = width(M.m_ptr[j]);
                        for (int i = 1; i < M.m_rows; i++) { //проходим по каждому ряду в столбце
                                if (max_width[j] < width(M.m_ptr[i * M.m_columns + j])) //если меньше макс, макс = width
                                        max_width[j] = width(M.m_ptr[i * M.m_columns + j]); //определяем самое "длинное" число в каждом столбце
                } }

                for (int i = 0; i < M.m_rows; i++) {
                        os << "|";
                        for (int j = 0; j < M.m_columns; j++) {
                                os << setw(max_width[j] + 1) << setprecision(10) << right << M.m_ptr[i * M.m_columns + j];
                        }
                        os << " |\n";
                }
                delete[] max_width;
        }
        else os << "( )\n";
        return os << endl;
};

Matrix transpose(const Matrix& m) noexcept{ //функция транспонирования матрицы
        if(m.empty()) {
                cout << "\nYou transpose empty matrix(no changes).\n";
                return Matrix();
        }
        Matrix T(m.m_columns, m.m_rows);
        for(int i = 0; i < m.m_columns; i++){
                for(int j = 0; j < m.m_rows; j++){
                        T.m_ptr[j*m.m_columns+i] = m.m_ptr[i*m.m_rows+j];
                }
        }
        return T;
}

Matrix concatenate(const Matrix& m1, const Matrix& m2){ //функция конкатенации
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

Matrix power(const Matrix &m, int degree){ //Возведение матрицы в степень
        if(m.m_columns != m.m_rows) throw runtime_error("\nError: matrix is not square shape.\n");
        if(m.m_columns == 0 || m.m_rows == 0) throw runtime_error("\nError: matrix is empty.\n");
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

Matrix &Matrix::gauss_forward() noexcept { //Прямой ход метода Гаусса
        int shift = 0;
        for (int i = 0; i < m_rows - 1; i++) { // m_rows - 1, т. к. не проверяем послелнюю строку
            if (shift + i >= m_columns) break; //shift - смещение для треугольного вида
            if (is_equal(m_ptr[i * m_columns + i + shift], 0)) { //маленькая ошибка из-за типа double 
                for (int p = i + 1; p < m_rows; p++){     //идем по всем следующим за i строкам
                    if (!is_equal(m_ptr[p * m_columns + i + shift], 0)) {
                        for (int g = i + shift; g < m_columns; g++) { //если первый элемент в следующей строке не 0
                            double element = m_ptr[i * m_columns + g]; //то меняем строки местами
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
                double coefficient = m_ptr[j * m_columns + i + shift] / m_ptr[i * m_columns + i + shift];
                if (coefficient == 0) continue; //если элемент в строке уже ноль
                for (int k = i; k < m_columns; k++) {
                    m_ptr[j * m_columns + k + shift] -= m_ptr[i * m_columns + k + shift] * coefficient; //хотя бы 1 элемент в строке теперь 0
                }
            }
        }
        return *this;
    }

Matrix &Matrix::gauss_backward() {
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
        tmp.m_ptr = new double[tmp.m_columns * tmp.m_rows];
        for (int i = 0; i < tmp.m_columns * tmp.m_rows; i++){
            if (i % tmp.m_columns == 0 && i != 0)
                real_shift += shift;
            tmp.m_ptr[i] = this->m_ptr[i + real_shift];
        }
        if (!tmp.det()){
            throw runtime_error("\nError! Determinant = 0.\n");
        }

        for (int i = m_rows - 1; i >= 0; i--){
            double norm_coef = this->m_ptr[i * this->m_columns + i];
            this->m_ptr[i * this->m_columns + i] /= norm_coef;
            for (int j = this->m_rows; j < this->m_columns; j++){
                this->m_ptr[i * this->m_columns + j] /= norm_coef;
            }
            if (i != 0){
                for (int j = i - 1; j >= 0; j--){
                    double coef = this->m_ptr[j * this->m_columns + i];
                    this->m_ptr[j * this->m_columns + i] -= this->m_ptr[i * this->m_columns + i] * coef;
                    for (int k = this->m_rows; k < this->m_columns; k++){
                        this->m_ptr[j * this->m_columns + k] -= this->m_ptr[i * this->m_columns + k] * coef;
                    }
                }
            }
        }
        return *this;
    }

double Matrix::det() const{
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
                            double val = tmp.m_ptr[i * this->m_columns + g];
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
                double coef = tmp.m_ptr[j * this->m_columns + i] / tmp.m_ptr[i * this->m_columns + i];
                for (int k = i; k < m_columns; k++){
                    tmp.m_ptr[j * this->m_columns + k] -= tmp.m_ptr[i * this->m_columns + k] * coef;
                }
            }
        }
        double res = tmp.m_ptr[0];
        for (int i = tmp.m_columns; i < this->m_rows * this->m_columns; i++) {
            if (i % tmp.m_columns == i / tmp.m_columns)
                res *= tmp.m_ptr[i];
        }
        return res * pow((-1), degree);
}

Matrix invert(const Matrix &mat){
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

Matrix solve(const Matrix &A, const Matrix &f) {
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


int Matrix::rank() const noexcept{ //ранг матрицы
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

Matrix Matrix::operator - () const noexcept{ //унарный минус  
    Matrix result(m_rows, m_columns);
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_columns; j++) {
            result(i, j) = -(*this)(i, j);
        }
    }
    return result;
}

double Matrix::minor(int num_row, int num_col){ //минор матрицы (индексация с 0)
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

Matrix vertical_concatenate(const Matrix &m1, const Matrix &m2){ //вертикальная конкатенация
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


}