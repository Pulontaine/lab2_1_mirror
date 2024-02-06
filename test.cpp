#include "test.h"
#include "Matrix.h"
#include <iostream>
#include <stdexcept>
using namespace std;
using namespace linalg;
#define T int


//======================================  test_for_constructors
void constructor_by_default(){
	cout << "\n\n________________________________________________________\n";
	cout << "linalg::Matrix m0;\n";
	cout << "Take 0 arguments and make empty matrix.\n";
	linalg::Matrix<T> m0;
	cout << "\nm0.empty(): " << m0.empty() << endl;
}

// void constructor_one_arg(){ 
//     cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0(m_rows);\n";
// 	cout << "Take (int m_rows) and make a m_rows-high column matrix.\n";
// 	linalg::Matrix m0(10);
// 	cout << "Exceptions:\nm_rows <= 0.\n";
// 	try{
// 		linalg::Matrix m1(-1);
// 	}
// 	catch (runtime_error exc) {
// 		cerr << exc.what();
// 	}
// }

// void constructor_two_arg(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0(m_rows, m_columns);\n";
// 	cout << "Take (int m_rows, int m_columns) and make a m_rows-high and m_columns-wide matrix.\n";
// 	linalg::Matrix m0(10, 10);
// 	cout << "Exceptions:\n1. m_rows <= 0.\n2. m_columns <= 0.\n";
// 	try{
// 		linalg::Matrix m1(-1, -3);
// 	}
// 	catch (runtime_error exc) {
// 		cerr << exc.what();
// 	}
// }

// void constructor_init_list(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0(std::initializer_list<double>);\n";
// 	cout << "Take (std::initializer_list<double>) and make a matrix with elements in the column.\n";
// 	linalg::Matrix m0{1.5, 10.5, 14.5, 5.9};
// }

// void constructor_init_init_list(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0(std::initializer_list<std::initializer_list<double>>);\n";
// 	cout << "Take (std::initializer_list<std::initializer_list<double>>) and make a matrix with these elements.\n";
// 	linalg::Matrix m0{{1.5, 10.5}, {14.5, 5.9}};
// 	try{
// 		linalg::Matrix m1{{1.5, 10.5}, {14.5, 5.9, 4.5}};
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what();
// 	}
// }

// void constructor_copy(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0(const Matrix &other);\n";
// 	cout << "Take (const Matrix &other) and make copy of a matrix with these elements.\n";
// 	linalg::Matrix m0;
// 	linalg::Matrix m1(m0);
// }

// void constructor_move(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0(Matrix &&move);\n";
// 	cout << "Take (Matrix &&move) and make copy of a matrix with these elements.\n";
// 	linalg::Matrix m0;
// 	linalg::Matrix m1(move(m0));
// }


// //===============================  test_for_methods
// void rows_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.rows();\n";
// 	cout << "Return quantity of m_rows in matrix.\n";
// 	linalg::Matrix m0(3,4);
// 	cout << "linalg::Matrix m0(3,4);\n" << "Rows: " << m0.rows();
// }

// void columns_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.columns();\n";
// 	cout << "Return quantity of m_columns in matrix.\n";
// 	linalg::Matrix m0(3,4);
// 	cout << "linalg::Matrix m0(3,4);\n" << "\nColumns: " << m0.columns();
// }

// void reshape_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.reshape(m_rows, m_columns);\n";
// 	cout << "Return (Matrix& M). Change quantity of rows and columns.\n";
// 	cout << "Exception:\n1. Reshaped matrix will consist more/less elements.\n";
// 	linalg::Matrix m0 = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
// 	cout << m0 << endl;
// 	cout << "m0.reshape(16, 1):\n" << m0.reshape(16, 1) << endl;
// 	try{
// 		m0.reshape(16, 2);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}
// }

// void norm_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.norm();\n";
// 	cout << "Return double - Frobenius norm of the matrix\n";
// 	linalg::Matrix m0 = {{1, 2}, {3, 4}};
// 	cout << m0 << "\nNorm: " << m0.norm() << endl;
// }

// void trace_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.trace();\n";
// 	cout << "Return double - trace of the matrix.\n";
// 	cout << "Exception:\n1. Matrix is not square.\n";
// 	linalg::Matrix m1 = {{1, 2}, {3, 4}};
// 	cout << m1 << "\nTrace: " << m1.trace() << endl;
// 	linalg::Matrix m0 = {{1, 2, 3}, {3, 4, 5}};
// 	try{
// 		m0.trace();
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}

// }

// void det_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.det();\n";
// 	cout << "Return double - determinant of the matrix.\n\n";
// 	cout << "Exception:\n1. Matrix is not square.\n2. Matrix is empty.\n\n";
// 	linalg::Matrix m0 = {{10, 20}, {3, 4}};
// 	linalg::Matrix m1 = {{1, 20}, {3, 4, 5}};
// 	cout << "\nm0:\n" << m0 << "\nDet m0: " << m0.det();
// 	cout << "\nm1:\n" << m1 << "\nDet m1: \n";
// 	try{
// 		m1.det();
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}

// 	linalg::Matrix m2;
// 	cout << "\nm2:\n" << m2 << "\nDet m2: \n";
// 	try{
// 		m2.det();
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}

// }

// void gauss_f_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.gauss_forward();\n";
// 	cout << "Return upper triangular matrix via Forward Gauss method.\n\n";
// 	linalg::Matrix m0 = {{1, 2, 3}, {2, 4, 6}, {1, 3, 2}};
// 	cout << m0 << "\nm0.gauss_forward():\n" << m0.gauss_forward() << endl;
// }

// void gauss_b_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.gauss_backward();\n";
// 	cout << "Return (E|F) matrix via Backward Gauss method.\n";
// 	linalg::Matrix m0 = {{1, 2, 3}, {4, 5, 6}, {1, 3, 2}};
// 	cout << "\nm0:\n" << m0;
// 	cout << "\nm0.gauss_backward():\n";
// 	try{
// 		m0.gauss_backward();
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}
	
// }

// void rank_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "linalg::Matrix m0.rank();\n";
// 	cout << "Return int - rank of the matrix.\n";
// 	linalg::Matrix m1 = {{1, 20}, {3, 4}};
// 	cout << "\nm1:\n" << m1 << "Rank: " << m1.rank();	
// 	linalg::Matrix m2 = {{1, 2}, {2, 4}};
// 	cout << "\n\nm2:\n" << m2 << "Rank: " << m2.rank();	
// }


// //===========================  test_for_functions
// void test_concatenate(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "concatenate(Matrix& M, Matrix& M1);\n";
// 	cout << "Return matrix of concatenates M and M1.\n";
// 	cout << "Exceptions:\n1. Different rows.\n2. Empty matrix.\n\n";
// 	linalg::Matrix m1 = {{1, 0}, {1, 0}};
// 	linalg::Matrix m2 = {{2, 3}, {4, 5}};
// 	cout << "m1:\n" << m1 <<"\nm2:\n" << m2 <<"\nConcatenate matrix:\n" << concatenate(m1, m2) << endl;
// 	linalg::Matrix m3;
// 	try{
// 		concatenate(m3, m1);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << '\n';
// 	}
// 	linalg::Matrix m4(10, 2);
// 	try{
// 		concatenate(m4, m1);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << '\n';
// 	}
// }

// void test_transpose(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "transpose(const Matrix& M);\n";
// 	cout << "Return transposed matrix of M.\n\n";
// 	linalg::Matrix m1{{1,2,3}, {4,5,6}};
// 	cout << "m1:\n" << m1 <<"Transposed matrix m1:\n" << transpose(m1) << endl;
// }

// void test_invert(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "invert(const Matrix& M);\n";
// 	cout << "Return inver matrix of M.\n";
// 	cout << "Exception:\n1. Det == 0;\n\n";
// 	linalg::Matrix m1{{1,2,3}, {4,5,10}, {3,4,1}};
// 	cout << "m1:\n" << m1 << "Invert matrix:\n" << invert(m1) << endl;
// 	linalg::Matrix m2{{1,2}, {1,2}};
// 	cout << "\nm2:\n" << m2 << "\ninvert(m2):\n";
// 	try{
// 		invert(m2);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}
// }

// void test_power(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "power(const Matrix &m, int degree);\n";
// 	cout << "Return matrix in pow(degree).\n";
// 	cout << "Exception:\n1. matrix is not square shape.\n\n";

// 	linalg::Matrix m3;
// 	cout << "m3:\n" << m3 <<"matrix in power:\n";
// 	try{
// 		power(m3, 2);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}
// 	linalg::Matrix m2{{9,6,6}, {12,2,4}, {8,6,5}};
// 	cout << "m2:\n" << m2 <<"matrix in power 0:\n" << power(m2, 0) << endl;

// 	linalg::Matrix m1{{1,2,3}, {4,5,1}, {1,3,6}};
// 	cout << "m1:\n" << m1 <<"matrix in power 3:\n" << power(m1, 3);

// 	linalg::Matrix m4{{1,2}, {1,4}, {3,5}};
// 	cout << "\nm4:\n" << m4;
// 	cout << "\npower(m4, 2):\n";
// 	try{
// 		power(m4, 2);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}
// }

// void test_solve(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "solve(const Matrix &A, const Matrix &f);\n";
// 	cout << "Return solution of the equation Ax = f.\n\n";
// 	//cout << "Exceptions:\n1. "

// 	linalg::Matrix m0 = {{3, 2, -5}, {2, -1, 3}, {1, 2, -1}};
//     linalg::Matrix f0 = {-1, 13, 9};
// 	cout << "m0:\n" << m0 << "\nf0:\n" << f0 << "\nsolve(m0,f0):\n" << solve(m0, f0) << endl;

// 	linalg::Matrix m1 = {{1, 2}, {1, 2}};
//     linalg::Matrix f1 = {6, 8};
// 	cout << "m1:\n" << m1 << "\nf1:\n" << f1 << "\nsolve(m1,f1):\n";
// 	try{
// 		solve(m1,f1);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << '\n';
// 	}

// 	linalg::Matrix m2 = {{1, 3}, {5, 7}, {9, 9}};
//     linalg::Matrix f2 = {{5, 6, 7, 8}};
// 	cout << "m2:\n" << m2 << "\nf2:\n" << f2 << "\nsolve(m2,f2):\n";
// 	try{
//         solve(m2, f2);
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }

// 	linalg::Matrix m3 = {{5, 9}, {0, -3}};
//     linalg::Matrix f3 = {5, 8, 1, 2};
// 	cout << "m3:\n" << m3 << "\nf3:\n" << f3 << "\nsolve(m3,f3):\n";
// 	try{
//         solve(m3, f3);
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }

// 	linalg::Matrix m4, f4;
// 	cout << "m4:\n" << m4 << "\nf4:\n" << f4 << "\nsolve(m4,f4):\n";
//     try{
//         solve(m4, f4);
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }


// //===========================  test_for_operators
// void test_copy(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix& operator = (const Matrix& copy);\n";
// 	cout << "Copy matrix, assign rows and columns of lvalue.\n";
// 	linalg::Matrix m0 = {{1, 2}, {3,4}};
// 	linalg::Matrix m1 = {{3, 5}};
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1;
// 	m1 = m0;
// 	cout << "\nm1 = m0;\nm1: " << m1 << endl;
// }

// void test_move(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix& operator = (const Matrix&& moved);\n";
// 	cout << "Move matrix, assign rows and columns of lvalue.\n";
// 	linalg::Matrix m0 = {{1, 2}, {3,4}};
// 	linalg::Matrix m1 = {{2,6}};
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1 << endl;
// 	m1 = move(m0);
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1 << endl;
// }

// void test_stream(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "ostream& operator << (ostream& os, const Matrix& out);\n";
// 	cout << "Print elements of matrix.\n";
// 	linalg::Matrix m1;
//     cout << "m1:\n" << m1 << endl;
// 	linalg::Matrix m2 = {{1,2}, {8,9}};
//     cout << "m2:\n" << m2 << endl;
// }

// void test_call(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "double& Matrix::operator (i,j);\n";
// 	cout << "Return [i,j] element of the matrix as reference.\n";
// 	cout << "Exceptions:\n1. Matrix is empty.\n2. Index is out of range.\n";
// 	linalg::Matrix m2 = {{1,2}, {8,9}};
//     cout << "m2(0,1):\n" << m2(0,1) << endl;
// 	cout << "m2(4,5):\n";
// 	try{
// 		m2(4,5);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}
// 	linalg::Matrix m3;
// 	cout << "m3(1,1):\n";
// 	try{
// 		m3(1,1);
// 	}
// 	catch(runtime_error exc){
// 		cerr << exc.what() << endl;
// 	}
// }

// void test_plus(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix operator + (const Matrix& right);\n";
// 	cout << "Return rvalue. Sum each elements of two matrixes.\n";
// 	cout << "Exception:\n1. Different shapes.\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}};
//     linalg::Matrix m1 = {{-1, 2}, {5, 8}};
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1 << "\nm1+m0:\n" << m1+m0;
// 	linalg::Matrix m2 = {1, 2, 3};
//     cout << "\nm2:\n"<< m2;
// 	cout << "\nm0+m2:\n";
//     try{
//         m0 + m2;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_eqplus(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix& operator += (const Matrix& right);\n";
// 	cout << "Return lvalue. Sum each elements of two matrixes.\n";
// 	cout << "Exception:\n1. Different shapes.\n\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}};
//     linalg::Matrix m1 = {{-1, 2}, {5, 8}};
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1;
// 	m1+=m0;
// 	cout << "\nm1+=m0;\nm1:\n" << m1;
// 	linalg::Matrix m2 = {1, 2, 3};
//     cout << "\nm2:\n"<< m2;
// 	cout << "\nm0+=m2\n";
//     try{
//         m0+=m2;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_minus(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix operator - (const Matrix& right);\n";
// 	cout << "Return rvalue. Substract each elements of two matrixes.\n";
// 	cout << "Exception:\n1. Different shapes.\n\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}};
//     linalg::Matrix m1 = {{-1, 2}, {5, 8}};
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1 << "\nm1-m0:\n" << m1-m0;
// 	linalg::Matrix m2 = {1, 2, 3};
//     cout << "\nm2:\n"<< m2;
// 	cout << "\nm0-m2:\n";
//     try{
//         m0 - m2;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_eqminus(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix operator -= (const Matrix& right);\n";
// 	cout << "Return lvalue. Substract each elements of two matrixes.\n";
// 	cout << "Exception:\n1. Different shapes.\n\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}};
//     linalg::Matrix m1 = {{-1, 2}, {5, 8}};
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1;
// 	m1-=m0;
// 	cout << "m1-=m0" <<"\nm1 now:\n\n" << m1;
// 	linalg::Matrix m2 = {1, 2, 3};
//     cout << "\nm2:\n"<< m2;
// 	cout << "\nm0-=m2\n";
//     try{
//         m0-=m2;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_mul(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix operator * (const Matrix& right);\n";
// 	cout << "Multiplication of two matrixes.\n";
// 	cout << "Exception:\n1. Different shapes.\n2. Empty matrix.\n\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}, {6, 7}};
//     linalg::Matrix m1 = {{-1, 2, 5}, {1, 5, 8}};
// 	cout << "m0:\n" << m0 << "\nm1:\n" << m1 << "\nm1 * m0:\n" << m1*m0;
// 	linalg::Matrix m2 = {1, 2, 3};
//     cout << "\nm2:\n"<< m2;
// 	cout << "\nm0 * m2:\n";
//     try{
//         m0*m2;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// 	linalg::Matrix m3;
// 	cout << "\nm3:\n" << m3 << "\nm0 * m3:\n";
// 	try{
//         m0*m3;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_Num_Mat(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix operator * (double num, const Matrix& rval);\n";
// 	cout << "Multiplication of each element of matrix with double.\n";
// 	cout << "Exception:\n1. Empty matrix.\n\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}, {6, 7}};
// 	cout << "m0:\n" << m0 << "\n5 * m0:\n" << 5 * m0;
// 	linalg::Matrix m3;
// 	cout << "\nm3:\n" << m3 << "\n5 * m3:\n";
// 	try{
//         5*m3;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_Mat_Num(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix operator * (const double& num);\n";
// 	cout << "Multiplication of each element of matrix with double.\n";
// 	cout << "Exception:\n1. Empty matrix.\n\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}, {6, 7}};
// 	cout << "m0:\n" << m0 << "\nm0 * 5:\n" << m0 * 5;
// 	linalg::Matrix m3;
// 	cout << "\nm3:\n" << m3 << "\nm3 * 5:\n";
// 	try{
//         m3 * 5;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_eq_mul_Mat(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix& operator *= (const Matrix& right)\n";
// 	cout << "Multiplication of two matrixes.\n";
// 	cout << "Exception:\n1. Different shapes.\n2. Empty matrix.\n\n";
	
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}, {6, 7}};
// 	linalg::Matrix m1 = {{4, 5}, {2, 3}, {6, 7}};
// 	linalg::Matrix m3 = {{3, 5}, {3, 5}};

// 	cout << "m0:\n" << m0 << "\nm3:\n" << m3 << "\nm0 *= m3;\nm0:\n";
// 	m0 *= m3;
// 	cout << m0;

// 	cout << "\n\nm0:\n" << m0 << "\nm1:\n" << m1 << "\nm0 *= m1:\n";
// 	try{
//         m0*=m1;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_eq_mul_double(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix& operator *= (const double& num);\n";
// 	cout << "Multiplication of double and matrix.\n";
// 	cout << "Exception:\n1. Empty matrix.\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}, {6, 7}};
// 	double var = 5.6;
// 	cout << "\nm0:\n" << m0 << "\nm0*=5:\n"; 
// 	m0*=var;
// 	cout << m0;
// 	linalg::Matrix m1;
// 	try{
//         m1*=5;
//     }
//     catch (const runtime_error &e){
//         std::cerr << e.what() << '\n';
//     }
// }

// void test_eq(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "bool operator == (const Matrix& left, const Matrix& right);\n";
// 	cout << "Return bool: True, if matrixes are equal.\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}, {6, 7}};
// 	linalg::Matrix m1 = {{3, 7}, {2, 3}, {6, 7}};
// 	cout << "\nm0:\n" << m0 << "\nm1:\n" << m1;
// 	cout << "\nm0==m1:\n";
// 	bool q = m0==m1;
// 	cout << q;
// }

// void test_not_eq(){cout << "\n\n________________________________________________________\n";
// 	cout << "bool operator == (const Matrix& left, const Matrix& right);\n";
// 	cout << "Return bool: True, if matrixes are equal.\n";
// 	linalg::Matrix m0 = {{3, 7}, {2, 3}, {6, 7}};
// 	linalg::Matrix m1 = {{3, 4}, {2, 3}, {6, 7}};
// 	cout << "\nm0:\n" << m0 << "\nm1:\n" << m1;
// 	cout << "\nm0!=m1:\n";
// 	bool q = m0!=m1;
// 	cout << q;
// }


// //=================================== тесты для дополнительного задания

// void minor_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "double Matrix::minor(int num_row, int num_col);\n";
// 	cout << "Return det of a matrix withou num_row and num_col.\n";
// 	linalg::Matrix m0 = {{3, 7, 4}, {1, 2, 3}, {2, 6, 7}};
// 	linalg::Matrix m1 = {{3, 4}, {2, 3}, {6, 7}};

// 	cout << "\nm1:\n" << m1 << "\nm0.minor(-1,0):\n";
// 	try{
//         m1.minor(-1, 0);
//     }
//     catch (runtime_error exc){
//         cerr << exc.what() << '\n';
//     }

// 	cout << "\nm0:\n" << m0 << "\nm0.minor(1,5):\n";
//     try{
//         m0.minor(1, 5);
//     }
//     catch (runtime_error exc){
//         cerr << exc.what() << '\n';
//     }

//     cout << "\nm1:\n" << m1 << "\nm1.minor(1,0):\n";
//     try{
//         m1.minor(1, 0);
//     }
//     catch (runtime_error exc){
//         cerr << exc.what() << '\n';
//     }

// }

// void concat_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix vertical_concatenate(const Matrix &m1, const Matrix &m2);\n";
// 	cout << "Return concatenated matrix (by columns).\n";
// 	linalg::Matrix m1 = {{3, 4}, {2, 3}, {6, 7}};
// 	linalg::Matrix m2 = {{3, 2}, {1, 5}};

// 	cout << "\nm1:\n" << m1 << "\nm2:\n" << m2 << "\nvertical_concatenate(m1, m2):\n";
// 	cout << vertical_concatenate(m1, m1);


// 	linalg::Matrix m3 = {{1, 4, 7}, {1, 3, 8}};
// 	cout << "\nm1:\n" << m1 << "\nm3:\n" << m3 << "\nvertical_concatenate(m1, m3):\n";
// 	try{
//         vertical_concatenate(m1, m3);
//     }
//     catch (runtime_error exc){
//         cerr << exc.what() << '\n';
//     }
    
//     linalg::Matrix m4, m5;
// 	cout << "\nm4:\n" << m4 << "\nm5:\n" << m5 << "\nvertical_concatenate(m4, m5):\n";
//     try{
//         vertical_concatenate(m4, m5);
//     }
//     catch (runtime_error exc){
//         cerr << exc.what() << '\n';
//     }

// }

// void un_minus_test(){
// 	cout << "\n\n________________________________________________________\n";
// 	cout << "Matrix Matrix::operator - () const;\n";
// 	cout << "Return matrix after unary minus.\n";
// 	linalg::Matrix m1 = {{3, 4}, {2, 3}, {6, 7}};
// 	cout << "\nm1:\n" << m1;
// 	cout << "\nm0 = -m1;\nm0:\n";
// 	linalg::Matrix m0 = -m1;
// 	cout << m0;
// }




// void dops(){
// 	minor_test();
// 	concat_test();
// 	un_minus_test();
// }

void test_for_constructors(){
    constructor_by_default();
//     constructor_one_arg();
// 	constructor_two_arg();
// 	constructor_init_list();
// 	constructor_init_init_list();
// 	constructor_copy();
// 	constructor_move();
}

// void test_for_methods(){
// 	rows_test();
// 	columns_test();
// 	reshape_test();
// 	norm_test();
// 	trace_test();
// 	rank_test();
// 	gauss_b_test();
// 	gauss_f_test();
// }

// void test_for_functions(){
// 	test_concatenate();
// 	test_transpose();
// 	test_invert();
// 	test_power();
// 	test_solve();
// }

// void test_for_operators(){
// 	test_copy();
// 	test_move();
// 	test_stream();
// 	test_call();
// 	test_plus();
// 	test_eqplus();
// 	test_minus();
// 	test_eqminus();
// 	test_mul();
// 	test_Num_Mat();
// 	test_Mat_Num();
// 	test_eq_mul_Mat();
// 	test_eq_mul_double();
// 	test_eq();
// 	test_not_eq();
// }