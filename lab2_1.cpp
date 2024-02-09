#include <iostream>
#include "Matrix.h"
#include "test.h"
int main()
{
    // test_for_constructors();
    // test_for_methods();
    // test_for_functions();
    // test_for_operators();
    // dops();
    try
    {
        linalg::Matrix<string> m = {{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}} ;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }

    try
    {linalg::Matrix<string> m = {{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}} ;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }

    try
    {
    linalg::Matrix<int> m2 = {{1,2,3},{4,5,6}} ;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
    
    //m+=m2;
    //cout << m;



    return 0;
}