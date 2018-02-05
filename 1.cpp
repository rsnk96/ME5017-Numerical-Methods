/* 
Written by R S Nikhil Krishna, ME14B088, IIT Madras
Solution to Programming Assignment 1

This program performs Gaussian Elimination to obtain roots of Ax=b
b is handled as a column matrix, and hence a nx1 2D matrix

Both questions' solutions are included
To obtain solution of Q2, just use the setValues function
 */

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;
typedef vector<float> vf;
typedef vector<vf> vvf;

void printArray(vvf A)
{
    cout << "Your Array is:\n";
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A[i].size(); j++)
            cout << A[i][j] << '\t';
        cout << endl;
    }
}

vvf forward_elimination(vvf &A, vvf &b)
{
    vvf L(A.size(), vf(A.size(), 0));

    for (size_t j = 0; j < A[0].size(); j++)
    {
        L[j][j] = 1;
        for (size_t i = j + 1; i < A.size(); i++)
        {
            L[i][j] = A[i][j] / A[j][j];
            for (size_t k = j; k < A[i].size(); k++)
                A[i][k] -= L[i][j] * A[j][k];
            b[i][0] -= L[i][j] * b[j][0];
        }
    }
    return L;
}

void backward_substitution(vvf &A, vvf &b)
{
    for (int i = A.size() - 1; i >= 0; i--)
    {
        for (size_t j = i+1; j < A.size(); j++)
            b[i][0] -= A[i][j] * b[j][0];
        b[i][0] /= A[i][i];
    }
}

double calc_norm(vvf b)
{
    double sum = 0;
    for (size_t i = 0; i < b.size(); i++)
        sum += b[i][0] * b[i][0];
    cout << "\nSum of squares of x values is: " << sum <<'\n';
    return sum;
}

void setValues(vvf &A, vvf &b)
{
    A={
        {-4,1,1,0,0,0,0,0},
        {1,0,-4,1,1,0,0,0},
        {0,0,1,0,-4,1,1,0},
        {2,-4,0,1,0,0,0,0},
        {0,1,2,-4,0,1,0,0},
        {0,0,0,1,2,-4,0,1},
        {0,0,0,0,2,0,-9,1},
        {0,0,0,0,0,2,2,-9}
    };
    b={{-1000},{-500},{-500},{-500},{0},{0},{-2000},{-1500}};
}

int main(int argc, char **argv)
{
    cout << fixed;
    cout << setprecision(3);
    int n = 32;
    if (argc == 1)
    {
        cout << "Please enter the dimensions of the square matrix: ";
        cin >> n;
    }
    else
        n = atoi(argv[1]);

    vvf A(n, vf(n, 0));
    vvf b(n, vf(1, 1));
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            A[i][j] = max(i + 1, j + 1);
    setValues(A,b);

    vvf L = forward_elimination(A, b);
    cout << "\nForward Elimination done. ";
    printArray(A);

    cout << "\nThe L Matrix is computed. ";
    printArray(L);

    backward_substitution(A, b);
    cout << "\nBackward Substition done. ";
    printArray(b);

    calc_norm(b);
    return 0;
}