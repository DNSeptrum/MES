#include "obliczenia.h"

using namespace std;

double gauss_legendre_1D(int n, double (*F)(double)) 
{
    double x[5][5] = 
    { 
        {0},
        { -sqrt(1.0 / 3.0),sqrt(1.0 / 3.0)},
        {-sqrt(3.0/5.0),0,sqrt(3.0 / 5.0)},
        {( - sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))) ,-sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) ,sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))},
        {-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459} 
    };
    double w[5][5] = 
    {
        {2},
        {1,1},
        {5.0/9.0,8.0/9.0,5.0 / 9.0},
        {(18+sqrt(30))/36,(18 - sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 + sqrt(30)) / 36},
        {0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851} 
    };

    double integral = 0.0;
    double x1 = 0.0;

    for (int i = 0; i < n; i++) 
    {
          x1 = x[n - 1][i];
        cout << (*F)(x1) << endl;
          integral += w[n-1][i] * (*F)(x1);
    }

    return integral;
}

double gauss_legendre_2D(int n, double (*F)(double,double)) {
    double x[5][5] = 
    {
        {0},
        { -sqrt(1.0 / 3.0),sqrt(1.0 / 3.0)},
        {-sqrt(3.0 / 5.0),0,sqrt(3.0 / 5.0)},
        {(-sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))) ,-sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) ,sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))},
        {-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459}
    };
    double w[5][5] = 
    {
        {2},
        {1,1},
        {5.0 / 9.0,8.0 / 9.0,5.0 / 9.0},
        {(18 + sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 + sqrt(30)) / 36},
        {0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851}
    };

    double integral = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;

    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            x1 = x[n - 1][i];
            y1 = x[n - 1][j];
            integral += w[n - 1][i] * w[n - 1][j] * (*F)(x1, y1);
        }
    }
    return integral;
}

double f(double x) {
    return 5 * (pow(x, 2)) + 3 * x + 6;
}
double f(double x, double y) {
    return 5 * (pow(x, 2)) * (pow(y, 2)) + 3 * x * y + 6;
};