#include "obliczenia.h"

//ksi
double n1ksi(double eta) {
    return -0.25 * (1 - eta);
}
double n2ksi(double eta) {
    return 0.25 * (1 - eta);
}
double n3ksi(double eta) {
    return 0.25 * (1 + eta);
}
double n4ksi(double eta) {
    return -0.25 * (1 + eta);
}

//eta 
double n1eta(double ksi) {
    return -0.25 * (1 - ksi);
}
double n2eta(double ksi) {
    return -0.25 * (1 + ksi);
}
double n3eta(double ksi) {
    return 0.25 * (1 + ksi);
}
double n4eta(double ksi) {
    return 0.25 * (1 - ksi);
}

void gauss_legendre_2D_eta_ksi(int n, ElementUniversalny& wsk) {
    double x[4][4] =
    {
        {0},
        { -sqrt(1.0 / 3.0),sqrt(1.0 / 3.0)},
        {-sqrt(3.0 / 5.0),0,sqrt(3.0 / 5.0)},
        {(-sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))) ,-sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) ,sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))}
    };
    double w[4][4] =
    {
        {2},
        {1,1},
        {5.0 / 9.0,8.0 / 9.0,5.0 / 9.0},
        {(18 + sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 + sqrt(30)) / 36}
    };
    double y1 = 0.0;
    double x1 = 0.0;

    for (int i = 0; i < n; i++)
    {
        y1 = x[n - 1][i];

        for (int j = 0; j < n; j++)
        {
         x1 = x[n - 1][j];

          wsk.eta[0][j + i * n] = n1ksi(y1);
          wsk.eta[1][j + i * n] = n2ksi(y1);
          wsk.eta[2][j + i * n] = n3ksi(y1);
          wsk.eta[3][j + i * n] = n4ksi(y1);
          wsk.ksi[0][j + i * n] = n1eta(x1);
          wsk.ksi[1][j + i * n] = n2eta(x1);
          wsk.ksi[2][j + i * n] = n3eta(x1);
          wsk.ksi[3][j + i * n] = n4eta(x1);

        }
    }

}

void wys_eta( ElementUniversalny & ksieta,int n ){
    cout << endl<< "Wartosci Eta: " << endl;
    cout << "N1: ";
    for (int i = 0; i < n*n; i++)
    {
        cout << ksieta.eta[0][i] << " , ";
    }
    cout << endl;
    cout << "N2: ";
    for (int i = 0; i < n * n; i++)
    {

        cout << ksieta.eta[1][i] << " , ";
    }
    cout << endl;
    cout << "N3: ";
    for (int i = 0; i < n * n; i++)
    {

        cout << ksieta.eta[2][i] << " , ";
    }
    cout << endl;
    cout << "N4: ";
    for (int i = 0; i < n * n; i++)
    {

        cout << ksieta.eta[3][i] << " , ";
    }
}
void wys_ksi(ElementUniversalny& ksieta, int n) {
    cout << endl << "Wartosci Ksi: " << endl;
    cout << "N1: ";
    for (int i = 0; i < n * n; i++)
    {
        cout << ksieta.ksi[0][i] << " , ";
    }
    cout << endl;
    cout << "N2: ";
    for (int i = 0; i < n * n; i++)
    {

        cout << ksieta.ksi[1][i] << " , ";
    }
    cout << endl;
    cout << "N3: ";
    for (int i = 0; i < n * n; i++)
    {

        cout << ksieta.ksi[2][i] << " , ";
    }
    cout << endl;
    cout << "N4: ";
    for (int i = 0; i < n * n; i++)
    {

        cout << ksieta.ksi[3][i] << " , ";
    }
}