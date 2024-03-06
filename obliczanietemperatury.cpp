#include "obliczenia.h"

void obliczanie(Grind& grind, Global_Data& data) {

    double** macierz = new double*[data.Nodes_number];
    for (int i=0; i < data.Nodes_number; i++) {
        macierz[i] = new double[data.Nodes_number];
    }
	for (int k = 0; k < data.Nodes_number; k++) {
		for (int i = 0; i < data.Nodes_number; i++) {
			macierz[i][k] = 0;
		}
	}
	double* macierzHC = new double[data.Nodes_number];
	 grind.Temp = new double[data.Nodes_number];
	for (int i = 0; i < data.Nodes_number; i++) {
		macierzHC[i] = 0;
        grind.Temp[i] = 0;
	}

	for (int k = 0; k < data.Nodes_number; k++) {
		for (int i = 0; i < data.Nodes_number; i++) {
            grind.Temp[k] += grind.FinalHC[i][k]*grind.FinalP[i];

		}
	}
    for (int i = 0; i < data.Nodes_number; i++) {
        grind.Temp0[i] = grind.Temp[i];

    }

}

const double eps = 1e-12;

// Funkcja dokonuje rozkładu LU macierzy A
//----------------------------------------
bool ludist(int n, double** A)
{
    int i, j, k;

    for (k = 0; k < n - 1; k++)
    {
        if (fabs(A[k][k]) < eps) return false;

        for (i = k + 1; i < n; i++)
            A[i][k] /= A[k][k];

        for (i = k + 1; i < n; i++)
            for (j = k + 1; j < n; j++)
                A[i][j] -= A[i][k] * A[k][j];
    }

    return true;
}

// Funkcja wyznacza wektor X na podstawie A i X [ i ] 
//---------------------------------------------------
bool lusolve(int k, int n, double** A, double** X)
{
    int    i, j;
    double s;

    for (i = 1; i < n; i++)
    {
        s = 0;

        for (j = 0; j < i; j++) s += A[i][j] * X[j][k];

        X[i][k] -= s;
    }

    if (fabs(A[n - 1][n - 1]) < eps) return false;

    X[n - 1][k] /= A[n - 1][n - 1];

    for (i = n - 2; i >= 0; i--)
    {
        s = 0;

        for (j = i + 1; j < n; j++) s += A[i][j] * X[j][k];

        if (fabs(A[i][i]) < eps) return false;

        X[i][k] = (X[i][k] - s) / A[i][i];
    }

    return true;
}

int  odwracanie(int n, Grind& grind, Global_Data& data)
{
    double** A, ** X;
   // int n,
       int i, j;
    bool ok;
    n = data.Nodes_number;


    A = new double* [n];
    X = new double* [n];

    for (i = 0; i < n; i++)
    {
        A[i] = new double[n];
        X[i] = new double[n];
    }


    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) A[i][j] = grind.FinalHC[i][j];


    if (ludist(n, A))
    {


        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++) X[i][j] = 0;
            X[i][i] = 1;
        }


        ok = true;
        for (i = 0; i < n; i++)
            if (!lusolve(i, n, A, X))
            {
                ok = false;
                break;
            }
    }
    else ok = false;


    if (ok)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                //   cout << setw(10) << X[i][j] << " ";
             //  cout << endl;
                grind.FinalHC[i][j] = X[i][j];
        }
    }
    else cout << "DZIELNIK ZERO\n";


    for (i = 0; i < n; i++)
    {
        delete[] A[i];
        delete[] X[i];
    }
    delete[] A;
    delete[] X;

    return 0;
}

void min(Grind grind,Global_Data data) {
    double min = grind.Temp[0];
    double pom;
    for (int i = 1; i < data.Nodes_number; i++) 
    {
         pom = grind.Temp[i];
        if (min > pom)
            min = pom;
    }
    cout << fixed;
    cout << "min temperatura: " << setprecision(10) << min << endl;
}

void max(Grind grind, Global_Data data) {
    double max = grind.Temp[0];
    double pom;
    for (int i = 1; i < data.Nodes_number; i++) 
    {
        pom = grind.Temp[i];
        if (max < pom)
            max = pom;
    }
    cout << fixed;
    cout << "max temperatura: " << setprecision(10) << max << endl;
}