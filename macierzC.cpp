#include "obliczenia.h"

void macierzC(ElementUniversalny& ksieta, Elements& id, Grind& nodes, Global_Data& data, int n) {
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
		{1,1,1,1},
		{5.0 / 9.0,8.0 / 9.0,5.0 / 9.0},
		{(18 + sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 + sqrt(30)) / 36},
		{0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851}
	};
	double result[4] = { 0,0,0,0 };
	double detj = 0;
	double N[4];
	double ksi = 0.0, eta =0.0;
	double C[4][4];
	int wagi1 = 0, wagi2 = 0;

	int pom1 = 0;
	int pom2 = 0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			id.C[i][j] = 0;
		}
	}

	// jedno przejście 1 punkt
	for (int k = 0; k < n*n; k++) {

		// zamieniany 3 z 4 punktem
		
		if (k < n) { pom1 = 0; };
		if (k < 2 * n && k >= n) { pom1 = 1; };
		if (k < 3 * n && k >= 2 * n) { pom1 = 2; };
		if (k < 4 * n && k >= 3 * n) { pom1 = 3; };

		if (k < n) { pom2 = k; }
		if (k < 2 * n && k >= n) { pom2 = k - n; }
		if (k < 3 * n && k >= 2 * n) { pom2 = k - n * 2; }
		if (k < 4 * n && k >= 3 * n) { pom2 = k - n * 3; }
		
		ksi = x[n - 1][pom2];
		eta = x[n - 1][pom1];
		
		N[0] = 0.25 * (1 - ksi) * (1 - eta);
		N[1] = 0.25 * (1 + ksi) * (1 - eta);
		N[2] = 0.25 * (1 + ksi) * (1 + eta);
		N[3] = 0.25 * (1 - ksi) * (1 + eta);

		for (int i = 0; i < 4; i++)
		{
			int wynelem = int(id.ID[i]);
			result[0] += ksieta.eta[i][k] * nodes.Node[wynelem - 1].x;
			result[1] += ksieta.eta[i][k] * nodes.Node[wynelem - 1].y;
			result[2] += ksieta.ksi[i][k] * nodes.Node[wynelem - 1].x;
			result[3] += ksieta.ksi[i][k] * nodes.Node[wynelem - 1].y;
			
		}
		detj = result[0] * result[3] - result[1] * result[2];


		if (k < n) { wagi1 = 0; wagi2 = 0; }
		if (k < 2 * n && k >= n) { wagi1 = 1; wagi2 = n; }
		if (k < 3 * n && k >= 2 * n) { wagi1 = 2; wagi2 = n * 2; }
		if (k < 4 * n && k >= 3 * n) { wagi1 = 3; wagi2 = n * 3; }
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				C[i][j] = N[i] * N[j] * data.SpecificHeat*data.Density*detj * w[n - 1][wagi1] * w[n - 1][k - wagi2];
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				id.C[i][j] += C[i][j];
			}
		}

		for (int i = 0; i < 4; i++)
		{
			result[i] = 0;
		}

	}

}

void wys_C(Elements& id) {
	cout << "\nmacierz C";
	for (int i = 0; i < 4; i++) {
		cout << endl;
		for (int j = 0; j < 4; j++) {
			cout << id.C[i][j] << ", ";
		}
	}
}