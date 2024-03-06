#include "obliczenia.h"
	
void macierzHbc(ElementUniversalny& ksieta, Elements& id, Grind& nodes, Global_Data& data, int n) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {

			id.Hbc[i][j] = 0.0;
		}
	}
	id.P[0] = 0;
	id.P[1] = 0;
	id.P[2] = 0;
	id.P[3] = 0;
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
	int boki[4] = { 0,0,0,0 };
	// sprawdzanie boków
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < data.Nodes_number; i++) {
			if (id.ID[j] == nodes.Node[i].Bc) { boki[j] = 1; }
		}
	}

	double ksi = 0, eta = 0;
	double N[4];
	double WektorP[4];
	double H[4][4];
	double detj = 0;
	double idx[4] = { 0.0,0.025,0.025,0.0 }, idy[4] = { 0.0,0.0,0.025,0.025 };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] = 0;
		}
	}
	for (int j = 0; j < 4; j++) {
		WektorP[j] = 0;
	}
	for (int k = 0; k < 4; k++) {
		for (int s = 0; s < n; s++) {

			// jakobian
			if (k == 0) {
				ksi = x[n - 1][s];
				eta = -1;
				detj = sqrt(pow((nodes.Node[int(id.ID[0] - 1)].x - nodes.Node[int(id.ID[1]) - 1].x), 2) + pow((nodes.Node[int(id.ID[0]) - 1].y - nodes.Node[int(id.ID[1]) - 1].y), 2)) / 2;
			}
			if (k == 1) {
				ksi = 1;
				eta = x[n - 1][s];
				detj = sqrt(pow((nodes.Node[int(id.ID[1] - 1)].x - nodes.Node[int(id.ID[2]) - 1].x), 2) + pow((nodes.Node[int(id.ID[1]) - 1].y - nodes.Node[int(id.ID[2]) - 1].y), 2)) / 2;
			}
			if (k == 2) {
				ksi = x[n - 1][s];
				eta = 1;
				detj = sqrt(pow((nodes.Node[int(id.ID[2] - 1)].x - nodes.Node[int(id.ID[3]) - 1].x), 2) + pow((nodes.Node[int(id.ID[2]) - 1].y - nodes.Node[int(id.ID[3]) - 1].y), 2)) / 2;
			}

			if (k == 3) {
				ksi = -1;
				eta = x[n - 1][s];
				detj = sqrt(pow((nodes.Node[int(id.ID[3] - 1)].x - nodes.Node[int(id.ID[0]) - 1].x), 2) + pow((nodes.Node[int(id.ID[3]) - 1].y - nodes.Node[int(id.ID[0]) - 1].y), 2)) / 2;
			}
			N[0] = 0.25 * (1 - ksi) * (1 - eta);
			N[1] = 0.25 * (1 + ksi) * (1 - eta);
			N[2] = 0.25 * (1 + ksi) * (1 + eta);
			N[3] = 0.25 * (1 - ksi) * (1 + eta);

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					H[i][j] += N[i] * N[j] * w[n - 1][s] * data.Alfa * detj;
				}
			}
			for (int i = 0; i < 4; i++) {
				WektorP[i] += N[i] * w[n - 1][s] * data.Tot * data.Alfa * detj;
			}

		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {



				if (boki[0] == 1 && boki[1] == 1 && k == 0) {
					id.Hbc[i][j] += H[i][j];
				}
				if (boki[1] == 1 && boki[2] == 1 && k == 1) {
					id.Hbc[i][j] += H[i][j];
				}

				if (boki[2] == 1 && boki[3] == 1 && k == 2) {
					id.Hbc[i][j] += H[i][j];
				}
				if (boki[3] == 1 && boki[0] == 1 && k == 3) {
					id.Hbc[i][j] += H[i][j];
				}
			}
		}

		for (int j = 0; j < 4; j++)
		{



			if (boki[0] == 1 && boki[1] == 1 && k == 0) {
				id.P[j] += WektorP[j];
			}
			if (boki[1] == 1 && boki[2] == 1 && k == 1) {
				id.P[j] += WektorP[j];
			}

			if (boki[2] == 1 && boki[3] == 1 && k == 2) {
				id.P[j] += WektorP[j];
			}
			if (boki[3] == 1 && boki[0] == 1 && k == 3) {
				id.P[j] += WektorP[j];
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				H[i][j] = 0;
			}
		}
		for (int j = 0; j < 4; j++) {
			WektorP[j] = 0;
		}

	}
}

void wys_Hbc(Elements& id) {
	cout << "\nmacierz Hbc";
	for (int i = 0; i < 4; i++) {
		cout << endl;
		for (int j = 0; j < 4; j++) {
			cout << id.Hbc[i][j] << ", ";
		}
	}
}
void wys_P(Elements& id) {
	cout << "\nwektor P";
		cout << endl;
		for (int i = 0; i < 4; i++) {
			cout << id.P[i] << ", ";
		}
	}