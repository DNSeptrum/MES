#include "obliczenia.h"

void jakobian(ElementUniversalny& ksieta, Elements& id, Grind& nodes,int n,Global_Data& data)
{
	inintmacierz(id);
	double w[5][5] =
	{
		{2},
		{1,1},
		{5.0 / 9.0,8.0 / 9.0,5.0 / 9.0},
		{(18 + sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 - sqrt(30)) / 36,(18 + sqrt(30)) / 36},
		{0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851}
	};
	int wagi1 = 0, wagi2 = 0;
	ElementUniversalny pom, wynikczas;
	double detj;
	double result[4] = { 0,0,0,0 };
	// k = liczba punktów
	for (int k = 0; k < n * n; k++)
	{

		// obliczanie jakobianu
		for (int i = 0; i < 4; i++) 
		{
			int wynelem = int(id.ID[i]);
		result[0] += ksieta.eta[i][k] * nodes.Node[wynelem-1].x;
		result[1] += -ksieta.eta[i][k] * nodes.Node[wynelem-1].y;
		result[2] += -ksieta.ksi[i][k] * nodes.Node[wynelem-1].x;
		result[3] += ksieta.ksi[i][k] * nodes.Node[wynelem -1].y;
		
		}
	detj = result[0] * result[3] - result[1] * result[2];

		for (int i = 0; i < 4; i++) {
			result[i] = result[i] * 1 / detj;
		}

		for (int i = 0; i < 4; i++) 
		{
				wynikczas.eta[i][k] = result[3] * ksieta.eta[i][k] + (result[1] * ksieta.ksi[i][k]);
				wynikczas.ksi[i][k]= (result[2] * ksieta.eta[i][k]) + result[0] * ksieta.ksi[i][k];
		}

		for (int i = 0; i < 4; i++) {
			result[i] = 0;	
		}
	if (k < n) { wagi1 = 0; wagi2 = 0; }
	if (k < 2 * n && k >= n) { wagi1 = 1; wagi2 = n; }
	if (k < 3 * n && k >= 2 * n) { wagi1 = 2; wagi2 = n * 2; }
	if (k < 4 * n && k >= 3 * n) { wagi1 = 3; wagi2 = n * 3; }
		// obliczanie macierzy
		for (int i = 0; i < 4; i++) 
		{
			for (int j = 0; j < 4; j++) 
			{
				pom.eta[i][j] = wynikczas.eta[i][k] * wynikczas.eta[j][k];
				pom.ksi[i][j] = wynikczas.ksi[i][k] * wynikczas.ksi[j][k];
				id.H[i][j] += w[n - 1][wagi1] * w[n - 1][k - wagi2] * data.Conductivity * detj * (pom.eta[i][j] + pom.ksi[i][j]);
			}
			
		}
			
	}

}

void wys_H(Elements& wypisz) {
	cout << "\nwartosc macierzy H:";
	for (int i = 0; i < 4; i++) {
			cout << endl;
		for (int j = 0; j < 4; j++)
		{
			cout << wypisz.H[i][j] << ", ";
		}
	}
}

void inintmacierz(Elements& wypisz) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {

			wypisz.H[i][j] = 0.0;
		}
	}
}