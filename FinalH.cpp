#include "obliczenia.h"

void FH(Elements& id) {
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			id.FH[i][j] = 0;
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			id.FH[i][j] = id.H[i][j]+id.Hbc[i][j];
		}
	}
	
}

void wys_FH(Elements& id,int k) {
	cout << "\nFinal_H[" << k << "] ";
	for (int i = 0; i < 4; i++) {
		cout << endl;
		for (int j = 0; j < 4; j++) {
			cout << id.FH[i][j] << " , ";
		}
	}
}

void FinalH(Grind& grind, Global_Data& data) {
	
	//inicjalizacja Global_H
	grind.FinalH = new double*[data.Nodes_number];
	for (int i = 0; i < data.Nodes_number; i++) {
		grind.FinalH[i] = new double[data.Nodes_number];
	}

	for (int i = 0; i < data.Nodes_number; i++) {
		for (int j = 0; j < data.Nodes_number; j++)
			grind.FinalH[i][j] = 0;
	}

	int pom1 = 0;
	int pom2 = 0;
	for (int k = 0; k<data.Elements_number; k++) {	
		for (int i = 0; i < 4; i++) {
			pom1 = int(grind.Element[k].ID[i]);
			for (int j = 0; j < 4; j++) {
				
				pom2 = int(grind.Element[k].ID[j]);
				grind.FinalH[pom1-1][pom2-1] += grind.Element[k].FH[i][j];
			}

		}
	}
}

void wys_GlobalH(Grind& grind, Global_Data& data) {
	cout << "\nGlobal_H";
	for (int i = 0; i < data.Nodes_number; i++) {
		cout << endl;
		for (int j = 0; j < data.Nodes_number; j++) {
			cout << grind.FinalH[i][j] << " , ";
		}
	}
}

void FinalC(Grind& grind, Global_Data& data) {

	//inicjalizacja Global_H
	grind.FinalC = new double* [data.Nodes_number];
	for (int i = 0; i < data.Nodes_number; i++) {
		grind.FinalC[i] = new double[data.Nodes_number];
	}

	for (int i = 0; i < data.Nodes_number; i++) {
		for (int j = 0; j < data.Nodes_number; j++)
			grind.FinalC[i][j] = 0;
	}

	int pom1 = 0;
	int pom2 = 0;
	for (int k = 0; k < data.Elements_number; k++) {
		for (int i = 0; i < 4; i++) {
			pom1 = int(grind.Element[k].ID[i]);
			for (int j = 0; j < 4; j++) {

				pom2 = int(grind.Element[k].ID[j]);
				grind.FinalC[pom1 - 1][pom2 - 1] += grind.Element[k].C[i][j];
			}

		}
	}
}

void wys_GlobalC(Grind& grind, Global_Data& data) {
	cout << "\nGlobal_C";
	for (int i = 0; i < data.Nodes_number; i++) {
		cout << endl;
		for (int j = 0; j < data.Nodes_number; j++) {
			cout << grind.FinalC[i][j] << " , ";
		}
	}
}

void FinalHC(Grind& grind, Global_Data& data) {

	//inicjalizacja Global_H
	grind.FinalHC = new double* [data.Nodes_number];
	for (int i = 0; i < data.Nodes_number; i++) {
		grind.FinalHC[i] = new double[data.Nodes_number];
	}

	for (int i = 0; i < data.Nodes_number; i++) {
		for (int j = 0; j < data.Nodes_number; j++)
			grind.FinalHC[i][j] = 0;
	}

		for (int i = 0; i < data.Nodes_number; i++) {

			for (int j = 0; j < data.Nodes_number; j++) {

				grind.FinalHC[i][j] += grind.FinalH[i][j]+(grind.FinalC[i][j]/data.SimulationStepTime);
			}

		}
}

void wys_GlobalHC(Grind& grind, Global_Data& data) {
	cout << "\nGlobal_H+C/dt";
	for (int i = 0; i < data.Nodes_number; i++) {
		cout << endl;
		for (int j = 0; j < data.Nodes_number; j++) {
			cout << grind.FinalHC[i][j] << " , ";
		}
	}
}