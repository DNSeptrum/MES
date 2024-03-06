#include "obliczenia.h"

void initTemp(Grind& grind, Global_Data& data) {
	
		grind.Temp0 = new double[data.Nodes_number];
		for (int i = 0; i < data.Nodes_number; i++) {
			grind.Temp0[i] = data.InitialTemp;
		}
	
}

void FinalP(Grind& grind, Global_Data& data) {

	grind.FinalP = new double [data.Nodes_number];

	for (int i = 0; i < data.Nodes_number; i++) {
			grind.FinalP[i] = 0;
	}

	int pom1 = 0;
	int pom2 = 0;
	for (int k = 0; k < data.Elements_number; k++) {
		for (int i = 0; i < 4; i++) {
			pom1 = int(grind.Element[k].ID[i]);

				grind.FinalP[pom1 - 1] += grind.Element[k].P[i];
		}

		}
	for (int k = 0; k < data.Nodes_number; k++) {
		for (int i = 0; i < data.Nodes_number; i++) {
			grind.FinalP[k] += grind.FinalC[k][i] / data.SimulationStepTime * grind.Temp0[i];
		}
	}
	
}

void wys_GlobalP(Grind& grind, Global_Data& data) {
	cout << "\nGlobal_P";
	for (int i = 0; i < data.Nodes_number; i++) {
			cout << grind.FinalP[i] << " , ";
	}
}