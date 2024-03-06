#pragma once
#include <iostream>
#include <cstring>
#include <fstream>
using namespace std;

struct Global_Data {

	int SimulationTime;
	int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
	int Nodes_number;
	int Elements_number;
};

struct Nodes {
	double x;
	double y;
	int Bc;
};

struct Elements {
	double ID[4];
	double H[4][4];
	double Hbc[4][4];
	double P[4];
	double FH[4][4];
	double C[4][4];
};

struct Grind {
	Nodes* Node;
	Elements* Element;
	double** FinalH;
	double** FinalC;
	double* FinalP;
	double** FinalHC;
	double* Temp;
	double* Temp0;
};

void wczytajZPliku(const std::string& nazwaPliku, Global_Data& data, Grind& grind);
void wys_Node(Grind& grind, Global_Data& data);
void wys_Element(Grind& grind, Global_Data& data);
void wys_data(Global_Data& data);