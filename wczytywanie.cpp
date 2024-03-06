#include "wczytywanie.h"
using namespace std;

void wys_Node(Grind& grind, Global_Data& data) {
	cout << "Nodes:" << endl;
	for (int i =0; i< data.Nodes_number; i++)
	cout << "Dla Node["<< i << "] " << "x = " << grind.Node[i].x <<" y = "<< grind.Node[i].y << " BC= " << grind.Node[i].Bc << endl;
}

void wys_Element(Grind& grind, Global_Data& data) {
	cout << "Elements:" << endl;
	for (int i = 0; i < data.Elements_number; i++)
		cout << "Dla Element[" << i << "] wartosci wynosza:" << grind.Element[i].ID[0] << ","<< grind.Element[i].ID[1] << ',' << grind.Element[i].ID[2] << "," << grind.Element[i].ID[3] << endl;
}

void wys_data(Global_Data& data) {
	cout << "Global data:" << endl;
	cout <<"SimulationTime: " << data.SimulationTime << endl;
	cout <<"SimulationStepTime: " << data.SimulationStepTime << endl;
	cout << "Conductivity: " << data.Conductivity << endl;
	cout << "Alfa: " << data.Alfa << endl;
	cout << "Tot: " << data.Tot << endl;
	cout << "InitialTemp: " << data.InitialTemp << endl;
	cout << "Density: " << data.Density << endl;
	cout << "SpecificHeat: " << data.SpecificHeat << endl;
	cout << "Nodes_number: " << data.Nodes_number << endl;
	cout << "Elements_number: " << data.Elements_number << endl;
}

void wczytajZPliku(const std::string& nazwaPliku, Global_Data& data,Grind& grind)
{
	std::ifstream plik(nazwaPliku);

	if (!plik.is_open())
	{
		std::cerr << "Nie mozna otworzyc pliku." << std::endl;
		return;
	}

	std::string token;
	while (plik >> token)
	{
		if (token == "SimulationTime") {
			plik >> data.SimulationTime;
		}
		else if (token == "SimulationStepTime") {
			plik >> data.SimulationStepTime;
		}
		else if (token == "Conductivity") {
			plik >> data.Conductivity;
		}
		else if (token == "Alfa") {
			plik >> data.Alfa;
		}
		else if (token == "Tot") {
			plik >> data.Tot;
		}
		else if (token == "InitialTemp") {
			plik >> data.InitialTemp;
		}
		else if (token == "Density") {
			plik >> data.Density;
		}
		else if (token == "SpecificHeat") {
			plik >> data.SpecificHeat;
		}
		else if (token == "Nodes" && plik >> token && token == "number") {
			plik >> data.Nodes_number;
		}
		else if (token == "Elements" && plik >> token && token == "number") {
			plik >> data.Elements_number;
		}
		
		else if (token == "*Node") 
		{
			grind.Node = new Nodes[data.Nodes_number];
			string pom;
			int nodeIndex = 0;
			while (nodeIndex <= data.Nodes_number)
			{
				grind.Node[nodeIndex].Bc = 0;
				plik >> pom >> grind.Node[nodeIndex].x >> pom >> grind.Node[nodeIndex].y; //>> grind.Node[nodeIndex].y;
				++nodeIndex;
				if (nodeIndex >= data.Nodes_number)
				{
					break;
				}
			}
		}
		
		else if (token == "*Element,")
		{
			string pom;
			plik >> pom; 
			int elementIndex = 0;
			grind.Element = new Elements[data.Elements_number];
			while (plik >> token)
			{
				
				for (int i = 0; i < 4; ++i)
				{
					if (i != 3)	plik >> grind.Element[elementIndex].ID[i] >> pom;
					else plik >> grind.Element[elementIndex].ID[i];
				}
				++elementIndex;
				if (elementIndex >= data.Elements_number)
				{
					break;
				}
			}
		}
		else if (token == "*BC") {
			string pom;
				int pom2 =0 ;
			int nodeIndex = 0;
			while (nodeIndex <= data.Nodes_number)
			{
				plik >> pom2 >> pom; 
				++nodeIndex;
			
				if (pom2 >= 0 && pom2 < data.Nodes_number+1) {
					grind.Node[pom2-1].Bc = pom2;
				}
				if (nodeIndex >= data.Nodes_number || pom2 == NULL)
				{
					break;
				}
			}
		}
	}
	plik.close();
}