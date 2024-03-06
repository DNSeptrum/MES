#include <cmath>
#include <iostream>
#include<iomanip>
using namespace std;
#include "wczytywanie.h"

struct ElementUniversalny {
public:
    long double ksi[4][16];
    long double eta[4][16];
}; 

struct MacierzH {
    double H[4][4];
};

double gauss_legendre_1D(int n, double (*F)(double));
double gauss_legendre_2D(int n, double (*F)(double, double));


double f(double x);
double f(double x, double y);

double n1ksi(double eta);
double n2ksi(double eta);
double n3ksi(double eta);
double n4ksi(double eta);

double n1eta(double ksi);
double n2eta(double ksi);
double n3eta(double ksi);
double n4eta(double ksi);

void gauss_legendre_2D_eta_ksi(int n, ElementUniversalny& wsk);

void jakobian(ElementUniversalny& ksieta, Elements& id, Grind& nodes, int n, Global_Data& data);

void inintmacierz(Elements& wypisz);

void wys_H(Elements& wypisz);

void wys_ksi(ElementUniversalny& ksieta,int n);

void wys_eta(ElementUniversalny& ksieta,int n);

void macierzHbc(ElementUniversalny& ksieta, Elements& id, Grind& nodes, Global_Data& data, int n);

void wys_Hbc(Elements& id);

void wys_P(Elements& id);

void FH(Elements& id);
void wys_FH(Elements& id,int k);

void FinalH(Grind& grind, Global_Data& data);

void wys_GlobalH(Grind& grind, Global_Data& data);

void macierzC(ElementUniversalny& ksieta, Elements& id, Grind& nodes, Global_Data& data, int n);

void wys_C(Elements& id);

void FinalP(Grind& grind, Global_Data& data);

void wys_GlobalP(Grind& grind, Global_Data& data);

void FinalC(Grind& grind, Global_Data& data);

void wys_GlobalC(Grind& grind, Global_Data& data);

void FinalHC(Grind& grind, Global_Data& data);

void wys_GlobalHC(Grind& grind, Global_Data& data);

void obliczanie(Grind& grind, Global_Data& data);

int odwracanie(int n, Grind& grind, Global_Data& data);

void min(Grind grind, Global_Data data);

void max(Grind grind, Global_Data data);

void initTemp(Grind& grind, Global_Data& data);

int zapis(int n,Grind& grind, Global_Data& data);