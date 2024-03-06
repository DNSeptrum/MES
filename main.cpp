#include "obliczenia.h"

int main() {

    	Global_Data Gdata;
    	Grind grind;

    	//wczytajZPliku("Test1.txt", Gdata, grind);
          wczytajZPliku("Test2MixGrid.txt", Gdata, grind);
        //wczytajZPliku("Test3_31_31_kwadrat.txt", Gdata, grind);
        //wczytajZPliku("Test4_31_31_trapez.txt", Gdata, grind);

    	//wys_data(Gdata);
    	//wys_Node(grind, Gdata);
    	//wys_Element(grind, Gdata);

        int n =2;
        ElementUniversalny ksieta;

        gauss_legendre_2D_eta_ksi(n,ksieta);
     // wys_eta(ksieta,n);
     // wys_ksi(ksieta, n);

        for (int i = 0; i < Gdata.Elements_number; i++) {
         
            jakobian(ksieta, grind.Element[i], grind, n, Gdata);
            macierzHbc(ksieta, grind.Element[i], grind, Gdata, n);
            FH(grind.Element[i]);
            macierzC(ksieta, grind.Element[i], grind, Gdata, n);
         // wys_P(grind.Element[i]);
         // wys_H(grind.Element[i]);
         // wys_Hbc(grind.Element[i]);
         // wys_FH(grind.Element[i], i);     
         // wys_C(grind.Element[i]);
        }

        initTemp(grind, Gdata);      
        FinalH(grind, Gdata);     
        FinalC(grind, Gdata);
        FinalHC(grind, Gdata);
   
        odwracanie(n, grind,Gdata);
        
        int pom = Gdata.SimulationTime / Gdata.SimulationStepTime;
        for (int i = 0; i < pom; i++) {
            FinalP(grind, Gdata);
            obliczanie(grind, Gdata);
            min(grind, Gdata);
            max(grind, Gdata);
          //  zapis(i, grind, Gdata);
        }
        
     // wys_GlobalHC(grind, Gdata);
     // wys_GlobalH(grind, Gdata);
     // wys_GlobalP(grind, Gdata);
     // wys_GlobalC(grind, Gdata);

        return 0;  
}