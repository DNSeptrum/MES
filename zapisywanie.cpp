#include "obliczenia.h"
using namespace std;

#include <fstream>
#include <sstream>
#include <filesystem>

int zapis(int n,Grind& grid, Global_Data& data) {

        std::string nazwaPliku = "mesh" + std::to_string(n) + ".vtk";
        std::string sciezkaPliku = "paraview/" + nazwaPliku;

        std::ofstream plik(sciezkaPliku);

        if (plik.is_open()) {

            std::string naglowek = "# vtk DataFile Version 2.0\nUnstructured Grid Example\nASCII\nDATASET UNSTRUCTURED_GRID\n";

            plik << naglowek;


            plik << "\nPOINTS " << data.Nodes_number << " float\n";

            for (int j = 0; j < data.Nodes_number; ++j) {
                plik << grid.Node[j].x << " " << grid.Node[j].y << " 0\n";
            }
            //Cells
            plik << "\nCELLS " << data.Elements_number << " " << data.Elements_number * 5 << "\n";

            for (int j = 0; j < data.Elements_number; ++j) {
                plik << "4";
                for (int i = 0; i < 4; i++) {
                    plik << " " << grid.Element[j].ID[i]-1 ;
                }
                
                plik << "\n";
            }
            //CELLS_TYPE
            plik << "\nCELL_TYPES" << " " << data.Elements_number << "\n";
            for (int j = 0; j < data.Elements_number; ++j) {
                plik << "9\n";
            }
            //TEMP


            //    POINT_DATA 16
            //    SCALARS Temp float 1
            //    LOOKUP_TABLE default
            plik << "\nPOINT_DATA " << data.Nodes_number <<"\nSCALARS Temp float 1\nLOOKUP_TABLE default\n";
            for (int j = 0; j < data.Nodes_number;j++) {
                plik << grid.Temp0[j] << "\n";
            }

            plik.close();

            std::cout << "Utworzono plik: " << sciezkaPliku << std::endl;
        }
        else {
            std::cerr << "Nie udało się utworzyć pliku: " << sciezkaPliku << std::endl;
        }

    return 0;
}