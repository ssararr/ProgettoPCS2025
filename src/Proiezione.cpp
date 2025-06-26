//Proiezione.cpp

#include <iostream>
#include <fstream>
#include <iomanip> // per setprecision
#include <cmath>
#include "Eigen/Eigen"
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;

//Si accede alle coordinate di ogni vertice Cell0DsCoordinates e si normalizzano
void ProiezioneSfera(PolyhedraMesh& meshtriangolata) 
{
	for (unsigned int i=0; i < meshtriangolata.NumCell0Ds; i++)
    {
        //precisione di 16 cifre decimali
        cout << setprecision(16);
		//coordinate del vertice i
		double x = meshtriangolata.Cell0DsCoordinates(0,i);
		double y = meshtriangolata.Cell0DsCoordinates(1,i);
		double z = meshtriangolata.Cell0DsCoordinates(2,i);
		
		//calcolo la norma
		double norma = sqrt(x * x + y * y + z * z);
		
		//modifico le coordinate
		meshtriangolata.Cell0DsCoordinates(0,i)= meshtriangolata.Cell0DsCoordinates(0,i)/norma;
		meshtriangolata.Cell0DsCoordinates(1,i)= meshtriangolata.Cell0DsCoordinates(1,i)/norma;
		meshtriangolata.Cell0DsCoordinates(2,i)= meshtriangolata.Cell0DsCoordinates(2,i)/norma;
		
	}
}