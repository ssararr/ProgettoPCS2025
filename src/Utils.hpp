//POLYHEDRA UTILS HPP

#pragma once

#include <iostream>
#include "PolyhedraMesh.hpp"

using namespace std;

<<<<<<< Updated upstream
namespace PolyhedraLibrary
{

/* Importa la mesh poligonale e verifica se la mesh è corretta, ritorna true se la lettura è andata a buon fine, false altrimenti */
bool ImportMesh(PolyhedraMesh& mesh);

/* Importa le proprietà delle celle 0D dal file Cell0Ds.csv ritorna il risultato della lettura, true se ha successo, false altrimenti */
bool ImportCell0Ds(PolyhedraMesh& mesh);

/* Importa le proprietà delle celle 1D dal file Cell1Ds.csv ritorna il risultato della lettura, true se ha successo, false altrimenti */
bool ImportCell1Ds(PolyhedraMesh& mesh);

/* Importa le proprietà delle celle 2D dal file Cell2Ds.csv ritorna il risultato della lettura, true se ha successo, false altrimenti */
bool ImportCell2Ds(PolyhedraMesh& mesh);

/* Importa le proprietà delle celle 3D dal file Cell3Ds.csv ritorna il risultato della lettura, true se ha successo, false altrimenti */
=======
namespace PolygonalLibrary
{

// Funzione che prende una mesh 
// e la riempie leggendo da file, se non riesce restituisce falso

bool ImportMesh(PolyhedraMesh& mesh);


bool ImportCell0Ds(PolyhedraMesh& mesh);


bool ImportCell1Ds(PolyhedraMesh& mesh);


bool ImportCell2Ds(PolyhedraMesh& mesh);


>>>>>>> Stashed changes
bool ImportCell3Ds(PolyhedraMesh& mesh);

}
