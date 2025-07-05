//POLYHEDRA UTILS HPP

#pragma once
#include <iostream>
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;


void OrderFaces(const vector<unsigned int>& unordered_faces, 
               vector<unsigned int>& ordered_faces, 
               const PolyhedraMesh& mesh); //funzione che ordina le facce di una mesh in modo che siano adiacenti tra loro, partendo da una faccia iniziale

namespace PolyhedraLibrary {             //riempie la nuova mesh con i vertici duali, gli spigoli duali e le facce duali

bool MeshDuale(PolyhedraMesh& mesh);

}
