//POLYHEDRA UTILS HPP

#pragma once
#include <iostream>
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;

double distanceBtw(const vector<double>& point1, const vector<double>& point2); //funzione che conta la distanza tra due vettori di coordinate	3D

vector<vector<double>> eigenMatrixToVectorVector(const MatrixXd& matrix); //funzione che converte un oggetto MatrixXd in un vettore di vettori di double

map<int, vector<int>> findClosestDualVertices(       //funzione che, dati due vector vector contenenti le coordinate dei vertici originali e dei vertici duali, restituisce una mappa in cui associa ad ogni id di un vertice originale gli id dei vertici duali più vicini
	const vector<vector<double>>& originalVertices,
	const vector<vector<double>>& dualVertices,
	double tolerance = 1e-6
);

namespace PolyhedraLibrary {             //riempie la nuova mesh con i vertici duali, gli spigoli duali e le facce duali

bool MeshDuale(PolyhedraMesh& mesh);

}
