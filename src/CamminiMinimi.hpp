//Cammini Minimi .hpp

#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <iomanip> // per setprecision
#include <cmath>
#include <limits>
#include <algorithm>
#include <Eigen/Dense>
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;
using ListaAdiacenza = vector<vector<pair<unsigned int, double>>>;
//using ListaAdiacenza = unordered_map<unsigned int, vector<pair<unsigned int, double>>>;

ListaAdiacenza CreaListaAdiacenza(const PolyhedraMesh& mesh);

pair<vector<unsigned int>, double> Dijkstra(const ListaAdiacenza& ListaAd, unsigned int nodo_iniziale, unsigned int nodo_finale);