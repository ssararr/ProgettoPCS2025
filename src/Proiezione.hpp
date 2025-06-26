//output.hpp

#pragma once
#include <iostream>
#include <fstream>
#include <iomanip> // per setprecision
#include <cmath>
#include "Eigen/Eigen"
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;

bool ProiezioneSfera(PolyhedraMesh& meshtriangolata);