// Triangolazione.hpp

#pragma once
#include <iostream>
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedraLibrary{

    PolyhedraMesh TriangolazioneI(PolyhedraMesh& mesh, unsigned int b);

}
