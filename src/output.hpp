//output.hpp

#pragma once
#include <string> 
#include "PolyhedraMesh.hpp" 

using namespace std;
using namespace Eigen;

bool outputFile(
    const PolyhedraMesh& triangolata,
    const string& outputFile0,
    const string& outputFile1,
    const string& outputFile2,
    const string& outputFile3
);
