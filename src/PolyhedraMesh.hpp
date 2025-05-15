#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedraLibrary {

    struct PolyhedronMesh {

        unsigned int NumCell0Ds;
        unsigned int NumCell1Ds;
        unsigned int NumCell2Ds;

        // Vettori contenenti Id, ovvero numeri interi
        vector<unsigned int> Cell0DsId;
        vector<unsigned int> Cell1DsId;
        vector<unsigned int> Cell2DsId;

        MatrixXd Cell0DsCoordinates;    // Matrice di coordinate dei punti, contiene doubles
        MatrixXi Cell1DsExtrema;    // Matrice di estremi dei segmenti, nel file sono tutti int

        // Vettore esterno contiene i poligoni, vettori interni contengono vertici dei poligoni 
        // (hanno dim diversa perché i poligoni non sono necessariamente uguali)
        vector<vector<unsigned int>> Cell2DsVertices; 
        vector<vector<unsigned int>> Cell2DsEdges;

        // Uso i dizionari per associare i markers a una lista a cui poi aggiugerò i dati
        map<unsigned int, list<unsigned int>> Cell0Ds_markers;
        map<unsigned int, list<unsigned int>> Cell1Ds_markers;
        map<unsigned int, list<unsigned int>> Cell2Ds_markers;
        
    };
}

