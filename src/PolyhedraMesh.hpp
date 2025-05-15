//POLYHEDRA MESH.HPP
#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedraLibrary {

    struct PolyhedraMesh {

        unsigned int NumCell0Ds;
        unsigned int NumCell1Ds;
        unsigned int NumCell2Ds;
        unsigned int NumCell3Ds;

        // Vettori contenenti Id, ovvero numeri interi
        vector<unsigned int> Cell0DsId;
        vector<unsigned int> Cell1DsId;
        vector<unsigned int> Cell2DsId;
        unsigned int Cell3DsId;

        MatrixXd Cell0DsCoordinates;    // Matrice di coordinate dei punti, contiene doubles
        MatrixXi Cell1DsExtrema;    // Matrice di estremi dei segmenti, nel file sono tutti int

        // Vettore esterno contiene i poligoni, vettori interni contengono vertici dei poligoni 
        // (hanno dim diversa perché i poligoni non sono necessariamente uguali)
        vector<vector<unsigned int>> Cell2DsVertices; 
        vector<vector<unsigned int>> Cell2DsEdges;
        vector<unsigned int> Cell2DsNumVertices; //Salvo il numero di vertici di ogni faccia (poligono)
        vector<unsigned int> Cell2DsNumEdges; //Salvo il numero di spigoli di ogni faccia (poligono)

        vector<unsigned int> Cell3DsFaces; //Salvo gli ID delle facce (poligoni) che compongono il poliedro
        vector<unsigned int> Cell3DsEdges; //Salvo gli ID degli spigoli che compongono il poliedro
        vector<unsigned int> Cell3DsVertices; //Salvo gli ID dei vertici che compongono il poliedro
        unsigned int Cell3DsNumFaces; 
        unsigned int Cell3DsNumEdges;
        unsigned int Cell3DsNumVertices;                
    };
}

