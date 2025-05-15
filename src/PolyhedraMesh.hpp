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

        // Vettori di ID
        vector<unsigned int> Cell0DsId;
        vector<unsigned int> Cell1DsId;
        vector<unsigned int> Cell2DsId;
        unsigned int Cell3DsId;

        MatrixXd Cell0DsCoordinates;    // Matrice di coordinate dei vertici, contiene doubles
        MatrixXi Cell1DsExtrema;    // Matrice contenetnte id di estremi dei segmenti

        // Vettori con Numero di vertici/edges/facce
        vector<unsigned int> Cell2DsNumVertices;
        vector<unsigned int> Cell2DsNumEdges;
        unsigned int Cell3DsNumVertices;
        unsigned int Cell3DsNumEdges;
        unsigned int Cell3DsNumFaces;



        // Salviamo gli ID dei vertici, dei lati costituenti le facce
        vector<vector<unsigned int>> Cell2DsVertices; 
        vector<vector<unsigned int>> Cell2DsEdges;
        vector<unsigned int> Cell2DsNumVertices; //Salvo il numero di vertici di ogni faccia (poligono)
        vector<unsigned int> Cell2DsNumEdges; //Salvo il numero di spigoli di ogni faccia (poligono)

              
        // Vettori contenenti ID dei vertici, archi e facce del poliedro
        vector<unsigned int> Cell3DsEdges;
        vector<unsigned int> Cell3DsVertices;
        vector<unsigned int> Cell3DsFaces;

    };
}

