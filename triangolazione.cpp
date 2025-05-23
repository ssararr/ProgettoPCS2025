
#include <iostream>
//#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;
#include "PolyhedraMesh.hpp"
// #include "Utils.hpp"
// #include "UCDUtilities.hpp"


// Funzione di triangolazione del poliedro tipo I
// generiamo uno spazio convesso tramite la combinazione convessa aA+bB+cC, con a+b+c=1, A,B,C vertici del triangolo della faccia
// la combinazione convessa viene discretizzata per generare (b+1)(b+2)/2 punti, dividendo la combinazione per b: (i/b)*A+(j/b)*B+(b-i-j/b)*C, i+j<=b
// i nuovi vertici sono (b+1)(b+2)/2
// i nuovi lati sono 3b(b+1)/2
// le nuove facce sono b^2

PolyhedraMesh TriangolazioneI(PolyhedraMesh& mesh, unsigned int b)
{
    NewNumVertices = (b + 1) * (b + 2) / 2; //numero di vertici
    NewNumEdges = 3 * b * (b + 1) / 2; //numero di edges

    MatrixXd newCell0DsCoordinates = MatrixXd::Zero(3, NewNumVertices); 
    MatrixXi newCell1DsExtrema = MatrixXi::Zero(2, NewNumEdges); 










    //Copilot
    PolyhedraMesh triangolo;
    triangolo.NumCell0Ds = (b + 1) * (b + 2) / 2; //numero di vertici
    triangolo.NumCell1Ds = 3 * b * (b + 1) / 2; //numero di edges
    triangolo.NumCell2Ds = b * b; //numero di facce

    triangolo.Cell0DsId.reserve(triangolo.NumCell0Ds); //riservo spazio per gli Id, cioè riservo spazio per il numero di vertici (punti 0D)
    triangolo.Cell1DsId.reserve(triangolo.NumCell1Ds); //riservo spazio per gli Id, cioè riservo spazio per il numero di spigoli (linee)
    triangolo.Cell2DsId.reserve(triangolo.NumCell2Ds); //riservo spazio per gli Id, cioè riservo spazio per il numero di facce (poligoni)

    triangolo.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, triangolo.NumCell0Ds); //riservo spazio per le coordinate x,y,z in una matrice 4xvertici
    triangolo.Cell1DsExtrema = Eigen::MatrixXi(2, triangolo.NumCell1Ds); //matrice di dimensione 2 x numero di spigoli per salvare gli estremi 
    triangolo.Cell2DsVertices.resize(triangolo.NumCell2Ds); //vettore di dimensione numero di facce per salvare i vertici delle facce
    triangolo.Cell2DsEdges.resize(triangolo.NumCell2Ds); //vettore di dimensione numero di facce per salvare gli spigoli delle facce

    unsigned int id = 0;
    for (unsigned int i = 0; i <= b; i++)
        for (unsigned int j = 0; j <= b - i; j++)
        {
            unsigned int k = b - i - j;
            triangolo.Cell0DsCoordinates(0, id) = i / static_cast<double>(b);
            triangolo.Cell0DsCoordinates(1, id) = j / static_cast<double>(b);
            triangolo.Cell0DsCoordinates(2,