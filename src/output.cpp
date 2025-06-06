//output.cpp

#include <iostream>
#include <fstream>
#include <iomanip> // per setprecision
#include "Eigen/Eigen"
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;

// !! nel main: 
// PolyhedraMesh triangolata = TriangolazioneI(PolyhedraMesh& mesh, unsigned int b) funzione per triangolare
// outputFile0(triangolata, "Cell0Ds.txt") !!

// Output: funzione che stampi in output Cell0Ds.txt  
bool outputFile(const PolyhedraMesh& triangolata, const string& OutputFile0, const string& OutputFile1, const string& OutputFile2, const string& OutputFile3) 
//passo la nuova mesh, generata con PolyhedraMesh triangolata = TriangolazioneI(PolyhedraMesh& mesh, unsigned int b) per riferimento costante, tanto non devo modificarla
//e passo il nome del file di output come stringa, così da poterlo cambiare facilmente
{
    // file di output0
    ofstream Cell0Ds;
    
    // Apro i file di output
    Cell0Ds.open(OutputFile0);
    
    // Controllo se i file sono stati aperti correttamente
    if (!Cell0Ds.is_open())
    {
        cerr << "Error opening output file." << endl;
        return false;
    }

    // Scrivo i dati nei file
    Cell0Ds << "ID X Y Z" << endl; // Header per le coordinate dei punti 0D
    /*
    NumCell0Ds unsigned int numero di vertici ==> n righe 
    vector<unsigned int> Cell0DsId ==> ID dei vertici
    MatrixXd Cell0DsCoordinates ==> matrice 3xNumCell0Ds con le coordinate dei vertici
    */

    Cell0Ds << fixed << setprecision(16);
    
    // Ciclo per scrivere le coordinate dei punti 0D
    for (unsigned int i=0; i<triangolata.NumCell0Ds; i++) //n righe
    {
        Cell0Ds.unsetf(ios::fixed); // disattiva il formato fixed, perché gli id mi servono interi (potrebbe essere superfluo, perché gli id sono già interi)
        Cell0Ds << triangolata.Cell0DsId[i] << " ";//ID del vertice
        Cell0Ds << fixed << setprecision(16)
                << triangolata.Cell0DsCoordinates(0, i) << " " //coordinata x
                << triangolata.Cell0DsCoordinates(1, i) << " " //coordinata y
                << triangolata.Cell0DsCoordinates(2, i) << endl; //coordinata z
    }

    // Chiudo il file
    Cell0Ds.close();


    // file di output1
    ofstream Cell1Ds;
    
    // Apro i file di output
    Cell1Ds.open(OutputFile1);
    
    // Controllo se i file sono stati aperti correttamente
    if (!Cell1Ds.is_open())
    {
        cerr << "Error opening output file." << endl;
        return false;
    }

    // Scrivo i dati nei file
    Cell1Ds << "ID Origin End" << endl; // Header per le coordinate dei punti 1D
    /*
    NumCell1Ds unsigned int numero di edges ==> n righe 
    vector<unsigned int> Cell1DsId ==> ID degli edges
    MatrixXi Cell1DsExtrema ==> matrice 2 x NumCell1Ds con gli estremi degli edges 
    */

    //non ho bisogno di setprecision, dato che gli ID dei vertici sono interi

    // Ciclo per scrivere le coordinate degli edges
    for (unsigned int i=0; i<triangolata.NumCell1Ds; i++) //n righe
    {
        Cell1Ds << triangolata.Cell1DsId[i] << " " //ID dell'edge
                << triangolata.Cell1DsExtrema(0, i) << " " //ID del vertice di partenza
                << triangolata.Cell1DsExtrema(1, i) << endl; //ID del vertice di arrivo
    }

    // Chiudo il file
    Cell1Ds.close();


    // file di output2
    ofstream Cell2Ds;
    
    // Apro i file di output
    Cell2Ds.open(OutputFile2);
    
    // Controllo se i file sono stati aperti correttamente
    if (!Cell2Ds.is_open())
    {
        cerr << "Error opening output file." << endl;
        return false;
    }

    // Scrivo i dati nei file
    Cell2Ds << "ID NumVertices Vertices NumEdges Edges" << endl; // Header per le coordinate dei punti 2D
    /*
    NumCell2Ds unsigned int numero di facce ==> n righe 
    vector<unsigned int> Cell2DsId ==> ID delle facce
    
    vector<unsigned int> Cell2DsNumVertices ==> vettore con numero di vertici
    vector<unsigned int> Cell2DsNumEdges ==> vettore con numero di edges

    vector<vector<unsigned int>> Cell2DsVertices ==> vettore di vettori con gli ID dei vertici delle facce
    vector<vector<unsigned int>> Cell2DsEdges ==> vettore di vettori con gli ID degli edges delle facce 
    
    vector<unsigned int> Cell2DsNumVertices ==> numero di vertici di ogni faccia (poligono)
    vector<unsigned int> Cell2DsNumEdges ==> numero di spigoli di ogni faccia (poligono)
    */
    for (unsigned int i=0; i<triangolata.NumCell2Ds; i++) //n righe
    {
        Cell2Ds << triangolata.Cell2DsId[i] << " " //ID della faccia
                << triangolata.Cell2DsNumVertices(i) << " "; //numero di vertici della faccia (sempre 3)

        //for (auto vert in triangolata.Cell2DsVertices[i], vert < triangolata.NumCell2Ds, i++) 
        for (auto vert : triangolata.Cell2DsVertices[i]) 
            Cell2Ds << vert << " "; //scrivo gli ID dei vertici della faccia, che sono un vettore di ID       
        Cell2Ds << triangolata.Cell2DsNumEdges(i) << " "; //numero di edges della faccia (sempre 3)
        for (auto edg : triangolata.Cell2DsEdges[i]) //ID degli edges della faccia 
            Cell2Ds << edg << " "; 
        Cell2Ds << endl;
    }

    // Chiudo il file
    Cell2Ds.close();


    // file di output3
    ofstream Cell3Ds;
    
    // Apro i file di output
    Cell3Ds.open(OutputFile3);
    
    // Controllo se i file sono stati aperti correttamente
    if (!Cell3Ds.is_open())
    {
        cerr << "Error opening output file." << endl;
        return false;
    }

    // Scrivo i dati nei file
    Cell3Ds << "ID NumVertices Vertices NumEdges Edges NumFaces Faces" << endl; // Header per le coordinate dei punti 3D
   
    /*
    NumCell3Ds unsigned int numero di facce ==> n righe 
    unsigned int Cell3DsId ==> ID del poliedro (sempre 0)

    unsigned int Cell3DsNumVertices ==> numero di vertici
    unsigned int Cell3DsNumEdges ==> numero di edges
    unsigned int Cell3DsNumFaces ==> numero di facce

    vector<unsigned int> Cell3DsEdges  ==> vettore con gli ID dei vertici del poliedro
    vector<unsigned int> Cell3DsVertices ==> vettore con gli ID degli edges del poliedro
    vector<unsigned int> Cell3DsFaces ==> vettore con gli ID delle facce del poliedro
    */
    Cell3Ds << triangolata.Cell3DsId << " " //ID del poliedro
            << triangolata.Cell3DsNumVertices << " "; //numero di vertici del poliedro
    for(auto v : triangolata.Cell3DsVertices) //Id dei vertici
        Cell3Ds << v << " ";
    Cell3Ds << triangolata.Cell3DsNumEdges << " "; //numero di edges del poliedro
    for(auto v : triangolata.Cell3DsEdges) //Id degli edges
        Cell3Ds << v << " ";//Id degli edges
    Cell3Ds << triangolata.Cell3DsNumFaces << " "; //numero di facce del poliedro
    for(auto v : triangolata.Cell3DsFaces) //Id delle facce
        Cell3Ds << v << " ";//Id delle facce
    Cell3Ds << endl; //fine della riga (superfluo sec. me)

    // Chiudo il file
    Cell3Ds.close();

    return true;
}





