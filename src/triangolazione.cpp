// TRIANGOLAZIONE.CPP

#include <iostream>
#include "Eigen/Eigen"
#include "PolyhedraMesh.hpp"
#include "Utils.hpp"
#include <map> // dizionari
#include <algorithm>  // per min e max
// #include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;


// Funzione per debuggare
void printMatrix(const Eigen::MatrixXi& m) {
    std::cout << "Matrix:\n" << m << std::endl;
}



// Funzione di triangolazione del poliedro tipo I
// generiamo uno spazio convesso tramite la combinazione convessa aA+bB+cC, con a+b+c=1, A,B,C vertici del triangolo della faccia
// la combinazione convessa viene discretizzata per generare (b+1)(b+2)/2 punti, dividendo la combinazione per b: (i/b)*A+(j/b)*B+(b-i-j/b)*C, i+j<=b


namespace PolyhedraLibrary{

PolyhedraMesh TriangolazioneI(PolyhedraMesh& mesh, unsigned int b)
{
    unsigned int NewNumEdges;
    unsigned int NewNumVertices;
    
    // Tetraedro
    if(mesh.NumCell0Ds == 4){ 
        NewNumVertices = 2*b*b + 2; // numero di vertici DI TUTTO IL POLIEDRO 
        NewNumEdges = 6*b*b; // numero di edges DI TUTTO IL POLIEDRO
    }
    
    // Ottaedro
    if(mesh.NumCell0Ds == 6){
        NewNumVertices = 4*b*b + 2;
        NewNumEdges = 12*b*b; 
    }

    // Icosaedro
    if(mesh.NumCell0Ds == 12){
        NewNumVertices = 10*b*b + 2;
        NewNumEdges = 30*b*b;
    }

    double b1 = b;  // Necessario perché, quando si fa i/b avrei divisione intera mentre mi serve che i/b sia un double
    double tol = 1e-6;

    MatrixXd NewCell0DsCoordinates = MatrixXd::Zero(3, NewNumVertices); 
    MatrixXi NewCell1DsExtrema = MatrixXi::Zero(2, NewNumEdges);

    // currentpoint segna a quale punto siamo arrivati, Id serve ad assegnare il corretto id ad ogni punto e Nuovo specifica se il punto esiste già o meno
    unsigned int currentpoint = 0;
    unsigned int currentedge = 0;
    unsigned int Id = 0;
    unsigned int Id_edge1 = 0;
    unsigned int Id_edge2 = 0;
    unsigned int Id_edge3 = 0;
    unsigned int Id_faccia = 0;

    map<pair<unsigned int, unsigned int>, unsigned int> edgeMap;


    bool Nuovo = true;


    vector<vector<unsigned int>> NewFacesVertices;
    vector<vector<unsigned int>> NewFacesEdges;
    vector<unsigned int> NewCell3DsFaces;
    vector<unsigned int> NewCell3DsVertices;
    vector<unsigned int> NewCell3DsEdges;
    
    // Ciclo per ogni terna di vertici (faccia) all'interno della lista delle facce descritte da vertici (Cell2DsVertices)
    for(const auto& faccia : mesh.Cell2DsVertices){

        MatrixXi MatriceTriangolazione(b+1, b+1);
        MatriceTriangolazione.setConstant(-1); 

        // Doppio ciclo per costruire i nuovi punti
        for (unsigned int i = 0; i < b+1; i++){
            for (unsigned int j = 0; j <= b - i; j++){
                
                double alpha = i/b1;
                double beta = j/b1;
                double gamma = (b - i - j)/b1;

                double x_A = mesh.Cell0DsCoordinates(0, faccia[0]);
                double y_A = mesh.Cell0DsCoordinates(1, faccia[0]);
                double z_A = mesh.Cell0DsCoordinates(2, faccia[0]);

                double x_B = mesh.Cell0DsCoordinates(0, faccia[1]);
                double y_B = mesh.Cell0DsCoordinates(1, faccia[1]);
                double z_B = mesh.Cell0DsCoordinates(2, faccia[1]);

                double x_C = mesh.Cell0DsCoordinates(0, faccia[2]);
                double y_C = mesh.Cell0DsCoordinates(1, faccia[2]);
                double z_C = mesh.Cell0DsCoordinates(2, faccia[2]);


                // Costruisco le coordinate x, y, z usando le coordinate baricentriche discretizzate
                double x = alpha * x_A + beta * x_B + gamma * x_C;
                double y = alpha * y_A + beta * y_B + gamma * y_C;
                double z = alpha * z_A + beta * z_B + gamma * z_C;
                
                Nuovo = true;

                // Controllo se il punto esiste già, ovvero se è già nella matrice delle coordinate riempita fino ad ora
                for(unsigned int n = 0; n < currentpoint; n++){
                    if(abs(NewCell0DsCoordinates(0, n)-x) <= tol && abs(NewCell0DsCoordinates(1, n) - y) <= tol && abs(NewCell0DsCoordinates(2, n) - z) <= tol){
                        Nuovo = false;
                        Id = n;
                        break;  // Importante: quando trovo il punto duplicato evito di continuare a ciclare!
                    }
                }

                // Se il punto è nuovo, il suo Id sarà quello del punto a cui siamo arrivati (dato dal contatore currentpoint), altrimenti sarà lo stesso del punto preesistente
                if(Nuovo){
                    Id = currentpoint;
                    NewCell0DsCoordinates(0, Id) = x;
                    NewCell0DsCoordinates(1, Id) = y;
                    NewCell0DsCoordinates(2, Id) = z;
                    NewCell3DsVertices.push_back(Id);
                    currentpoint++;
                }

                // Riempio la matrice di triangolazione che usiamo per unire tra loro i punti
                MatriceTriangolazione(i,j) = Id;

            }
        }

        // Controllando la matrice di triangolazione, creo i due tipi di facce e li aggiungo alla lista di facce NewFaces
        for(unsigned int i=0; i<b; i++){
            for(unsigned int j=0; j<b; j++){
                int v1 = MatriceTriangolazione(i,j);
                int v2 = MatriceTriangolazione(i+1,j);
                int v3 = MatriceTriangolazione(i,j+1);

                if(v1 != -1 && v2 != -1 && v3 != -1){

                    // Uso static_cast<T> per trasformare in quest'occorrenza il tipo int in unsigned int (si può fare solo con tipi "simili")
                    // perché NewFacesVertices è un vector<vector<unsigned int>>
                    NewFacesVertices.push_back({static_cast<unsigned int>(v1), static_cast<unsigned int>(v2), static_cast<unsigned int>(v3)});
                    NewCell3DsFaces.push_back(Id_faccia++);



                    // Il check di esistenza degli edge ora avviene con dizionari: il dizionario è del tipo {v1, v2} : Id_edge, viene fatto 
                    // per ogni combinazione possibile {v1, v2}, {v2, v3}, {v3, v1}

                    // Ordine crescente degli estremi così da non rendere importante l'ordine in cui sono forniti
                    unsigned int a = min(v1, v2);
                    unsigned int b = max(v1, v2);
                    pair<unsigned int, unsigned int> edgeKey = {a, b};

                    // Iteratore sulle chiavi (ovvero le coppie) del dizionario
                    auto it = edgeMap.find(edgeKey);
                    if (it != edgeMap.end()) {
                        // Già esiste -> usa ID esistente
                        Id_edge1 = it->second;
                    } else {
                        // Nuovo edge -> aggiungi alla matrice e crea nuovo elemento nella mappa
                        Id_edge1 = currentedge;
                        NewCell1DsExtrema(0, Id_edge1) = v1;
                        NewCell1DsExtrema(1, Id_edge1) = v2;
                        NewCell3DsEdges.push_back(Id_edge1);
                        edgeMap[edgeKey] = Id_edge1;
                        currentedge++;
                    }


                    a = min(v2, v3);
                    b = max(v2, v3);
                    edgeKey = {a, b};

                    it = edgeMap.find(edgeKey);
                    if (it != edgeMap.end()) {
                        Id_edge2 = it->second;}
                    else {
                        Id_edge2 = currentedge;
                        NewCell1DsExtrema(0, Id_edge2) = v2;
                        NewCell1DsExtrema(1, Id_edge2) = v3;
                        NewCell3DsEdges.push_back(Id_edge2);
                        edgeMap[edgeKey] = Id_edge2;
                        currentedge++;
                    }


                    a = min(v1, v3);
                    b = max(v1, v3);
                    edgeKey = {a, b};

                    it = edgeMap.find(edgeKey);
                    if (it != edgeMap.end()) {
                        Id_edge3 = it->second;}
                    else {
                        Id_edge3 = currentedge;
                        NewCell1DsExtrema(0, Id_edge3) = v1;
                        NewCell1DsExtrema(1, Id_edge3) = v3;
                        NewCell3DsEdges.push_back(Id_edge3);
                        edgeMap[edgeKey] = Id_edge3;
                        currentedge++;
                    }



                    // Aggiungo la faccia caratterizzata dagli edges
                    NewFacesEdges.push_back({Id_edge1, Id_edge2, Id_edge3});
                
                }


                // Costruzione del triangolo "inverso"
                if(j!=0){
                    int v3 = MatriceTriangolazione(i+1, j-1);

                    if(v1 != -1 && v2 != -1 && v3 != -1){

                        NewFacesVertices.push_back({static_cast<unsigned int>(v1), static_cast<unsigned int>(v2), static_cast<unsigned int>(v3)});
                        NewCell3DsFaces.push_back(Id_faccia++);
                        
                        
                        unsigned int a = min(v1, v2);
                        unsigned int b = max(v2, v1);
                        pair<unsigned int, unsigned int> edgeKey = {a, b};

                        auto it = edgeMap.find(edgeKey);
                        if (it != edgeMap.end()) {
                            Id_edge1 = it->second;}
                        else {
                            Id_edge1 = currentedge;
                            NewCell1DsExtrema(0, Id_edge1) = v1;
                            NewCell1DsExtrema(1, Id_edge1) = v2;
                            NewCell3DsEdges.push_back(Id_edge1);
                            edgeMap[edgeKey] = Id_edge1;
                            currentedge++;
                        }
                        
                        
                        a = min(v2, v3);
                        b = max(v2, v3);
                        edgeKey = {a, b};

                        it = edgeMap.find(edgeKey);
                        if (it != edgeMap.end()) {
                            Id_edge2 = it->second;}
                        else {
                            Id_edge2 = currentedge;
                            NewCell1DsExtrema(0, Id_edge2) = v2;
                            NewCell1DsExtrema(1, Id_edge2) = v3;
                            NewCell3DsEdges.push_back(Id_edge2);
                            edgeMap[edgeKey] = Id_edge2;
                            currentedge++;
                        }


                        a = min(v1, v3);
                        b = max(v1, v3);
                        edgeKey = {a, b};

                        it = edgeMap.find(edgeKey);
                        if (it != edgeMap.end()){
                            Id_edge3 = it->second;}
                         else{
                            Id_edge3 = currentedge;
                            NewCell1DsExtrema(0, Id_edge3) = v1;
                            NewCell1DsExtrema(1, Id_edge3) = v3;
                            NewCell3DsEdges.push_back(Id_edge3);
                            edgeMap[edgeKey] = Id_edge3;
                            currentedge++;
                        }
                    
                        NewFacesEdges.push_back({Id_edge1, Id_edge2, Id_edge3});
                        
                    }

                }
                
            }
        }
    }

    vector<unsigned int> NewCell2DsNumEdges(NewFacesEdges.size(), 3);
    
    // Riempio la mesh del poligono triangolato
    mesh.NumCell0Ds = currentpoint;
    mesh.NumCell1Ds = currentedge;
    mesh.NumCell2Ds = NewFacesEdges.size();

    mesh.Cell0DsCoordinates = NewCell0DsCoordinates;
    mesh.Cell0DsId = NewCell3DsVertices;
    mesh.Cell1DsId = NewCell3DsEdges;
    mesh.Cell1DsExtrema = NewCell1DsExtrema;

    mesh.Cell2DsNumEdges = NewCell2DsNumEdges;
    mesh.Cell2DsNumVertices = NewCell2DsNumEdges;
    mesh.Cell2DsId = NewCell3DsFaces;
    mesh.Cell2DsVertices = NewFacesVertices;
    mesh.Cell2DsEdges = NewFacesEdges;

    mesh.Cell3DsVertices = NewCell3DsVertices;
    mesh.Cell3DsEdges = NewCell3DsEdges;
    mesh.Cell3DsFaces = NewCell3DsFaces;

    mesh.Cell3DsNumVertices = currentpoint;
    mesh.Cell3DsNumEdges = currentedge;
    mesh.Cell3DsNumFaces = NewFacesEdges.size();

    return mesh; 

}

}
