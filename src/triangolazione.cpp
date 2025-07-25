// TRIANGOLAZIONE.CPP

#include <iostream>
#include "Eigen/Eigen"
#include "PolyhedraMesh.hpp"
#include "Utils.hpp"
#include <map> // dizionari
#include <set>
#include <cmath>
#include <algorithm>  // per min e max e sort()


using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;

namespace PolyhedraLibrary{

    namespace FunzioniAusiliarie{

        // Funzione che verifica se un punto P giace sul segmento AB
        bool PuntoSuSpigolo(const Vector3d& A, const Vector3d& B, const Vector3d& P, double tol = 1e-12) {
            Vector3d AB = B - A;
            Vector3d AP = P - A;

            double ab_len = AB.norm();
            double cross_norm = AB.cross(AP).norm();
            double dot = AB.dot(AP);

            bool allineato = (cross_norm / ab_len) < tol;
            bool compreso = (dot >= -tol) && (dot <= AB.squaredNorm() + tol);

            return allineato && compreso;

        }

        // Funzione che aggiunge un vertice se non è già presente, aggiungendolo alla matrice di coordinate e aggiungendone l'Id alla lista di Id
        bool AggiungiVertice(const Vector3d& P, 
                    map<array<double, 3>, unsigned int>& VerticesMap, 
                    MatrixXd& Coordinates, 
                    vector<unsigned int>& Ids,
                    unsigned int& currentpoint){

                    auto risultato = VerticesMap.try_emplace({P(0), P(1), P(2)}, currentpoint);
                    if(risultato.second){
                        unsigned int id = risultato.first->second;
                        Coordinates.col(id) = P;
                        Ids.push_back(id);
                        currentpoint++;
                        return true;
                    }
                    else return false;
        }

        // Funzione che prova ad aggiungere edges se non sono già presenti
        bool AggiungiEdge(const unsigned int id_P, 
                const unsigned int id_Q, 
                map<pair<unsigned int, unsigned int>, unsigned int>& EdgesMap,
                MatrixXi& Extrema,
                vector<unsigned int>& Ids, 
                unsigned int& currentedge){
                
                auto risultato = EdgesMap.try_emplace({min(id_P, id_Q), max(id_P, id_Q)}, currentedge);
                if(risultato.second){
                    Extrema(0, currentedge) = min(id_P, id_Q);
                    Extrema(1, currentedge) = max(id_P, id_Q);
                    Ids.push_back(currentedge);
                    currentedge++;
                    return true;
                }
                else return false;
        }

        // Funzione che aggiunge le facce sia etichettate per vertici che per edges: prova ad aggiungere faccia per vertici, se successo allora recupera
        // gli id dei lati corrispondenti dalla mappa EdgesMap e aggiunge anche la faccia per lati, poi incrementa il contatore una volta sola
        bool AggiungiFaccia(const unsigned int id_A, 
            const unsigned int id_B, 
            const unsigned int id_C,
            map<array<unsigned int, 3>, unsigned int>& FacesVerticesMap,
            const map<pair<unsigned int, unsigned int>, unsigned int>& EdgesMap,
            vector<unsigned int>& Ids,
            vector<vector<unsigned int>>& FacesVertices,
            vector<vector<unsigned int>>& FacesEdges,
            unsigned int& currentface){
                array<unsigned int, 3> faccia = {id_A, id_B, id_C};
                sort(faccia.begin(), faccia.end()); // sort ordina in modo crescente, prende in input due puntatori, uno al primo e uno dopo l'ultimo elemento
                auto risultato = FacesVerticesMap.try_emplace(faccia, currentface);
                if(risultato.second){
                    vector<unsigned int> face(faccia.begin(), faccia.end());    // Converto in vector perché FacesVertices contiene vectors
                    FacesVertices.push_back(face);
                    unsigned int id_AB = EdgesMap.at({min(id_A, id_B), max(id_A, id_B)});   // Uso .at() per l'accesso perché ho un const e non posso usare []
                    unsigned int id_BC = EdgesMap.at({min(id_C, id_B), max(id_C, id_B)});
                    unsigned int id_CA = EdgesMap.at({min(id_A, id_C), max(id_A, id_C)});
                    FacesEdges.push_back({id_AB, id_BC, id_CA});
                    Ids.push_back(currentface);
                    currentface++;
                    return true;
                }
                else return false;
        }
    }

// Funzione di triangolazione del poliedro tipo I
// generiamo uno spazio convesso tramite la combinazione convessa aA+bB+cC, con a+b+c=1, A,B,C vertici del triangolo della faccia
// la combinazione convessa viene discretizzata per generare (b+1)(b+2)/2 punti, dividendo la combinazione per b: (i/b)*A+(j/b)*B+(b-i-j/b)*C, i+j<=b

PolyhedraMesh TriangolazioneI(PolyhedraMesh& mesh, unsigned int b)
{
    using namespace FunzioniAusiliarie;

    unsigned int NewNumEdges;
    unsigned int NewNumVertices;
    
    switch(mesh.NumCell0Ds){
        case 4: // Tetraedro
            NewNumVertices = 2*b*b + 2; 
            NewNumEdges = 6*b*b;
            break;

        case 6: // Ottaedro
            NewNumVertices = 4*b*b + 2;
            NewNumEdges = 12*b*b;
            break;
        
        case 12: // Icosaedro
            NewNumVertices = 10*b*b + 2;
            NewNumEdges = 30*b*b;
            break;
    }

    double b1 = b;  // Necessario perché, quando si fa i/b avrei divisione intera mentre mi serve che i/b sia un double

    unsigned int currentpoint = 0;
    MatrixXd NewCell0DsCoordinates = MatrixXd::Zero(3, NewNumVertices);
    map<array<double, 3>, unsigned int> VerticesMap; 

    MatrixXi NewCell1DsExtrema = MatrixXi::Zero(2, NewNumEdges);

    // currentpoint segna a quale punto siamo arrivati, Id serve ad assegnare il corretto id ad ogni punto
    unsigned int currentedge = 0;
    unsigned int Id = 0;
    unsigned int Id_edge1 = 0;
    unsigned int Id_edge2 = 0;
    unsigned int Id_edge3 = 0;
    unsigned int Id_faccia = 0;

    map<pair<unsigned int, unsigned int>, unsigned int> edgeMap;

    vector<vector<unsigned int>> NewFacesVertices;
    vector<vector<unsigned int>> NewFacesEdges;
    vector<unsigned int> NewCell3DsFaces;
    vector<unsigned int> NewCell3DsVertices;
    for (unsigned int i=0; i < mesh.NumCell0Ds; i++){
        NewCell3DsVertices.push_back(i);
    }
    vector<unsigned int> NewCell3DsEdges;

    // Copia i vertici originali nei primi ID
    for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i) {
        Vector3d P = mesh.Cell0DsCoordinates.col(i);
        NewCell0DsCoordinates.col(i) = mesh.Cell0DsCoordinates.col(i);
        AggiungiVertice(P, VerticesMap, NewCell0DsCoordinates, NewCell3DsVertices, currentpoint);
    }
    
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

                Vector3d A = mesh.Cell0DsCoordinates.col(faccia[0]);
                VerticesMap.try_emplace({A(0), A(1), A(2)}, faccia[0]);
                Vector3d B = mesh.Cell0DsCoordinates.col(faccia[1]);
                VerticesMap.try_emplace({B(0), B(1), B(2)}, faccia[1]);
                Vector3d C = mesh.Cell0DsCoordinates.col(faccia[2]);
                VerticesMap.try_emplace({C(0), C(1), C(2)}, faccia[2]);

                // Costruisco il nuovo punto P usando le coordinate baricentriche discretizzate
                Vector3d P = alpha * A + beta * B + gamma * C;
                
                bool aggiunto = AggiungiVertice(P, VerticesMap, NewCell0DsCoordinates, NewCell3DsVertices, currentpoint);
                Id = aggiunto ? currentpoint-1 : VerticesMap[{P(0), P(1), P(2)}];

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

                    bool aggiunto = AggiungiEdge(v1, v2, edgeMap, NewCell1DsExtrema, NewCell3DsEdges, currentedge);
                    Id_edge1 = aggiunto ? currentedge - 1 : edgeMap[{min(v1, v2), max(v1, v2)}];    
                    // If inline: Id_edge1 = (condizione)? ( = se condizione verificata) : (= se condizione non verificata)

                    aggiunto = AggiungiEdge(v2, v3, edgeMap, NewCell1DsExtrema, NewCell3DsEdges, currentedge);
                    Id_edge2 = aggiunto ? currentedge - 1 : edgeMap[{min(v2, v3), max(v2, v3)}];

                    aggiunto = AggiungiEdge(v3, v1, edgeMap, NewCell1DsExtrema, NewCell3DsEdges, currentedge);
                    Id_edge3 = aggiunto ? currentedge - 1 : edgeMap[{min(v3, v1), max(v3, v1)}];

                    // Aggiungo la faccia caratterizzata dagli edges
                    NewFacesEdges.push_back({Id_edge1, Id_edge2, Id_edge3});
                
                }

                // Costruzione del triangolo "inverso"
                if(j!=0){
                    int v3 = MatriceTriangolazione(i+1, j-1);

                    if(v1 != -1 && v2 != -1 && v3 != -1){

                        NewFacesVertices.push_back({static_cast<unsigned int>(v1), static_cast<unsigned int>(v2), static_cast<unsigned int>(v3)});
                        NewCell3DsFaces.push_back(Id_faccia++);

                        bool aggiunto = AggiungiEdge(v1, v2, edgeMap, NewCell1DsExtrema, NewCell3DsEdges, currentedge);
                        Id_edge1 = aggiunto ? currentedge - 1 : edgeMap[{min(v1, v2), max(v1, v2)}];

                        aggiunto = AggiungiEdge(v2, v3, edgeMap, NewCell1DsExtrema, NewCell3DsEdges, currentedge);
                        Id_edge2 = aggiunto ? currentedge - 1 : edgeMap[{min(v2, v3), max(v2, v3)}];

                        aggiunto = AggiungiEdge(v3, v1, edgeMap, NewCell1DsExtrema, NewCell3DsEdges, currentedge);
                        Id_edge3 = aggiunto ? currentedge - 1 : edgeMap[{min(v3, v1), max(v3, v1)}];
 
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

    mesh.Cell0DsCoordinates = NewCell0DsCoordinates.leftCols(currentpoint);
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


PolyhedraMesh TriangolazioneII(PolyhedraMesh& mesh, unsigned int b){

    using namespace FunzioniAusiliarie;

    unsigned int NewNumE;
    unsigned int NewNumV;

    // Salvo gli spigoli del poliedro di partenza, perché serviranno per verificare i punti da inserire, in un vettore di coppie, ogni
    // elemento della coppia è un vettore Eigen di coordinate
    vector<pair<Vector3d, Vector3d>> EdgesOriginali;
    for(unsigned int n = 0; n < mesh.NumCell1Ds; n++){
        Vector3d P = mesh.Cell0DsCoordinates.col(mesh.Cell1DsExtrema(0,n));
        Vector3d Q = mesh.Cell0DsCoordinates.col(mesh.Cell1DsExtrema(1,n));
        EdgesOriginali.push_back({P,Q});
    }

    switch(mesh.NumCell0Ds){
        case 4: // Tetraedro
            NewNumV = 6*b*b + 6*b + 2;
            NewNumE = 18*b*b + 18*b + 2 + 1000;
            break;

        case 6: // Ottaedro
            NewNumV = 12*b*b + 12*b + 2;
            NewNumE = 36*b*b + 36*b + 2 + 1000;
            break;
        
        case 12: // Icosaedro
            NewNumV = 30*b*b + 30*b + 2;
            NewNumE = 90*b*b + 90*b + 2 + 1000;
            break;
    }

    // La triangolazione II si basa sulla I: triangolo di tipo I ogni faccia e poi, per ogni nuova faccia, effettuo la 
    // triangolazione II nel caso b = 1
    
    TriangolazioneI(mesh, b);
    unsigned int currentpoint = mesh.NumCell0Ds;
    unsigned int currentedge = 0;
    unsigned int currentface = 0;

    // Usiamo dizionari per verificare se gli elementi esistono già, le chiavi sono le coordinate o gli id di origine e fine o gli id 
    // caratteristici di una faccia, così da poter scorrerle rapidamente con try_emplace per verificare duplicati, e i valori sono gli id
    map<array<double, 3>, unsigned int> VerticesMap;
    map<pair<unsigned int, unsigned int>, unsigned int> EdgesMap;
    map<array<double, 3>, unsigned int> MidToBar;
    map<array<unsigned int, 3>, unsigned int> FacesVerticesMap; // Usiamo mappa solo basata su Id di vertici perché da questi risaliamo a edges
    bool UnireBar;

    // conservativeResize aumenta la dimensione mantenendo gli elementi già presenti (i punti della triangolazione I)
    mesh.Cell0DsCoordinates.conservativeResize(3, NewNumV);
    MatrixXi NewCell1DsExtrema = MatrixXi::Zero(2, NewNumE);
    vector<unsigned int> NewCell1DsId;
    vector<unsigned int> NewCell2DsId;
    vector<vector<unsigned int>> NewFacesVertices;
    vector<vector<unsigned int>> NewFacesEdges;


    for(const auto& faccia : mesh.Cell2DsVertices){
        
        Vector3d A = mesh.Cell0DsCoordinates.col(faccia[0]);
        VerticesMap.try_emplace({A(0), A(1), A(2)}, faccia[0]);  // try_emplace verifica se è già esistente la chiave, altrimenti inserisce coppia chiave-valore
        unsigned int id_A = VerticesMap[{A(0), A(1), A(2)}];
        
        Vector3d B = mesh.Cell0DsCoordinates.col(faccia[1]);
        VerticesMap.try_emplace({B(0), B(1), B(2)}, faccia[1]);
        unsigned int id_B = VerticesMap[{B(0), B(1), B(2)}];

        Vector3d C = mesh.Cell0DsCoordinates.col(faccia[2]);
        VerticesMap.try_emplace({C(0), C(1), C(2)}, faccia[2]);
        unsigned int id_C = VerticesMap[{C(0), C(1), C(2)}];


        // Calcolo punti medi dei lati e baricentro e provo a inserirli nella mappa, se riesco a inserirli allora sono nuovi punti e
        // incremento il contatore currentpoint, gli ID dei punti della nuova triangolazione saranno tutti maggiori degli ID dei precedenti
        // poiché currentpoint parte dall'ultimo Id della triangolazione I

        // Non controllo se il baricentro esiste già perché è sempre interno al triangolo quindi mai in comune
        Vector3d baricentro = (A+B+C)/3;
        unsigned int id_bar;
        bool aggiunto = AggiungiVertice(baricentro, VerticesMap, mesh.Cell0DsCoordinates, mesh.Cell0DsId, currentpoint);
        if(aggiunto) id_bar = currentpoint - 1;
        else id_bar = VerticesMap[{baricentro(0), baricentro(1), baricentro(2)}];
        AggiungiEdge(id_bar, id_A, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
        AggiungiEdge(id_bar, id_B, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
        AggiungiEdge(id_bar, id_C, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
        
        // Per verificare se un punto va aggiunto o meno, verifico se giace su un edge presente all'inizio, in tal caso lo aggiungo, 
        // altrimenti si tratta di un punto interno non presente nella triangolazione II

        Vector3d m_AB = (A+B)/2;
        Vector3d m_BC = (B+C)/2;
        Vector3d m_CA = (C+A)/2;
        array<Vector3d, 3> midpoints = {m_AB, m_BC, m_CA};

        for(auto& m : midpoints){
            
            UnireBar = true;
            for(const auto& edge : EdgesOriginali){
                if(PuntoSuSpigolo(edge.first, edge.second, m)){
                    unsigned int id = currentpoint;
                    aggiunto = AggiungiVertice(m, VerticesMap, mesh.Cell0DsCoordinates, mesh.Cell0DsId, currentpoint);
                    if(!aggiunto) id = VerticesMap[{m(0), m(1), m(2)}];

                    if(PuntoSuSpigolo(A, B, m)){
                        AggiungiEdge(id, id_A, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
                        AggiungiEdge(id, id_B, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
                        AggiungiEdge(id, id_bar, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);

                        AggiungiFaccia(id, id_A, id_bar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                        AggiungiFaccia(id, id_B, id_bar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                    }

                    if(PuntoSuSpigolo(B, C, m)){
                        AggiungiEdge(id, id_B, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
                        AggiungiEdge(id, id_C, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
                        AggiungiEdge(id, id_bar, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);

                        AggiungiFaccia(id, id_C, id_bar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                        AggiungiFaccia(id, id_B, id_bar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                    }

                    if(PuntoSuSpigolo(C, A, m)){
                        AggiungiEdge(id, id_C, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
                        AggiungiEdge(id, id_A, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);
                        AggiungiEdge(id, id_bar, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge);

                        AggiungiFaccia(id, id_A, id_bar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                        AggiungiFaccia(id, id_C, id_bar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                    }

                    UnireBar = false;
                    break;  // Se il punto va aggiunto, smetto di verificare se vada aggiunto
                }
            }

            // Questi if si occupano di gestire l'unione tra due baricentri di facce vicine, questo accade se il lato in comune è un lato interno
            // e non uno degli spigoli originali. In tal caso, non viene aggiunto il punto medio di quel lato, e viene salvata in un dizionario la 
            // chiave {coordinate punto medio} con valore id_baricentro della faccia. Se provando ad inserire le coordinate del punto medio vedo che
            // in realtà esse sono già nel dizionario, vuol dire che ho già provato ad aggiungere il punto medio da un'altra faccia (comunicante), quindi che nel dizionario
            // la coppia è già presente ma ha come valore il baricentro di quest'altra faccia. Devo quindi recuperare gli id dei due baricentri, id_bar e id_altrobar
            // e creare l'edge tra i due.
        
            if(UnireBar){
                auto risultato = MidToBar.try_emplace({m(0), m(1), m(2)}, id_bar);
                if(!risultato.second){
                    unsigned int id_altrobar = MidToBar[{m(0), m(1), m(2)}];
                    if(AggiungiEdge(id_bar, id_altrobar, EdgesMap, NewCell1DsExtrema, NewCell1DsId, currentedge)){

            // Se ho creato l'edge con successo, devo creare le facce che hanno come vertici i due baricentri e un terzo punto, che sarà un punto della faccia di partenza.
            // Le facce saranno due e condivideranno un lato, per cui ci sarà un vertice della faccia di partenza da escludere: per capire quale, 
            // guardo quale vertice è allineato ad entrambi i baricentri e lo escludo, dopodiché aggiungo le facce composte dai due baricentri e dai punti rimasti. 

                        if(PuntoSuSpigolo(A, m, baricentro) || PuntoSuSpigolo(A, baricentro, m)){
                            AggiungiFaccia(id_bar, id_B, id_altrobar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                            AggiungiFaccia(id_bar, id_C, id_altrobar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                        }
                        else if(PuntoSuSpigolo(B, m, baricentro) || PuntoSuSpigolo(B, baricentro, m)){
                            AggiungiFaccia(id_bar, id_A, id_altrobar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                            AggiungiFaccia(id_bar, id_C, id_altrobar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                        }
                        else if(PuntoSuSpigolo(C, m, baricentro) || PuntoSuSpigolo(C, baricentro, m)){
                            AggiungiFaccia(id_bar, id_A, id_altrobar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                            AggiungiFaccia(id_bar, id_B, id_altrobar, FacesVerticesMap, EdgesMap, NewCell2DsId, NewFacesVertices, NewFacesEdges, currentface);
                        }
                    }
                }
            }
        }
    }

    mesh.NumCell0Ds = currentpoint;
    mesh.NumCell1Ds = currentedge;
    mesh.NumCell2Ds = currentface;

    // Per difficoltà nelle formule di calcolo degli edges del poliedro triangolato, uso una dimensione più ampia di default e poi "taglio" tenendo solo le colonne necessarie
    mesh.Cell1DsExtrema = NewCell1DsExtrema.leftCols(currentedge);
    mesh.Cell1DsId = NewCell1DsId;

    mesh.Cell2DsId = NewCell2DsId;
    mesh.Cell2DsEdges = NewFacesEdges;
    mesh.Cell2DsVertices = NewFacesVertices;

    vector<unsigned int> NewCell2DsNum(currentface, 3); 
    mesh.Cell2DsNumEdges = NewCell2DsNum;
    mesh.Cell2DsNumVertices = NewCell2DsNum;

    mesh.Cell3DsVertices = mesh.Cell0DsId;
    mesh.Cell3DsEdges = mesh.Cell1DsId;
    mesh.Cell3DsFaces = mesh.Cell2DsId;

    mesh.Cell3DsNumVertices = currentpoint;
    mesh.Cell3DsNumEdges = currentedge;
    mesh.Cell3DsNumFaces = currentface;


    return mesh;
}

}