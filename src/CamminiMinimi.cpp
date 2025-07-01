//Cammini Minimi.cpp

#include <iostream>
#include <vector>
//#include <unordered_map> 
#include <queue> //per priority queue
#include <iomanip> // per setprecision
#include <cmath> // per sqrt()
#include <limits> //per infinity
#include <algorithm>
#include <Eigen/Dense>
#include "PolyhedraMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;


using ListaAdiacenza = vector<vector<pair<unsigned int, double>>>;

//costruisco una lista di adiacenza (un grafo) dal PolyhedraMesh
//dove per ogni nodo ho una lista di nodi adiacenti con i relativi pesi
ListaAdiacenza CreaListaAdiacenza(const PolyhedraMesh& mesh) 
{
    ListaAdiacenza grafo(mesh.NumCell0Ds); //inizializzo con tanti elementi quanti sono i vertici

    for (unsigned int i = 0; i < mesh.NumCell1Ds; ++i) //per ogni edge, prendo i vertici che ne identificano l'inizio e la fine
    {
        unsigned int vertice1 = mesh.Cell1DsExtrema(0, i); //id del vertice di origine dell'edge in Cell1Ds
        unsigned int vertice2 = mesh.Cell1DsExtrema(1, i); //id del vertice di fine dell'edge in Cell1Ds

        //distanza tra v1 e v2 (peso dell'arco)
        cout << setprecision(16);
        const double dx = mesh.Cell0DsCoordinates(0, vertice1) - mesh.Cell0DsCoordinates(0, vertice2);
        const double dy = mesh.Cell0DsCoordinates(1, vertice1) - mesh.Cell0DsCoordinates(1, vertice2);
        const double dz = mesh.Cell0DsCoordinates(2, vertice1) - mesh.Cell0DsCoordinates(2, vertice2);
        double peso = sqrt(dx*dx + dy*dy + dz*dz);

        // Aggiungi entrambe le direzioni (il grafo non è diretto)
        grafo[vertice1].emplace_back(vertice2, peso);
        grafo[vertice2].emplace_back(vertice1, peso);
    }
    
    return grafo;
}

//con pair, Dijkstra restituisce la coppia (percorso (cioé vertici visitati), distanza totale)
//percorso = vettore di unsigned int 
//distanza totale = double
pair<vector<unsigned int>, double> Dijkstra(const ListaAdiacenza& grafo, unsigned int nodo_iniziale, unsigned int nodo_finale) 
{
    //inizializzo le strutture dati necessarie per l'algoritmo di Dijkstra
    // ! grafo.size() è numCell0Ds, cioé il numero di vertici
    vector<double> distanze(grafo.size(), numeric_limits<double>::infinity()); //vettore delle distanze minime, inizializzato ad infinito
    vector<unsigned int> predecessori(grafo.size(), -1); //vettore per tenere traccia dei predecessori di ogni nodo nel percorso (ovvero dei nodi visitati)
    priority_queue<pair<double, unsigned int>, 
                vector<pair<double, unsigned int>>,
                greater<pair<double, unsigned int>>> coda_prioritaria; //coda di priorità per i nodi da visitare, ordinata per distanza crescente
    //la coda ha coppie di tipo (distanza, nodo), dove la distanza è un double e il nodo è un unsigned int

    
    cout << setprecision(16); //imposto la precisione per i numeri in virgola mobile (forse inutile in tutto il codice, perché sono double?)
    predecessori[nodo_iniziale] = nodo_iniziale; //il predecessore del nodo iniziale è se stesso
    distanze[nodo_iniziale] = 0.0; //distanza del nodo iniziale da se stesso è 0  
    coda_prioritaria.emplace(0.0, nodo_iniziale); //inserisco il nodo iniziale nella coda di priorità con distanza 0
    //enqueue in priority queue è emplace

    /* for (i = 0, i < grafo.size(); ++i) //per ogni nodo del grafo
    {
        coda_prioritaria.enqueue({distanze[i], i}); //inserisco nella coda di priorità il nodo con la sua distanza
    } */ //ridondante e spreca memoria

    while (!coda_prioritaria.empty())
    {
        //ogni elemento della priority queue è una coppia (distanza, nodo), 
        //quindi coda_prioritaria.top() restituisce il prossimo elemento della coda (una coppia dist_attuale, u)

        //dequeue in priority queue è top/pop
        auto [dist_attuale, u] = coda_prioritaria.top(); //leggo l'elemento con la distanza minima dalla coda di priorità
        coda_prioritaria.pop(); //rimuovo l'elemento dalla coda di priorità
         
        
        if (u == nodo_finale) break; //sono già arrivata a destinazione

        if (dist_attuale > distanze[u]) continue; //se la distanza attuale è maggiore della distanza minima trovata finora, salto questo nodo
 
        for (const auto& [v, pesi] : grafo[u]) 
        {
            if (distanze[v] > distanze[u] + pesi) 
            {
                distanze[v] = distanze[u] + pesi;
                predecessori[v] = u;
                coda_prioritaria.emplace(distanze[v], v);
            }
        }
    } 

    //ricompongo il percorso
    vector<unsigned int> percorso;
    if (distanze[nodo_finale] == numeric_limits<double>::infinity()) 
    {
        return {percorso, -1.0}; //se non ho trovato nessun percorso
    }
    
    //at = inizio dal nodo finale e risalgo i predecessori fino al nodo iniziale
    for (unsigned int at = nodo_finale; at != nodo_iniziale; at = predecessori[at])
    {
        percorso.push_back(at);
    } 

    percorso.push_back(nodo_iniziale); //aggiungo il nodo iniziale al percorso

    reverse(percorso.begin(), percorso.end()); //inverto il percorso per avere l'ordine corretto (dall'inizio alla fine, mentre l'avevo costruito dalla fine all'inizio)


    return {percorso, distanze[nodo_finale]};
    }

