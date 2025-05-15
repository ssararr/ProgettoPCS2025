//POLYHEDRA UTILS CPP
#include "Utils.hpp"
#include <iostream>
#include <fstream> //libreria per la gestione dei file
#include <sstream> //libreria per manipolare le stringhe
#include <string> //libreria per le stringhe
//(la libreria Eigen è già inclusa in Utils.hpp e PolyhedraMesh.hpp)

namespace PolyhedraLibrary 
{
bool ImportMesh(PolyhedraMesh& mesh, string Poliedro) //funzione che importa una mesh poligonale. se ci sono errori e una funzione fallisce viene restituito false
{

    if(!ImportCell0Ds(mesh), string Poliedro) //importo la funzione per i dati relativi ai PUNTI 
        return false;

    if(!ImportCell1Ds(mesh), string Poliedro) //importo la funzione per i dati relativi agli SPIGOLI
        return false;

    if(!ImportCell2Ds(mesh), string Poliedro) //importo la funzione per i dati relativi ai POLIGONI
        return false;
    
    if(!ImportCell3Ds(mesh), string Poliedro) //importo la funzione per i dati relativi ai POLIEDRI
        return false;
}



// ***************************************************************************
bool ImportCell0Ds(PolyhedraMesh& mesh, string Poliedro) //funzione per importare i punti 0D
//la funzione prende in input la mesh e il nome del poliedro
{
    ifstream file;
    file.open("./" << Poliedro << "_Cell0Ds.csv"); //apro il file

    if(file.fail()) //controllo se il file è aperto
        return false;

    //leggo il file riga per riga e inserisco le righe in una lista
    list<string> listLines; //creo una lista di stringhe
    string line;
    while (getline(file, line)) //leggo il file riga per riga
        listLines.push_back(line); //inserisco le righe nella lista
    file.close(); //chiudo il file

    //rimuovo la prima riga
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size(); //numero di righe, cioè numero di vertici

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl; //se non ci sono righe
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds); //riservo spazio per gli Id, cioè riservo spazio per il numero di vertici (punti 0D)

    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(4, mesh.NumCell0Ds); //riservo spazio per le coordinate x,y,z in una matrice 4xvertici

    for (const string& line : listLines) //per ogni riga della lista, passata per riferimento
    {
        istringstream converter(line); //creo un oggetto di tipo istringstream per convertire la stringa in un flusso di dati
        unsigned int id; //Id del punto
        Vector3d coord; //coordinate del punto

        converter >> id >> mesh.Cell0DsCoordinates(0, id) >> mesh.Cell0DsCoordinates(1, id) >> mesh.Cell0DsCoordinates(2, id); //leggo l'id e le coordinate x,y,z 
        mesh.Cell0DsId.push_back(id); //inserisco l'id in Cell0DsId 
    }
    
    return true;
}

// ***************************************************************************

bool ImportCell1Ds(PolyhedraMesh& mesh, string Poliedro) //implemento la funzione che deve raccogliere i dati relativi agli SPIGOLI
{
    ifstream file;
    file.open("./" << Polidero << "_Cell1Ds.csv"); //apro il file da cui leggiamo i dati

    if(file.fail()) //gestiamo l'apertura del file
        return false;

    list<string> listLines; //creo una lista di stringhe

    string line;
    while (getline(file, line)) //leggo il file linea per linea
        listLines.push_back(line); //aggiungo ogni linea del file come elemento della lista

    file.close(); //chiudo il file

    // rimuovo l'header della lista
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size(); //imposto il numero di spigoli uguale alla dimensione della lista 

    if (mesh.NumCell1Ds == 0) //nel caso la dimensione sia nulla gestisco l'errore
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    // Alloco memoria per gli ID e le coordinate dei punti
    mesh.Cell1DsId.reserve(mesh.NumCell1Ds); //vettore di dimensione numero di spigoli per salvare gli ID degli spigoli 
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds); //matrice di dimensione 2x numero di spigoli per salvare gli estremi 

    for (const string& line : listLines) //per ogni riga della lista, presa per riferimento 
    {
        istringstream converter(line); //oggetto istringstream che mi permette di manipolare le righe di stringhe

        unsigned int id;               
        Vector2i vertices; //definisco un vettore di dimensione 2x1 per gli estremi
        converter >>  id;                
        converter >>  mesh.Cell1DsExtrema(0, id) /*coordinatax*/;        
        converter >>  mesh.Cell1DsExtrema(1, id) /*coordinatay*/;
        mesh.Cell1DsId.push_back(id); //aggiungo id all'ultima posizione del vettore Cell1DsId
    }

    return true; //restituisco true se la funzione è andata a buon fine
}

// ***************************************************************************
