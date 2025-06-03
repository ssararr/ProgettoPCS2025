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

    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds); //riservo spazio per le coordinate x,y,z in una matrice 4xvertici

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
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds); //matrice di dimensione 2 x numero di spigoli per salvare gli estremi 

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
bool ImportCell2Ds(PolyhedraMesh& mesh, string Poliedro)
{
    ifstream file;
    file.open("./" << Poliedro << "_" << "Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))

        listLines.push_back(line);  // Salva nella lista listlines quello che ha letto in line con getline

    file.close();

    listLines.pop_front();  // Toglie header

    mesh.NumCell2Ds = listLines.size();     // Riempio NumCell2Ds con il numero di celle quindi di facce

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "Non esistono celle 2D" << endl;
        return false;
    }

    // Riservo spazio in memoria per riempire i vettori grandi tante quante celle
    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);    
    mesh.Cell2DsNumVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsNumEdges.reserve(mesh.NumCell2Ds);

    // Vectors dinamici contenenti i vertices e gli edges 
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    // Cicla su ogni linea in listlines
    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int NumVertices;
        unsigned int NumEdges;
        double tmp;
        vector<unsigned int> vertices;
        vector<unsigned int> edges;

        converter >> id >> NumVertices

        // Ciclo per aggiungere ogni vertice alla lista dei vertici
        vertices.reserve(NumVertices);  // Il numero di vertici dipende dalla variabile letta dal file
        for(unsigned int i = 0; i < NumVertices; i++){  // Il ciclo itera in base al numero di vertici riportati nel .csv
            converter >> tmp;
            vertices.push_back(tmp);
        }

        converter >> NumEdges;
        // Ciclo per aggiungere ogni arco alla lista degli archi
        edges.reserve(NumEdges);
        for(unsigned int i = 0; i < NumEdges; i++)
            {
<<<<<<< Updated upstream
            converter >> tmp;
=======
            converter >> separator >> tmp;
>>>>>>> Stashed changes
            edges.push_back(tmp);}

        // Aggiungo alle liste di id e NumVertices, NumEdges i dati letti
        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DNumVertices.push_back(NumVertices);
        mesh.Cell2DNumEdges.push_back(NumEdges);

        // Aggiungo alle liste contenenti ID di vertici e archi le liste di vertici e archi
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);

    }

    return true;
}

// ***************************************************************************
bool ImportCell3Ds(PolyhedraMesh& mesh, string Poliedro) //implemento la funzione che deve raccogliere i dati relativi agli SPIGOLI
{
    ifstream file;
    file.open("./" << Polidero << "_Cell3Ds.csv"); //apro il file da cui leggiamo i dati

    if(file.fail()) //gestiamo l'apertura del file
        return false;

    list<string> listLines; //creo una lista di stringhe

    string line;
    while (getline(file, line)) //leggo il file linea per linea
        listLines.push_back(line); //aggiungo ogni linea del file come elemento della lista

    file.close(); //chiudo il file

    // rimuovo l'header della lista
    listLines.pop_front();

    mesh.NumCell3Ds = listLines.size(); //imposto il numero di spigoli uguale alla dimensione della lista 

    if (mesh.NumCell3Ds == 0) //nel caso la dimensione sia nulla gestisco l'errore
    {
        cerr << "There is no cell 3D" << endl;
        return false;
    }

    // Alloco memoria per gli ID e le coordinate dei punti
    mesh.Cell3DsId.reserve(mesh.NumCell3Ds); //vettore di dimensione numero di spigoli per salvare gli ID degli spigoli

    for (const string& line : listLines) //per ogni riga della lista, presa per riferimento 
    {
        istringstream converter(line); //oggetto istringstream che mi permette di manipolare le righe di stringhe

        unsigned int id;               
        unsigned int NumVertices;
        unsigned int NumEdges;
        unsigned int NumFaces;

        converter >>  id;     

        converter >>  NumVertices;
        mesh.Cell3dVertices.reserve(NumVertices); //riservo lo spazio pari al numero totale di poliedri
        for (unsigned int i = 0; i < NumVertices; i++)
        {
            unsigned int vertex;
            converter >> vertex;
            mesh.Cell3dVertices.push_back(vertex); //aggiungo i vertici al vettore di vertici diel poliedro
        }

        converter >>  NumEdges;
        mesh.Cell3dEdges.reserve(NumEdges); //riservo lo spazio pari al numero totale di poliedri
        for (unsigned int j = 0; j < NumEdges; j++)
        {
            unsigned int edge;
            converter >> edge;
            mesh.Cell3dEdges.push_back(edge); //aggiungo gli spigoli al vettore di spigoli del poliedro
        }

        converter >>  NumFaces;
        mesh.Cell3dFaces.reserve(NumFaces); //riservo lo spazio pari al numero totale di poliedri
        for (unsigned int k = 0; k < NumFaces; k++)
        {
            unsigned int face;
            converter >> face;
            mesh.Cell3dFaces.push_back(face); //aggiungo le facce al vettore di facce del poliedro
        }

        mesh.Cell3DsId.push_back(id);
        mesh.Cell3DsNumVertices = NumVertices;
        mesh.Cell3DsNumEdges = NumEdges;
        mesh.Cell3DsNumFaces = NumFaces;
    }
    
    return true; //restituisco true se la funzione è andata a buon fine
}

}
