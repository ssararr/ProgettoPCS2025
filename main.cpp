// MAIN POLYHEDRA

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include "PolyhedraMesh.hpp"
#include "Utils.hpp"
#include "triangolazione.hpp"
#include "Duale.hpp"
#include "output.hpp"
#include "src/Proiezione.hpp"
#include "CamminiMinimi.hpp"
#include "ExportParaview/UCDUtilities.hpp"
#include <fstream>

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;
using namespace Gedim;

// La funzione prende la coppia {p, q} e una variabile booleana che serve
//per capire se dobbiamo fare il duale (duale = true) o no (duale = false)

string ClassificaPoliedro(unsigned int& p, unsigned int& q, bool& duale){
    if(p == 3){
        duale = false;
        if(q == 3){
            return "Tetraedro";
        }

        else if(q == 4){
            return "Ottaedro";
        }

        else if(q == 5){
            return "Icosaedro";
        }

        else{
            return "Errore";
        }
    }

    // Cubo, ma scambiamo subito p, q e ci segnamo di dover passare al duale
    else if(p == 4 && q == 3){
        duale = true;
        unsigned int tmp = p;
        p = q;
        q = tmp;
        return "Ottaedro";
    }

    // Dodecaedro, ma scambiamo subito p, q e ci segnamo di dover passare al duale
    else if(p == 5 && q == 3){
        duale = true;
        unsigned int tmp = p;
        p = q;
        q = tmp;
        return "Icosaedro";
    }

    // Se i numeri non identificano un solido platonico, restituiamo "Errore" che sarà gestito dopo
    else{
        return "Errore";
    }
}


int main(int argc, char *argv[]){
    
    // Salvo l'input nelle variabili p, q, b, c
    vector<unsigned int> vec;
    vec.reserve(argc-1);
    bool duale  = false;
    if(argc == 5) //input (p q b c) per sola triangolazione
    {
        for(int i = 1; i < argc; i++){

            istringstream convert(argv[i]);
            convert >> vec[i-1]; // converto l'argomento (una stringa) in un intero     
        }
        
        unsigned int p = vec[0];
        unsigned int q = vec[1];
        unsigned int b = vec[2];
        unsigned int c = vec[3];

        string Poliedro = ClassificaPoliedro(p, q, duale);

        cout << "Poliedro: " << Poliedro << "\nDuale: " << duale << endl;

        PolyhedraMesh mesh;

        if(!ImportMesh(mesh, Poliedro)){
            cerr << "La mesh non viene riempita correttamente" << endl;
        }

        // Triangolazione del poliedro
        if (b>0 && c==0) {
            TriangolazioneI(mesh, b);
            cout << "Triangolazione di tipo I del " << Poliedro << " eseguita con b = " << b << endl;
        }
        else if (b==0 && c>0) {
            TriangolazioneI(mesh, c);
            cout << "Triangolazione di tipo I del " << Poliedro << " eseguita con c = " << c << endl;
        }
        else if (b>0 && c==b) {
            TriangolazioneII(mesh, b);
            cout << "Triangolazione di tipo II del " << Poliedro << " eseguita con b = " << b << " e c = " << c << endl;
        }
        else{
            cout << "I dati non sono validi per una triangolazione" << endl;
        }

        // Duale
        if(duale){ 
            MeshDuale(mesh);
            if(Poliedro == "Ottaedro"){
                Poliedro = "Cubo";                
            }
            else if(Poliedro == "Icosaedro"){
                Poliedro = "Dodecaedro";
            }
            cout << "Mesh duale del " << Poliedro << " creata." << endl;
        }
  
        
        //proiezione su una sfera
        ProiezioneSfera(mesh); 


        //output della mesh triangolata
        outputFile(mesh, "Cell0Ds.txt", "Cell1Ds.txt", "Cell2Ds.txt", "Cell3Ds.txt"); 
                       

        //esporto le mesh
	    UCDUtilities utilities;
	    utilities.ExportPoints("./Cell0Ds.inp", 
	    					   mesh.Cell0DsCoordinates);
	
	    utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema);

        utilities.ExportPolygons("./Cell2Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell2DsVertices);
                
        return 0;
    }


    else if(argc == 7) //input (p q b c nodo_iniziale nodo_finale) per trovare i cammini minimi
    {
        for(int i = 1; i < argc; i++){

            istringstream convert(argv[i]);
            convert >> vec[i-1]; // converto l'argomento (una stringa) in un intero     
        }
        
        unsigned int p = vec[0];
        unsigned int q = vec[1];
        unsigned int b = vec[2];
        unsigned int c = vec[3];
        unsigned int nodo_iniziale = vec[4];
        unsigned int nodo_finale = vec[5];

        string Poliedro = ClassificaPoliedro(p, q, duale);

        cout << "Poliedro: " << Poliedro << "\nDuale: " << duale << endl;

        PolyhedraMesh mesh;

        if(!ImportMesh(mesh, Poliedro)){
            cerr << "La mesh non viene riempita correttamente" << endl;
        }

        // Triangolazione del poliedro
        if (b>0 && c==0) {
            TriangolazioneI(mesh, b);
            cout << "Triangolazione di tipo I del " << Poliedro << " eseguita con b = " << b << endl;
        }
        else if (b==0 && c>0) {
            TriangolazioneI(mesh, c);
            cout << "Triangolazione di tipo I del " << Poliedro << " eseguita con c = " << c << endl;
        }
        else if (b>0 && c==b) {
            TriangolazioneII(mesh, b);
            cout << "Triangolazione di tipo II del " << Poliedro << " eseguita con b = " << b << " e c = " << c << endl;
        }
        else{
            cout << "I dati non sono validi per una triangolazione" << endl;
        }

        // Duale
        if(duale){ 
            MeshDuale(mesh);
            if(Poliedro == "Ottaedro"){
                Poliedro = "Cubo";                
            }
            else if(Poliedro == "Icosaedro"){
                Poliedro = "Dodecaedro";
            }
            cout << "Mesh duale del " << Poliedro << " creata." << endl;
        }

        //cammini minimi
        if (nodo_iniziale >= mesh.NumCell0Ds || nodo_finale >= mesh.NumCell0Ds) 
        {
            cerr << "Nodi iniziale o finale non validi." << endl;
            return 1;
        }
        else if (nodo_iniziale == nodo_finale) 
        {
            cout << "Il nodo iniziale e il nodo finale sono gli stessi. Distanza totale: 0." << endl;
            cout << "Il percorso individuato passa per il nodo: " << nodo_iniziale << endl;
            return 0;
        }
        else 
        {
            ListaAdiacenza grafo = CreaListaAdiacenza(mesh);
            auto [percorso, distanza_totale] = Dijkstra(grafo, nodo_iniziale, nodo_finale);
            if (percorso.empty()) //-1.0
            {
                cout << "Nessun percorso trovato." << endl;
            }
            else 
            {
                cout << "Il percorso individuato ha distanza totale: " << distanza_totale << " e passa per i nodi: ";
                for (unsigned int nodo : percorso)
                {
                    cout << nodo << " ";
                }
                cout << endl;
            }
            
            //proiezione su una sfera
            //ProiezioneSfera(mesh);


            //output della mesh triangolata
            outputFile(mesh, "Cell0Ds.txt", "Cell1Ds.txt", "Cell2Ds.txt", "Cell3Ds.txt"); 


            //esporto le mesh 
            
            //Proprietà 1. per i punti del cammino minimo
            vector<double> PuntiCammino(mesh.NumCell0Ds, 0.0); //cioè su paraview, nella scala tra blu (0) e rosso (1), i punti inizialmente non fanno parte del cammino minimo
            for (const auto& punto : percorso)
                PuntiCammino[punto] = 1.0; //1 ai punti che fanno parte del cammino minimo
                                
            UCDProperty<double> PropPuntiCammino;
            PropPuntiCammino.Label = "Cammino Minimo: Punti"; 
            //PropCamminoMinimo.UnitLabel = "";
            PropPuntiCammino.Size = PuntiCammino.size(); 
            PropPuntiCammino.NumComponents = 1; 
            PropPuntiCammino.Data = PuntiCammino.data(); 

            vector<UCDProperty<double>> PropPunti;
            PropPunti.push_back(PropPuntiCammino);  


            //Proprietà 2. per i segmenti del cammino minimo
            //inizialmente tutti i segmenti hanno valore 0.0 (non fanno parte del cammino minimo)
            //poi assegno 1.0 ai segmenti che fanno parte del cammino minimo
            vector<double> SegmentiCammino(mesh.Cell1DsExtrema.cols(), 0.0); //inizializzo il vettore dei segmenti del cammino minimo, che all'inizio contiene tanti segmenti quanti ce ne sono nella mesh, con 0.0 
            

            map<pair<unsigned int, unsigned int>, unsigned int> mappaSegmenti; //mappa per trovare rapidamente gli indici dei segmenti
            for (unsigned int e = 0; e < mesh.Cell1DsExtrema.cols(); e++) //for e = edge in Cell1DsExtrema, ovvero per ogni segmento nella matrice di estremi dei segmenti
            {
                //v = vertice
                unsigned int v0 = mesh.Cell1DsExtrema(0, e); //primo vertice del segmento della matrice Cell1DsExtrema  
                unsigned int v1 = mesh.Cell1DsExtrema(1, e); //secondo vertice del segmento della matrice Cell1DsExtrema
                mappaSegmenti[make_pair(v0, v1)] = e; //associo il segmento (v0, v1) all'indice e corrispondente nella mappa
            }
            //in questo modo, la mappa permette di trovare rapidamente l'indice del segmento dato una coppia di vertici (v0, v1)

            //assegno 1.0 ai segmenti che fanno parte del cammino
            for (unsigned int i = 0; i < percorso.size() - 1; i++) //size è il numero di nodi nel cammino
            {
                 //v_c = vertice cammino
                unsigned int v_c0 = percorso[i]; //nodo corrente
                unsigned int v_c1 = percorso[i+1]; //nodo successivo nel cammino
                
                //per ogni coppia di nodi consecutivi nel cammino, cerco il segmento corrispondente
                //cerco il segmento in entrambe le direzioni:
                auto it = mappaSegmenti.find(make_pair(v_c0, v_c1)); //cerco il segmento (v_c0, v_c1), ordinati in questo ordine, nella mappa
                if (it != mappaSegmenti.end()) //se il segmento (v_c0, v_c1) esiste nella mappa
                {
                    //it->second accede al valore del segmento nella mappa, che è l'indice del segmento nella matrice Cell1DsExtrema
                    SegmentiCammino[it->second] = 1.0; //se il segmento (v_c0, v_c1) esiste, assegno 1.0 al segmento corrispondente
                } 
                else 
                {
                    //provo nella direzione inversa
                    it = mappaSegmenti.find(make_pair(v_c1, v_c0)); //cerco il segmento (v_c1, v_c0) nella mappa
                    if (it != mappaSegmenti.end()) 
                    {
                        SegmentiCammino[it->second] = 1.0;
                    }
                }
            }

            UCDProperty<double> PropSegmentiCammino;
            PropSegmentiCammino.Label = "Cammino Minimo: Segmenti";  // Stesso nome per consistenza
            PropSegmentiCammino.Size = SegmentiCammino.size();
            PropSegmentiCammino.NumComponents = 1;
            PropSegmentiCammino.Data = SegmentiCammino.data();

            vector<UCDProperty<double>> PropSegmenti;
            PropSegmenti.push_back(PropSegmentiCammino);


            UCDUtilities utilities;
            utilities.ExportPoints("./Cell0Ds.inp", 
                                    mesh.Cell0DsCoordinates,
                                    PropPunti);
            
            utilities.ExportSegments("./Cell1Ds.inp",
                                    mesh.Cell0DsCoordinates,
                                    mesh.Cell1DsExtrema,
                                    PropPunti,
                                    PropSegmenti); 

            utilities.ExportPolygons("./Cell2Ds.inp",
                                    mesh.Cell0DsCoordinates,
                                    mesh.Cell2DsVertices); 
                                    
            return 0;
        }
    }

    else if(argc <= 1){

        cout << "Nessun input rilevato!" << endl;
        return 1;

    }

    else{ 
        cout << "Non sono stati inseriti input a sufficienza" << endl;
        return 1; 
    }
}


