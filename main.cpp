// MAIN POLYHEDRA

#include <iostream>
#include <sstream>
#include <string>
#include "PolyhedraMesh.hpp"
#include "Utils.hpp"
#include "triangolazione.hpp"
#include "Duale.hpp"
#include "output.hpp"
// #include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;


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
    if(argc > 3)
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

        outputFile(mesh, "Cell0Ds.txt", "Cell1Ds.txt", "Cell2Ds.txt", "Cell3Ds.txt"); 

        return 0;

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


