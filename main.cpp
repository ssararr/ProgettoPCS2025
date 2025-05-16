// MAIN POLYHEDRA

#include <iostream>
#include <sstream>
#include <string>
// #include "PolyhedraMesh.hpp"
// #include "Utils.hpp"
// #include "UCDUtilities.hpp"

using namespace std;
//using namespace Eigen;
//using namespace PolyhedraLibrary;


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
    unsigned int arr[argc-1] = {};
    bool duale  = false;
    if(argc > 1)
    {
        for(unsigned int i = 1; i < argc; i++){

            istringstream convert(argv[i]);
            convert >> arr[i-1]; // converto l'argomento (una stringa) in un intero
        
        }
        
    }

    unsigned int p = arr[0];
    unsigned int q = arr[1];
    unsigned int b = arr[2];
    unsigned int c = arr[3];

    string Poliedro = ClassificaPoliedro(p, q, duale);

    cout << "Poliedro: " << Poliedro << "\nDuale: " << duale << endl;



    // Identifico il poliedro



}


