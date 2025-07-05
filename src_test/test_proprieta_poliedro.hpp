// POLYHEDRA TEST

#pragma once

#include <iostream>
#include <map>
#include <string>
#include <gtest/gtest.h>
#include <algorithm>
#include "PolyhedraMesh.hpp"
#include "Utils.hpp"
#include "triangolazione.hpp"
#include "Duale.hpp"
#include "CamminiMinimi.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;


unsigned int VerticiAttesi(unsigned int b, const PolyhedraMesh& Poliedro, bool TipoII){
    if(!TipoII){
        unsigned int T = b*b;
        switch(Poliedro.NumCell0Ds){
            case(4):
               return 2*T+2;
               break;
            case(6):
                return 4*T+2;
                break;
            case(12):
                return 10*T+2;
                break; 
        }
    }
    else{
        return Poliedro.NumCell0Ds + Poliedro.NumCell1Ds*(2*b-1) + Poliedro.NumCell2Ds * ((3.0*b*b)/2.0 - (3.0*b/2.0) + 1); // Formula per vertici in triangolazione II
    }
    return 0;
}
unsigned int EdgesAttesi(unsigned int b, const PolyhedraMesh& Poliedro, bool TipoII){
    if(!TipoII){
        unsigned int T = b*b;
        switch(Poliedro.NumCell0Ds){
            case(4):
               return 6*T;
               break;
            case(6):
                return 12*T;
                break;
            case(12):
                return 30*T;
                break; 
        }
    }
    else{
        return Poliedro.NumCell1Ds * 2 * b + Poliedro.NumCell2Ds * ((9.0*b*b)/2.0 + (3.0*b/2.0)); // Formula per edges in triangolazione II
    }
    return 0;
}
unsigned int FacceAttese(unsigned int b, const PolyhedraMesh& Poliedro, bool TipoII){
    if(!TipoII){
        unsigned int T = b*b;
        switch(Poliedro.NumCell0Ds){
            case(4):
               return 4*T;
               break;
            case(6):
                return 8*T;
                break;
            case(12):
                return 20*T;
                break; 
        }
    }
    else{
        return Poliedro.NumCell2Ds * (3*b*b + 3*b);
    }
    return 0;
}

// Restituisco il numero di vertici ad una data valenza, solo per TriangolazioneI perché i dati sono errati per TriangolazioneII
unsigned int NumVerticiConVal(const PolyhedraMesh& Poliedro, unsigned int Valenza) {
    unsigned int NumVertici = 0;
    for(const auto& id_vertice : Poliedro.Cell0DsId) {
        unsigned int Valenza_attuale = 0;
        for(const auto& faccia : Poliedro.Cell2DsVertices){
            if(find(faccia.begin(), faccia.end(), id_vertice) != faccia.end())
                Valenza_attuale++;
        }
        if(Valenza_attuale == Valenza) 
            NumVertici++;
    }

    return NumVertici;
}


namespace TestPolyhedraLibrary {

// Test sulle triangolazioni

    TEST(TestTriangolazione, TestTetraedroI)
    {
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Tetraedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 5;
        unsigned int V = VerticiAttesi(b, Originale, false);
        unsigned int E = EdgesAttesi(b, Originale, false);
        unsigned int F = FacceAttese(b, Originale, false);

        TriangolazioneI(Originale, b);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(E, Originale.NumCell1Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);
        
        // Verifico la valenza dei vertici
        unsigned int Valenza1 = 3;
        unsigned int Valenza2 = 6;
        EXPECT_EQ(8, NumVerticiConVal(Originale, Valenza1));
        EXPECT_EQ(2*(b*b-1), NumVerticiConVal(Originale, Valenza2));
    }

    TEST(TestTriangolazione, TestTetraedroII)
    {
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Tetraedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int V = VerticiAttesi(b, Originale, true);
        unsigned int E = EdgesAttesi(b, Originale, true);
        unsigned int F = FacceAttese(b, Originale, true);

        TriangolazioneII(Originale, b);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(E, Originale.NumCell1Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);

        // Non verifico valenza perché non rispetta le formule
    }

    TEST(TestTriangolazione, TestOttaedroI)
    {
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Ottaedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int V = VerticiAttesi(b, Originale, false);
        unsigned int E = EdgesAttesi(b, Originale, false);
        unsigned int F = FacceAttese(b, Originale, false);

        TriangolazioneI(Originale, b);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(E, Originale.NumCell1Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);
        
        // Verifico la valenza dei vertici
        EXPECT_EQ(12, NumVerticiConVal(Originale, 4));
        EXPECT_EQ(4*(b*b-1), NumVerticiConVal(Originale, 6));
    }

    TEST(TestTriangolazione, TestOttaedroII){
        
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Ottaedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int V = VerticiAttesi(b, Originale, true);
        unsigned int E = EdgesAttesi(b, Originale, true);
        unsigned int F = FacceAttese(b, Originale, true);

        TriangolazioneII(Originale, b);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(E, Originale.NumCell1Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);
    }

    TEST(TestTriangolazione, TestIcosaedroI){
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Icosaedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int V = VerticiAttesi(b, Originale, false);
        unsigned int E = EdgesAttesi(b, Originale, false);
        unsigned int F = FacceAttese(b, Originale, false);

        TriangolazioneI(Originale, b);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(E, Originale.NumCell1Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);
        
        // Verifico la valenza dei vertici
        EXPECT_EQ(24, NumVerticiConVal(Originale, 5));
        EXPECT_EQ(10*(b*b-1), NumVerticiConVal(Originale, 6));
    }

    TEST(TestTriangolazione, TestIcosaedroII){
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Icosaedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int V = VerticiAttesi(b, Originale, true);
        unsigned int E = EdgesAttesi(b, Originale, true);
        unsigned int F = FacceAttese(b, Originale, true);

        TriangolazioneII(Originale, b);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(E, Originale.NumCell1Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);
    }

// Test Duale

    TEST(TestDuale, TestSuTriangolazioneI){
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Tetraedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int F = VerticiAttesi(b, Originale, false);    // Scambio facce e vertici perché è il duale
        unsigned int V = FacceAttese(b, Originale, false);

        TriangolazioneI(Originale, b);
        MeshDuale(Originale);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);
    }

    TEST(TestDuale, TestSuTriangolazioneII){
        PolyhedraMesh Originale;
	    if (!ImportMesh(Originale, "Tetraedro"))
		    FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int F = VerticiAttesi(b, Originale, true);    // Scambio facce e vertici perché è il duale
        unsigned int V = FacceAttese(b, Originale, true);

        TriangolazioneII(Originale, b);
        MeshDuale(Originale);
        
        EXPECT_EQ(V, Originale.NumCell0Ds);
        EXPECT_EQ(F, Originale.NumCell2Ds);
    }

// Test Cammini Minimi

    TEST(TestCamminiMinimi, TestSuTriangolazioneI){
        PolyhedraMesh Mesh;
        if(!ImportMesh(Mesh, "Icosaedro"))
            FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int nodo_iniziale = 76;
        unsigned int nodo_finale = 81;

        TriangolazioneI(Mesh, b);

        vector<unsigned int> PercorsoAtteso = {76, 86, 82, 81};
        double DistanzaAttesa = 0;
        for(unsigned int i = 1; i < PercorsoAtteso.size(); i++){
            double dist = (Mesh.Cell0DsCoordinates.col(PercorsoAtteso[i]) - Mesh.Cell0DsCoordinates.col(PercorsoAtteso[i-1])).norm();
            DistanzaAttesa += dist;
        }

        // Creo i cammini minimi
        ListaAdiacenza grafo = CreaListaAdiacenza(Mesh);
        auto [percorso, distanza_totale] = Dijkstra(grafo, nodo_iniziale, nodo_finale);

        EXPECT_EQ(percorso, PercorsoAtteso);
        EXPECT_NEAR(distanza_totale, DistanzaAttesa, 1e-6);
    }

    TEST(TestCamminiMinimi, TestSuTriangolazioneII){

        PolyhedraMesh Mesh;
        if(!ImportMesh(Mesh, "Icosaedro"))
            FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int nodo_iniziale = 76;
        unsigned int nodo_finale = 81;

        TriangolazioneII(Mesh, b);

        vector<unsigned int> PercorsoAtteso = {76, 367, 86, 358, 353, 81};
        double DistanzaAttesa = 0;
        for(unsigned int i = 1; i < PercorsoAtteso.size(); i++){
            double dist = (Mesh.Cell0DsCoordinates.col(PercorsoAtteso[i]) - Mesh.Cell0DsCoordinates.col(PercorsoAtteso[i-1])).norm();
            DistanzaAttesa += dist;
        }

        // Creo i cammini minimi
        ListaAdiacenza grafo = CreaListaAdiacenza(Mesh);
        auto [percorso, distanza_totale] = Dijkstra(grafo, nodo_iniziale, nodo_finale);

        EXPECT_EQ(percorso, PercorsoAtteso);
        EXPECT_NEAR(distanza_totale, DistanzaAttesa, 1e-6);
        
    }

    TEST(TestCamminiMinimi, TestSuDuale){
        PolyhedraMesh Mesh;
        if(!ImportMesh(Mesh, "Ottaedro"))
            FAIL() << "Errore nel riempimento della mesh";

        unsigned int b = 4;
        unsigned int nodo_iniziale = 12;
        unsigned int nodo_finale = 25;

        TriangolazioneI(Mesh, b);
        MeshDuale(Mesh);

        vector<unsigned int> PercorsoAtteso = {12, 9, 7, 2, 1, 23, 25};
        double DistanzaAttesa = 0;
        for(unsigned int i = 1; i < PercorsoAtteso.size(); i++){
            double dist = (Mesh.Cell0DsCoordinates.col(PercorsoAtteso[i]) - Mesh.Cell0DsCoordinates.col(PercorsoAtteso[i-1])).norm();
            DistanzaAttesa += dist;
        }

        // Creo i cammini minimi
        ListaAdiacenza grafo = CreaListaAdiacenza(Mesh);
        auto [percorso, distanza_totale] = Dijkstra(grafo, nodo_iniziale, nodo_finale);

        EXPECT_EQ(percorso, PercorsoAtteso);
        EXPECT_NEAR(distanza_totale, DistanzaAttesa, 1e-6);
    }

}
