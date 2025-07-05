//DUALE UTILS CPP
#include "Duale.hpp"
#include <iostream>
#include <fstream> //libreria per la gestione dei file
#include <sstream> //libreria per manipolare le stringhe
#include <string> //libreria per le stringhe
#include <map> //libreria per le mappe
#include <vector> //libreria per i vettori
#include <cmath> //libreria per le funzioni matematiche
#include <algorithm> //libreria per le funzioni di ricerca e ordinamento
#include <set> //Per std::set

//(la libreria Eigen dev'essere già inclusa in Duale.hpp e PolyhedraMesh.hpp)

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;

void OrderFaces(const vector<unsigned int>& unordered_faces, 
               vector<unsigned int>& ordered_faces, 
               const PolyhedraMesh& mesh) 
{
    // Controllo input vuoto
    if (unordered_faces.empty()) {
        return;
    }

    // Facce rimanenti da ordinare
    vector<unsigned int> remaining_faces = unordered_faces;

    // Inizia dalla prima faccia
    unsigned int current = remaining_faces[0];
    ordered_faces.push_back(current);
    remaining_faces.erase(remaining_faces.begin());

    // Ordina le facce rimanenti
    while (!remaining_faces.empty()) {
        const vector<unsigned int>& current_edges = mesh.Cell2DsEdges[current];
        bool found = false;
        size_t found_index = 0;

        // Cerca faccia adiacente
        for (size_t i = 0; i < remaining_faces.size(); i++) {
            const vector<unsigned int>& candidate_edges = mesh.Cell2DsEdges[remaining_faces[i]];

            // Cerca edge in comune
            for (unsigned int e1 : current_edges) {
                for (unsigned int e2 : candidate_edges) {
                    if (e1 == e2) {
                        found = true;
                        found_index = i;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }

        if (found) {
            // Aggiungi faccia adiacente
            current = remaining_faces[found_index];
            ordered_faces.push_back(current);
            remaining_faces.erase(remaining_faces.begin() + found_index);
        } else {
            // Se nessuna faccia adiacente trovata, prendi la successiva disponibile
            current = remaining_faces[0];
            ordered_faces.push_back(current);
            remaining_faces.erase(remaining_faces.begin());
        }
    }
}

namespace PolyhedraLibrary {

bool MeshDuale(PolyhedraMesh& mesh) 
{
    /////////////////////////////////////////////////////////////////////////////7

    // 1. VERIFICA DELLA MESH ORIGINALE
    if (mesh.NumCell3Ds != 1) {
        cerr << "ERRORE: La mesh deve contenere esattamente 1 cella 3D" << endl;
        return false;
    }

    unsigned int NumFaces = mesh.Cell3DsNumFaces;
    if (NumFaces == 0) {
        cerr << "ERRORE: Mesh senza facce!" << endl;
        return false;
    }

    cout << "\n=== INIZIO DUALIZZAZIONE ===" << endl;
    cout << "Facce originali: " << NumFaces << endl;
    cout << "Vertici originali: " << mesh.NumCell0Ds << endl;

    //////////////////////////////////////////////////////////////////////////////

    // 2. CREAZIONE VERTICI DUALI (BARICENTRI DELLE FACCE)
    vector<unsigned int> NewCell0DsId(NumFaces);
    MatrixXd NewCell0DsCoordinates(3, NumFaces);
    map<unsigned int, unsigned int> Faces_bar;

    for (unsigned int faceID = 0; faceID < NumFaces; faceID++) {
        if (mesh.Cell2DsNumVertices[faceID] == 0) {
            cerr << "ERRORE: Faccia " << faceID << " senza vertici!" << endl;
            return false;
        }

        Vector3d baricentro = Vector3d::Zero();
        for (const auto& id_vertex : mesh.Cell2DsVertices[faceID]) {
            if (id_vertex >= mesh.NumCell0Ds) {
                cerr << "ERRORE: Faccia " << faceID << " riferimento a vertice inesistente " 
                     << id_vertex << endl;
                return false;
            }
            baricentro += mesh.Cell0DsCoordinates.col(id_vertex);
        }
        baricentro /= mesh.Cell2DsNumVertices[faceID];

        NewCell0DsCoordinates.col(faceID) = baricentro;
        NewCell0DsId[faceID] = faceID;
        Faces_bar[faceID] = faceID;
    }

    // 3. COSTRUZIONE SPIGOLI DUALI (CONNESSIONI TRA FACCE ADIACENTI)
    map<unsigned int, vector<unsigned int>> edgeToFaces;
    for (unsigned int faceID = 0; faceID < NumFaces; faceID++) {
        for (const auto& edgeID : mesh.Cell2DsEdges[faceID]) {
            if (edgeID >= mesh.NumCell1Ds) {
                cerr << "ERRORE: Faccia " << faceID << " riferimento a spigolo inesistente " 
                     << edgeID << endl;
                return false;
            }
            edgeToFaces[edgeID].push_back(faceID);
        }
    }

    MatrixXi NewCell1DsExtrema(2, edgeToFaces.size());
    vector<unsigned int> NewCell1DsId(edgeToFaces.size());
    map<pair<unsigned int, unsigned int>, unsigned int> edgeMap;
    set<pair<unsigned int, unsigned int>> createdEdges;
    unsigned int edgeIndex = 0;

    for (const auto& entry : edgeToFaces) {
        if (entry.second.size() != 2) {
            cerr << "ATTENZIONE: Spigolo " << entry.first << " con " 
                 << entry.second.size() << " facce associate (attese 2)" << endl;
            continue;
        }

        unsigned int f1 = entry.second[0];
        unsigned int f2 = entry.second[1];
        auto edgeKey = make_pair(min(f1, f2), max(f1, f2));

        if (createdEdges.find(edgeKey) == createdEdges.end()) {
            NewCell1DsExtrema(0, edgeIndex) = f1;
            NewCell1DsExtrema(1, edgeIndex) = f2;
            NewCell1DsId[edgeIndex] = edgeIndex;
            edgeMap[edgeKey] = edgeIndex;
            createdEdges.insert(edgeKey);
            edgeIndex++;
        }
    }

    cout << "Spigoli duali creati: " << edgeIndex << "/" << edgeToFaces.size() << endl;

    // 4. COSTRUZIONE FACCE DUALI (PER OGNI VERTICE ORIGINALE)
    vector<unsigned int> NewCell2DsId;
    vector<vector<unsigned int>> NewCell2DsVertices;
    vector<vector<unsigned int>> NewCell2DsEdges;
    map<unsigned int, vector<unsigned int>> vertex_to_faces;
    set<unsigned int> problematic_vertices;
    map<vector<unsigned int>, unsigned int> face_check;

    // Costruisci mappa vertice->facce
    for (unsigned int f = 0; f < mesh.Cell2DsVertices.size(); ++f) {
        for (unsigned int v : mesh.Cell2DsVertices[f]) {
            vertex_to_faces[v].push_back(f);
        }
    }

    // Crea facce duali 
    for (const auto& [v_idx, face_ids] : vertex_to_faces) {
        vector<unsigned int> orderedFaces;
        OrderFaces(face_ids, orderedFaces, mesh);

        // Crea la faccia duale
        unsigned int valence = orderedFaces.size();
        vector<unsigned int> dualFaceEdges(valence);

        // Collega gli spigoli con verifica
        for (unsigned int k = 0; k < valence; k++) {
            unsigned int curr = orderedFaces[k];
            unsigned int next = orderedFaces[(k+1)%valence];
            auto edgeKey = make_pair(min(curr, next), max(curr, next));

            if (edgeMap.find(edgeKey) == edgeMap.end()) {
                cerr << "ERRORE: Spigolo mancante tra facce " << curr << " e " << next << endl;
                return false;
            }
            dualFaceEdges[k] = edgeMap[edgeKey];
        }

        // Controllo duplicati
        vector<unsigned int> sorted_face = orderedFaces;
        sort(sorted_face.begin(), sorted_face.end());
        
        if (face_check.count(sorted_face)) {
            cerr << "FACCIA DUPLICATA: Vertice " << v_idx << " uguale a " 
                 << face_check[sorted_face] << endl;
            continue;
        }
        face_check[sorted_face] = v_idx;

        NewCell2DsId.push_back(v_idx);
        NewCell2DsVertices.push_back(orderedFaces);
        NewCell2DsEdges.push_back(dualFaceEdges);
    }

    ///////////////////////////////////////////////////

    // 5. VERIFICHE FINALI
    cout << "\n=== VERIFICHE FINALI ===" << endl;
    cout << "Facce duali create: " << NewCell2DsId.size() << endl;
    cout << "Vertici problematici: " << problematic_vertices.size() << endl;    

    /////////////////////////////////////////////////////

    // 6. AGGIORNAMENTO MESH DUALE
    mesh.NumCell0Ds = NumFaces;
    mesh.Cell0DsId = NewCell0DsId;
    mesh.Cell0DsCoordinates = NewCell0DsCoordinates;

    mesh.NumCell1Ds = edgeIndex;
    mesh.Cell1DsId = NewCell1DsId;
    mesh.Cell1DsExtrema = NewCell1DsExtrema;

    mesh.NumCell2Ds = NewCell2DsId.size();
    mesh.Cell2DsId = NewCell2DsId;
    mesh.Cell2DsVertices = NewCell2DsVertices;
    mesh.Cell2DsEdges = NewCell2DsEdges;
    
    // Imposta numero di vertici e spigoli per faccia
    mesh.Cell2DsNumVertices.resize(mesh.NumCell2Ds);
    mesh.Cell2DsNumEdges.resize(mesh.NumCell2Ds);
    for (unsigned int i = 0; i < mesh.NumCell2Ds; i++) {
        mesh.Cell2DsNumVertices[i] = mesh.Cell2DsVertices[i].size();
        mesh.Cell2DsNumEdges[i] = mesh.Cell2DsEdges[i].size();
    }

    // Aggiorna cella 3D
    mesh.NumCell3Ds = 1;
    mesh.Cell3DsId = {0};
    mesh.Cell3DsNumVertices = {mesh.NumCell0Ds};
    mesh.Cell3DsNumEdges = {mesh.NumCell1Ds};
    mesh.Cell3DsNumFaces = {mesh.NumCell2Ds};
    mesh.Cell3DsVertices = {mesh.Cell0DsId};
    mesh.Cell3DsEdges = {mesh.Cell1DsId};
    mesh.Cell3DsFaces = {mesh.Cell2DsId};

    cout << "Dualizzazione completata con successo" << endl;
    return true;
}

} // namespace PolyhedraLibrary
