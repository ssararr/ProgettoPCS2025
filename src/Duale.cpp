//DUALE UTILS CPP
#include "Duale.hpp"
#include <iostream>
#include <fstream> //libreria per la gestione dei file
#include <sstream> //libreria per manipolare le stringhe
#include <string> //libreria per le stringhe
#include <map> //libreria per le mappe
#include <vector> //libreria per i vettori
#include <cmath> //libreria per le funzioni matematiche
#include <limits> //libreria per i limiti numerici
#include <algorithm> //libreria per le funzioni di ricerca e ordinamento
#include <numeric>  // Per std::iota
#include <set> //Per std::set
#include <queue> //Per std::queue

//(la libreria Eigen dev'essere già inclusa in Duale.hpp e PolyhedraMesh.hpp)

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;

void OrderFaces(const vector<unsigned int>& unordered_faces,
               vector<unsigned int>& ordered_faces,
               const PolyhedraMesh& mesh)
{
    ordered_faces.clear();
    
    // Caso base: nessuna faccia o una sola faccia
    if (unordered_faces.empty()) {
        return;
    }
    if (unordered_faces.size() == 1) {
        ordered_faces.push_back(unordered_faces[0]);
        return;
    }

    // Mappa per trovare rapidamente le facce adiacenti a ciascuna faccia
    map<unsigned int, vector<unsigned int>> face_adjacency;
    
    // Costruisci la mappa di adiacenza tra facce
    for (size_t i = 0; i < unordered_faces.size(); ++i) {
        unsigned int face1 = unordered_faces[i];
        const auto& edges1 = mesh.Cell2DsEdges[face1];
        
        for (size_t j = i+1; j < unordered_faces.size(); ++j) {
            unsigned int face2 = unordered_faces[j];
            const auto& edges2 = mesh.Cell2DsEdges[face2];
            
            // Cerca un edge in comune
            bool found_common = false;
            for (unsigned int e1 : edges1) {
                for (unsigned int e2 : edges2) {
                    if (e1 == e2) {
                        face_adjacency[face1].push_back(face2);
                        face_adjacency[face2].push_back(face1);
                        found_common = true;
                        break;
                    }
                }
                if (found_common) break;
            }
        }
    }

    // Verifica che tutte le facce siano connesse
    set<unsigned int> visited;
    queue<unsigned int> to_visit;
    to_visit.push(unordered_faces[0]);
    visited.insert(unordered_faces[0]);

    while (!to_visit.empty()) {
        unsigned int current = to_visit.front();
        to_visit.pop();
        
        for (unsigned int neighbor : face_adjacency[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                to_visit.push(neighbor);
            }
        }
    }

    // Se non tutte le facce sono connesse, usiamo un ordinamento semplice
    if (visited.size() != unordered_faces.size()) {
        ordered_faces = unordered_faces;
        return;
    }

    // Ordinamento ciclico delle facce
    ordered_faces.push_back(unordered_faces[0]);
    unsigned int prev_face = unordered_faces[0];
    unsigned int first_face = prev_face;
    
    while (ordered_faces.size() < unordered_faces.size()) {
        bool found_next = false;
        
        // Cerca la prossima faccia non ancora usata che è adiacente alla precedente
        for (unsigned int candidate : face_adjacency[prev_face]) {
            if (find(ordered_faces.begin(), ordered_faces.end(), candidate) == ordered_faces.end()) {
                ordered_faces.push_back(candidate);
                prev_face = candidate;
                found_next = true;
                break;
            }
        }
        
        // Se non troviamo una faccia adiacente, cerchiamo qualsiasi faccia rimanente
        if (!found_next) {
            for (unsigned int remaining : unordered_faces) {
                if (find(ordered_faces.begin(), ordered_faces.end(), remaining) == ordered_faces.end()) {
                    ordered_faces.push_back(remaining);
                    prev_face = remaining;
                    break;
                }
            }
        }
        
        // Controllo di sicurezza per evitare loop infiniti
        if (ordered_faces.size() == unordered_faces.size()) {
            break;
        }
    }

    // Verifica che l'ultima faccia sia adiacente alla prima per chiudere il ciclo
    bool last_connected_to_first = false;
    for (unsigned int neighbor : face_adjacency[ordered_faces.back()]) {
        if (neighbor == first_face) {
            last_connected_to_first = true;
            break;
        }
    }
    
    if (!last_connected_to_first) {
        // Se non sono connesse, riordiniamo per trovare una sequenza migliore
        vector<unsigned int> reordered;
        reordered.push_back(first_face);
        
        for (size_t i = 1; i < ordered_faces.size(); ++i) {
            bool found_better = false;
            
            for (size_t j = i; j < ordered_faces.size(); ++j) {
                for (unsigned int neighbor : face_adjacency[reordered.back()]) {
                    if (neighbor == ordered_faces[j]) {
                        reordered.push_back(ordered_faces[j]);
                        ordered_faces.erase(ordered_faces.begin() + j);
                        found_better = true;
                        break;
                    }
                }
                if (found_better) break;
            }
            
            if (!found_better && i < ordered_faces.size()) {
                reordered.push_back(ordered_faces[i]);
            }
        }
        
        ordered_faces = reordered;
    }
}



//definisco una funzione distanza che calcola la distanza euclidea tra due vettori di dimensione 3 che contengono le cooridnate x, y, z
double distanceBtw(const vector<double>& point1, const vector<double>& point2) {
	// Calcola la distanza euclidea tra due punti 3D rappresentati come vettori di coordinate
	if (point1.size() != 3 || point2.size() != 3) {
		throw invalid_argument("Entrambi i punti devono essere vettori di dimensione 3.");
	}
	double dx = point1[0] - point2[0];
	double dy = point1[1] - point2[1];
	double dz = point1[2] - point2[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

//AGGIUNTI CONTROLLI SULLE DIMENSIONI DEI VECTOR<VECTOR OTTENUTI (VOGLIAMO CHE SIANO 3D) E SULLA DIMENSIONE DELLA MATRICE (VOGLIAMO CHE ABBIA 3 COLONNE, ALTRIMENTI NON SONO COORDINATE 3D))
vector<vector<double>> eigenMatrixToVectorVector(const MatrixXd& matrix) {
    // Controllo che la matrice abbia esattamente 3 righe (x,y,z)
    if (matrix.rows() != 3) {
        throw invalid_argument("La matrice deve avere 3 righe (coordinate 3D). Ricevuto: " + 
                              to_string(matrix.rows()) + " righe.");
    }

    const int numPoints = matrix.cols(); // Numero di punti (colonne)
    vector<vector<double>> result;
    result.reserve(numPoints);

    for (int i = 0; i < numPoints; ++i) {
        // Estrae direttamente la colonna i-esima come vector<double>
        vector<double> point = {
            matrix(0, i), // x
            matrix(1, i), // y
            matrix(2, i)  // z
        };
        result.push_back(point);
    }
	return result;
}

map<int, vector<unsigned int>> findClosestDualVertices(
    const vector<vector<double>>& originalVertices,
    const vector<vector<double>>& dualVertices,
    unsigned int num_closest = 4) // Default a 4 (per cubo/ottaedro), passare 5 per dodecaedro/icosaedro
{
    map<int, vector<unsigned int>> vertexToDualMap;

    // Verifica preliminare
    if (num_closest == 0) {
        throw invalid_argument("num_closest deve essere almeno 1");
    }

    // Controlla che tutti i vertici siano 3D
    for (const auto& v : originalVertices) {
        if (v.size() != 3) throw invalid_argument("Vertici originali devono essere 3D");
    }
    for (const auto& v : dualVertices) {
        if (v.size() != 3) throw invalid_argument("Vertici duali devono essere 3D");
    }

    // Prealloca memoria per le distanze
    vector<pair<double, unsigned int>> distances;
    distances.reserve(dualVertices.size());

    for (unsigned int i = 0; i < originalVertices.size(); ++i) {
        const auto& original = originalVertices[i];
        distances.clear();

        // Calcola tutte le distanze usando distanceBtw
        for (unsigned int j = 0; j < dualVertices.size(); ++j) {
            const auto& dual = dualVertices[j];
            double dist = distanceBtw(original, dual);
            distances.emplace_back(dist, j);
        }

        // Seleziona i num_closest più vicini
        unsigned int n = min(num_closest, static_cast<unsigned int>(distances.size()));
        
        // Ordina parzialmente per ottenere i primi n elementi
        std::partial_sort(
            distances.begin(), 
            distances.begin() + n, 
            distances.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; }
        );

        // Estrai gli indici
        vector<unsigned int> closestIndices;
        closestIndices.reserve(n);
        for (unsigned int k = 0; k < n; ++k) {
            closestIndices.push_back(distances[k].second);
        }

        vertexToDualMap[i] = closestIndices;
    }

    return vertexToDualMap;
}

/*map<int, vector<unsigned int>> findClosestDualVertices(
    const vector<vector<double>>& originalVertices,
    const vector<vector<double>>& dualVertices) 
{
	double tolerance = 1e-10; // Definisco una tolleranza per considerare due distanze equivalenti (può essere modificata in base alle esigenze)
    map<int, vector<unsigned int>> vertexToDualMap;

	for(auto& original : originalVertices){

	}

    for (unsigned int i = 0; i < originalVertices.size(); ++i) {
        const auto& original = originalVertices[i]; // [x,y,z]
        double minDistance = numeric_limits<double>::max();
        vector<unsigned int> closestDualIds;

        for (unsigned int j = 0; j < dualVertices.size(); ++j) {
            const auto& dual = dualVertices[j]; // [x,y,z]
            
            double dist = distanceBtw(original, dual);

            if (dist < minDistance - tolerance) {
                minDistance = dist;
                closestDualIds = {j}; // Sostituisci con nuovo minimo
            } 
            else if (abs(dist - minDistance) <= tolerance) {
                closestDualIds.push_back(j); // Aggiungi a minimi esistenti
            }
        }

        vertexToDualMap[i] = closestDualIds;
    }

    return vertexToDualMap;
}*/


namespace PolyhedraLibrary {

bool MeshDuale(PolyhedraMesh& mesh) 
{
    // 1. INIZIALIZZAZIONE E VERIFICA DELLA MESH ORIGINALE
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

    // Verifica vertici senza facce
    for (unsigned int v = 0; v < mesh.NumCell0Ds; ++v) {
        if (vertex_to_faces[v].empty()) {
            cerr << "ATTENZIONE: Vertice " << v << " non ha facce associate!" << endl;
        }
    }

    // Crea facce duali con controlli rigorosi
    for (const auto& [v_idx, face_ids] : vertex_to_faces) {
        vector<unsigned int> orderedFaces;
        OrderFaces(face_ids, orderedFaces, mesh);

            // Semplifichiamo prendendo solo le prime 4 facce
           /* if (orderedFaces.size() > 6) {
                orderedFaces.resize(4); // Forza a quadrilatero
                cerr << "  Ridotto a quadrilatero" << endl;
            }
        
        // Controllo integrità ordinamento
        if (orderedFaces.size() != face_ids.size()) {
            cerr << "PROBLEMA: Vertice " << v_idx << " - facce ordinate " 
                 << orderedFaces.size() << " vs originali " << face_ids.size() << endl;
            problematic_vertices.insert(v_idx);
            continue;
        }*/

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

    // 5. VERIFICHE FINALI
    cout << "\n=== VERIFICHE FINALI ===" << endl;
    cout << "Facce duali create: " << NewCell2DsId.size() << endl;
    cout << "Vertici problematici: " << problematic_vertices.size() << endl;
    
    if (NewCell2DsId.size() != mesh.NumCell0Ds) {
        cerr << "DISCREPANZA: " << mesh.NumCell0Ds << " vertici ma " 
             << NewCell2DsId.size() << " facce duali create" << endl;
        
        // Trova vertici mancanti
        set<unsigned int> dual_vertices(NewCell2DsId.begin(), NewCell2DsId.end());
        for (unsigned int v = 0; v < mesh.NumCell0Ds; ++v) {
            if (dual_vertices.find(v) == dual_vertices.end()) {
                cerr << "Vertice " << v << " senza faccia duale" << endl;
            }
        }
    }

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

    /////////////////////////////////////////////////////////////////////////////////////

    // Verifica DOPO la creazione della mesh duale ma PRIMA dell'esportazione
    /*cout << "\n=== ANALISI FACCE DUALI PRIMA ESPORTAZIONE ===" << endl;

    // Conta i tipi di facce
    map<size_t, size_t> face_types;
    for (const auto& face : mesh.Cell2DsVertices) {
    face_types[face.size()]++;
    }

    // Stampa la distribuzione
    for (const auto& [vertex_count, num_faces] : face_types) {
    cout << "Facce con " << vertex_count << " vertici: " << num_faces << endl;
    if (vertex_count < 3) {
        cerr << "ERRORE: Trovata faccia con meno di 3 vertici!" << endl;
    }
    }

    // Verifica esplicita delle facce non triangolari/quadrangolari
    vector<size_t> problematic_faces;
    for (size_t i = 0; i < mesh.Cell2DsVertices.size(); ++i) {
    const auto& face = mesh.Cell2DsVertices[i];
    if (face.size() != 3 && face.size() != 4) {
        problematic_faces.push_back(i);
    }
    }

    if (!problematic_faces.empty()) {
    cerr << "\nATTENZIONE: Trovate " << problematic_faces.size() 
            << " facce non triangolari/quadrangolari" << endl;
    cerr << "Indici facce problematiche:";
    for (size_t i : problematic_faces) cerr << " " << i;
    cerr << endl;
    
    // Stampa i dettagli delle prime 5 facce problematiche
    cerr << "\nDettagli prime 5 facce problematiche:" << endl;
    for (size_t i = 0; i < min(problematic_faces.size(), size_t(5)); ++i) {
        size_t idx = problematic_faces[i];
        cerr << "Faccia " << idx << " (" << mesh.Cell2DsVertices[idx].size() 
                << " vertici):";
        for (unsigned int v : mesh.Cell2DsVertices[idx]) cerr << " " << v;
        cerr << endl;
    }
}*/

///////////////////////////////////////////////////////////////////////////////////////////////////

    cout << "Dualizzazione completata con successo" << endl;
    return true;
}

} // namespace PolyhedraLibrary
