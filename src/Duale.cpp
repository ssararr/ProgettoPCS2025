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

//(la libreria Eigen dev'essere già inclusa in Duale.hpp e PolyhedraMesh.hpp)

using namespace std;
using namespace Eigen;
using namespace PolyhedraLibrary;


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
	//ASSOCIO OGNI FACCIA AL SUO VERTICE DUALE (RIEMPIO NewCell0DsId e NewCell0DsCoordinates)
	unsigned int NumFaces = mesh.Cell3DsNumFaces;
	if (NumFaces == 0) {
		cerr << "Non sono presenti facce nella mesh." << endl;
		return false;
	}
	vector<unsigned int> NewCell0DsId;
	MatrixXd NewCell0DsCoordinates(3, NumFaces); // Inizializzo la matrice delle coordinate dei vertici duali

	for (unsigned int faceID = 0; faceID < NumFaces; faceID++) {
		if (mesh.Cell2DsNumVertices[faceID] == 0) {
			cerr << "La faccia faceID " << faceID << " non ha vertici." << endl;
			return false;
		}
		Vector3d baricentro = Vector3d::Zero(); // Inizializzo il baricentro a zero
		for (const auto& id_vertex : mesh.Cell2DsVertices[faceID]) {	// [[0,]]
			baricentro += mesh.Cell0DsCoordinates.col(id_vertex); // Sommo le coordinate dei vertici della faccia
		}		
		baricentro /= mesh.Cell2DsNumVertices[faceID]; // Calcolo il baricentro dividendo per il numero di vertici della faccia
		NewCell0DsCoordinates.col(faceID) = baricentro; // Assegno il baricentro alla colonna corrispondente nella matrice delle coordinate dei vertici duali
		NewCell0DsId.push_back(faceID); // Associo l'ID della faccia al nuovo vertice duale
		//(inutile?) dualVertexToFace[faceID] = baricentro; // Associo il baricentro alla faccia
	}
	
	
	//UNISCO OGNI BARICENTRO (VERTICI DUALI) AI BARICENTRI DELLE FACCE ADIACENTI (RIEMPIO NewCell1DsId e NewCell1DsExtrema)
	vector<pair<unsigned int, unsigned int>> dualEdges; // Vettore di coppie di ID delle facce che condividono un edge 

	map<unsigned int, vector<unsigned int>> edgeToFaces; // Mappa che associa ad ogni spigolo le facce che condividono tale spigolo

	//Per ogni faccia, scorro i suoi spigoli e aggiungo tale faccia nella lista delle facce che condividono tale spigolo
	for (unsigned int faceID = 0; faceID < NumFaces; faceID++) {		
		for (const auto& edgeID : mesh.Cell2DsEdges[faceID]) {
			edgeToFaces[edgeID].push_back(faceID);
		}
	}
	//Per ogni spigolo condiviso da esattamente DUE facce, aggiungo un dualEdge
	for (const auto& entry : edgeToFaces) {
		const auto& facesSharingEdge = entry.second; //accedo alla lista di facce che condividono tale edge (chiave della mappa) 
		if (facesSharingEdge.size() == 2) { // Se lo spigolo è condiviso da due facce
			dualEdges.emplace_back(min(facesSharingEdge[0], facesSharingEdge[1]), max(facesSharingEdge[0], facesSharingEdge[1])); // Aggiungo la coppia di ID delle facce come un nuovo spigolo duale
		}
	}
	
	MatrixXi NewCell1DsExtrema(2, dualEdges.size()); // Inizializzo la matrice degli estremi degli spigoli duali
	vector<unsigned int> NewCell1DsId(dualEdges.size()); // Vettore per gli ID degli spigoli duali
	unsigned int edgeIndex = 0;
	for (auto& entry : edgeToFaces)
	{
		auto& edge = entry.second;
		NewCell1DsExtrema(0, edgeIndex) = edge[0]; // Assegno il primo ID della coppia come primo estremo dello spigolo duale
		NewCell1DsExtrema(1, edgeIndex) = edge[1]; // Assegno il secondo ID della coppia come secondo estremo dello spigolo duale
		NewCell1DsId.push_back(edgeIndex); // Aggiungo l'ID dello spigolo duale
		edgeIndex++;
	}


	//DEFINISCO OGNI FACCIA DUALE COME L'INSIEME DEI VERTICI DUALI PIU' VICINI AD OGNI VERTICE ORIGINARIO
	vector<vector<double>> originalVertices = eigenMatrixToVectorVector(mesh.Cell0DsCoordinates); // Converto le coordinate dei vertici originali in un vettore di vettori
	vector<vector<double>> dualVertices = eigenMatrixToVectorVector(NewCell0DsCoordinates); // Converto le coordinate dei vertici duali in un vettore di vettori

	map<int, vector<unsigned int>> mappaVertexToDual;	

	//AGGIUNTI CONTROLLI SULLE DIMENSIONI DELLE mesh.Cell0DsCoordinates E NewCell0DsCoordinates PRIMA DEL RICHIAMO DI findClosestDualVertices
	cout << "Dimensioni Cell0DsId: " << mesh.Cell0DsId.size() << endl;
	cout << "Dimensioni NewCell0DsId: " << NewCell0DsId.size() << endl;
	cout << "Dimensioni Cell0DsCoordinates: " << mesh.Cell0DsCoordinates.rows() << "x" << mesh.Cell0DsCoordinates.cols() << endl;
    cout << "Dimensioni NewCell0DsCoordinates: " << NewCell0DsCoordinates.rows() << "x" << NewCell0DsCoordinates.cols() << endl;


	unsigned int num_closest = 4;
	mappaVertexToDual = findClosestDualVertices(originalVertices, dualVertices, num_closest); // Trovo i vertici duali più vicini per ogni vertice originale /*

	vector<vector<unsigned int>> NewCell2DsVertices;
	vector<unsigned int> NewCell2DsId;

	// 1. Costruisci mappa vertex → facce che lo contengono
	map<unsigned int, vector<unsigned int>> vertex_to_faces;

	for (unsigned int f = 0; f < mesh.Cell2DsVertices.size(); ++f) {
		for (unsigned int v : mesh.Cell2DsVertices[f]) {
			vertex_to_faces[v].push_back(f);
		}
	}

	// 2. Ordina ciclicamente i centroidi delle facce attorno a ciascun vertice
	for (const auto& [v_idx, face_ids] : vertex_to_faces) {
		const auto& origin = originalVertices[v_idx];
		vector<pair<double, unsigned int>> angle_face_pairs;

		// Calcolo normale media per costruire un piano locale
		vector<double> normal(3, 0.0);

		for (unsigned int f_id : face_ids) {
			const auto& centroid = dualVertices[f_id];
			normal[0] += (centroid[1] - origin[1]) * (centroid[2] - origin[2]);
			normal[1] += (centroid[2] - origin[2]) * (centroid[0] - origin[0]);
			normal[2] += (centroid[0] - origin[0]) * (centroid[1] - origin[1]);
		}

		// Vettori base locale ortogonali
		vector<double> e1 = {1, 0, 0};
		vector<double> e2 = {0, 1, 0};
		if (normal[0] != 0 || normal[1] != 0 || normal[2] != 0) {
			// ortonormalizzazione (semplificata)
			double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
			for (int i = 0; i < 3; ++i) normal[i] /= norm;
			e1 = {-normal[1], normal[0], 0}; // un vettore ortogonale (arbitrario)
			double len = sqrt(e1[0]*e1[0] + e1[1]*e1[1]);
			if (len > 1e-6)
				for (int i = 0; i < 3; ++i) e1[i] /= len;
			e2 = {
				normal[1]*e1[2] - normal[2]*e1[1],
				normal[2]*e1[0] - normal[0]*e1[2],
				normal[0]*e1[1] - normal[1]*e1[0]
			};
		}

		for (unsigned int f_id : face_ids) {
			const auto& centroid = dualVertices[f_id];
			std::vector<double> vec = {
				centroid[0] - origin[0],
				centroid[1] - origin[1],
				centroid[2] - origin[2]
			};
			double x = vec[0]*e1[0] + vec[1]*e1[1] + vec[2]*e1[2];
			double y = vec[0]*e2[0] + vec[1]*e2[1] + vec[2]*e2[2];
			double angle = atan2(y, x);
			angle_face_pairs.emplace_back(angle, f_id);
		}

		sort(angle_face_pairs.begin(), angle_face_pairs.end());

		vector<unsigned int> ordered;
		for (const auto& [angle, fid] : angle_face_pairs) {
			ordered.push_back(fid);
		}

		NewCell2DsId.push_back(v_idx);
		NewCell2DsVertices.push_back(ordered);
	}


	for (const auto& entry : mappaVertexToDual) { // Scorro la mappa dei vertici originali e i loro vertici duali più vicini
		unsigned int originalVertexId = entry.first; // ID del vertice originale
		const vector<unsigned int>& dualVertexIds = entry.second; // Indici dei vertici duali più vicini
		NewCell2DsId.push_back(originalVertexId); // Aggiungo l'ID del vertice originale alla lista degli ID delle facce duali
		NewCell2DsVertices.push_back(dualVertexIds); // Aggiungo gli ID dei vertici duali come un nuovo vettore di vertici per la faccia duale
	}

	vector<vector<unsigned int>> NewCell2DsEdges(mappaVertexToDual.size()); // Inizializza correttamente

	for (unsigned int i = 0; i < edgeToFaces.size(); ++i) //ogni i è un edge duale
	{
		const auto& edge = edgeToFaces[i]; // Prendo la coppia di facce che condividono lo spigolo
		
		unsigned int j = 0;
		for (const auto& entry : mappaVertexToDual) {
    		const auto& dualVerts = entry.second;
    		if (dualVerts.size() >= 2) {
        		if (find(dualVerts.begin(), dualVerts.end(), edge[0]) != dualVerts.end() &&
            		find(dualVerts.begin(), dualVerts.end(), edge[1]) != dualVerts.end()) {
            		NewCell2DsEdges[j].push_back(i); // i deve essere l'indice reale dello spigolo duale
        		}
   			 }
    		++j;
		}
	}

	//RIEMPIO CELL0
	mesh.NumCell0Ds = NumFaces; // Imposto il numero di celle 0D (vertici duali) uguale al numero di facce
	mesh.Cell0DsId = NewCell0DsId; // Assegno gli ID dei vertici duali al vettore degli ID dei vertici del poliedro duale
	mesh.Cell0DsCoordinates = NewCell0DsCoordinates; // Assegno le coordinate dei vertici duali alla matrice delle coordinate dei vertici del poliedro duale

	//RIEMPIO CELL1
	mesh.NumCell1Ds = edgeToFaces.size(); // Imposto il numero di celle 1D (spigoli duali) uguale al numero di spigoli duali
	mesh.Cell1DsId = NewCell1DsId; // Assegno gli ID degli spigoli duali al vettore degli ID degli spigoli del poliedro duale
	mesh.Cell1DsExtrema = NewCell1DsExtrema; // Assegno gli estremi degli spigoli duali alla matrice degli estremi degli spigoli del poliedro duale

	//RIEMPIO CELL2
	mesh.NumCell2Ds = mappaVertexToDual.size(); // Imposto il numero di celle 2D (facce duali) uguale al numero di vertici originali
	mesh.Cell2DsId = NewCell2DsId;
	mesh.Cell2DsVertices = NewCell2DsVertices;
	mesh.Cell2DsEdges = NewCell2DsEdges; // Assegno gli spigoli duali alla matrice degli spigoli del poliedro duale
	mesh.Cell2DsNumVertices.assign(mappaVertexToDual.size(), NewCell2DsVertices[0].size()); // Assegno il numero di vertici di ogni faccia duale (tutti hanno lo stesso numero di vertici)
	mesh.Cell2DsNumEdges.assign(mappaVertexToDual.size(), NewCell2DsVertices[0].size()); // Assegno il numero di spigoli di ogni faccia duale (tutti hanno lo stesso numero di spigoli)

	//RIEMIPIO CELL3
	mesh.NumCell3Ds = 1; // Il poliedro duale ha una sola cella 3D (il poliedro stesso)
	mesh.Cell3DsId = 0; // L'ID della cella 3D del poliedro duale è 0
	
	//SENZA TRIANGOLAZIONE VALGONO LE RIGHE SEGUENTI
	mesh.Cell3DsNumVertices = NumFaces; // Il numero di vertici del poliedro duale è uguale al numero di facce del poliedro originale
	mesh.Cell3DsNumEdges = edgeToFaces.size(); // Il numero di spigoli del poliedro duale è uguale al numero di spigoli duali
	mesh.Cell3DsNumFaces = mappaVertexToDual.size(); // Il numero di facce del poliedro duale è uguale al numero di vertici originali
	mesh.Cell3DsVertices = NewCell0DsId; // Assegno gli ID dei vertici del poliedro duale
	mesh.Cell3DsEdges = NewCell1DsId; // Assegno gli ID degli spigoli del poliedro duale
	mesh.Cell3DsFaces = NewCell2DsId; // Assegno gli ID delle facce del poliedro duale
    
}
}   
