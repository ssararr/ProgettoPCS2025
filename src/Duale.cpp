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

//Definisco una struct per la definizione di punti 3d in cui è presente il costruttore e la funzione che calcola la distanza euclidea tra due punti 3D
/*struct Point3D {
    double x, y, z;

    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {} // Costruttore per inizializzare le coordinate del punto

    // Calcola la distanza euclidea tra questo punto e un altro
    double distanceTo(const Point3D& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
};*/

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
}

/*if (matrix.rows() != 3) {
        throw invalid_argument("La matrice deve avere 3 righe (x, y, z).");
    }

    vector<vector<double>> result;
    result.reserve(matrix.cols());  // Numero di punti (colonne)

    for (int i = 0; i < matrix.cols(); ++i) {
        // Crea un vettore 3D direttamente dalle righe x, y, z della colonna i-esima
        vector<double> point = {
            matrix(0, i),  // x
            matrix(1, i),  // y
            matrix(2, i)   // z
        };
        result.push_back(point);
    }

    return result;
}*/

//Definisco la funzione findClosestDualVertices che prende in input i due vettori di punti3D che contengono
// le coordinate dei vertici originali e dei vertici duali, e restituisce in output una mappa 
// che associa ad ogni vertice originale gli indici dei vertici duali più vicini, entro una certa tolleranza.

/*const vector<vector<double>>& originalVertices,
    const vector<vector<double>>& dualVertices,
    double tolerance 
) {
    map<int, vector<unsigned int>> vertexToDualMap; //definisco la mappa che vorrò in output

    for (unsigned int i = 0; i < originalVertices.size(); ++i) {
        const vector<double> original = originalVertices[i]; // Prendo il vertice originale corrente (oggetto Point3D)
        double minDistance = numeric_limits<double>::max(); //inizializzo la minima distanza al valore massimo possibile, in modo tale che qualsiasi prima distanza misurata sia minore di quella inizializzata
        vector<unsigned int> closestDualIndices; // Vettore per memorizzare gli indici dei vertici duali più vicini a questo vertice originale

	//CONTROLLO DELLE DIMENSIONI DEI VERTICI PRIMA DEL RICHIAMO DI DISTANCE BETWEEN
	for (unsigned int i = 0; i < originalVertices.size(); ++i) {
		cout << "Vertice originale " << i << ": ";
		for (double coord : originalVertices[i]) cout << coord << " ";
		cout << endl;
	}

	for (unsigned int j = 0; j < dualVertices.size(); ++j) {
		cout << "Vertice duale " << j << ": ";
		for (double coord : dualVertices[j]) cout << coord << " ";
			cout << endl;
	}


        // Trova la distanza minima per questo vertice originale
        for (unsigned int j = 0; j < dualVertices.size(); ++j) {
            double dist = distanceBtw(original, dualVertices[j]); //uso la funzione distanceTo della struct Point3D per calcolare la distanza tra il vertice originale e il vertice duale corrente
            
            if (dist < minDistance - tolerance) {
                // Nuovo minimo trovato, resetto la lista e inserisco solo il nuovo vertice corrente
                minDistance = dist;
                closestDualIndices.clear();
                closestDualIndices.push_back(j);
            }
            else if (abs(dist - minDistance) < tolerance) {
                // Distanza equivalente al minimo, aggiungo alla lista tale vertice duale (poiché lavoriamo con poligoni regolari, il ruolo della tolleranza sarebbe evitabile, visto che tutti i vertici che cerchiamo avranno la stessa distanza minima da quello origianale. Valutare se eliminarla)
                closestDualIndices.push_back(j);
            }
        }

        vertexToDualMap[i] = closestDualIndices; // Associo l'indice del vertice originale agli indici dei vertici duali più vicini
    }

    return vertexToDualMap; //restituisco la mappa 
}*/

map<int, vector<unsigned int>> findClosestDualVertices(
    const vector<vector<double>>& originalVertices,
    const vector<vector<double>>& dualVertices) 
{
	double tolerance = 1e-6; // Definisco una tolleranza per considerare due distanze equivalenti (può essere modificata in base alle esigenze)
    map<int, vector<unsigned int>> vertexToDualMap;

    // Verifica preliminare delle dimensioni
    for (const auto& v : originalVertices) {
        if (v.size() != 3) {
            throw invalid_argument("Tutti i vertici originali devono essere 3D");
        }
    }
    for (const auto& v : dualVertices) {
        if (v.size() != 3) {
            throw invalid_argument("Tutti i vertici duali devono essere 3D");
        }
    }

    for (unsigned int i = 0; i < originalVertices.size(); ++i) {
        const auto& original = originalVertices[i]; // [x,y,z]
        double minDistance = numeric_limits<double>::max();
        vector<unsigned int> closestDualIndices;

        for (unsigned int j = 0; j < dualVertices.size(); ++j) {
            const auto& dual = dualVertices[j]; // [x,y,z]
            
            // Calcolo distanza euclidea manuale (ottimizzato per 3D)
            double dx = original[0] - dual[0];
            double dy = original[1] - dual[1];
            double dz = original[2] - dual[2];
            double dist = sqrt(dx*dx + dy*dy + dz*dz);

            if (dist < minDistance - tolerance) {
                minDistance = dist;
                closestDualIndices = {j}; // Sostituisci con nuovo minimo
            } 
            else if (abs(dist - minDistance) <= tolerance) {
                closestDualIndices.push_back(j); // Aggiungi a minimi esistenti
            }
        }

        vertexToDualMap[i] = closestDualIndices;
    }

    return vertexToDualMap;
}





namespace PolyhedraLibrary {

bool MeshDuale(PolyhedraMesh& mesh) 
{
	//ASSOCIO OGNI FACCIA AL SUO VERTICE DUALE (RIEMPIO NewCell0DsId e NewCell0DsCoordinates)
	unsigned int NumFaces = mesh.Cell3DsNumFaces;
	if (NumFaces == 0) {
		cerr << "Non sono presenti facce nella mesh." << endl;
		return false;
	}

	//(inutile?)map<unsigned int, vector<double>> dualVertexToFace; //Mappa che associa ad ogni id di ogni faccia del poliedro iniziale le coordinate del suo baricentro (ovvero il nuovo vertice duale)
	vector<unsigned int> NewCell0DsId(NumFaces); 	
	MatrixXd NewCell0DsCoordinates(3, NumFaces); // Inizializzo la matrice delle coordinate dei vertici duali

	for (unsigned int faceID = 0; faceID < NumFaces; faceID++) {
		//gestisci errore(?)
		if (mesh.Cell2DsNumVertices[faceID] == 0) {
			cerr << "La faccia faceID " << faceID << " non ha vertici." << endl;
			return false;
		}
		Vector3d baricentro = Vector3d::Zero(); // Inizializzo il baricentro a zero
		for (const auto& id_vertex : mesh.Cell2DsVertices[faceID]) {
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
			dualEdges.emplace_back(facesSharingEdge[0], facesSharingEdge[1]); // Aggiungo la coppia di ID delle facce come un nuovo spigolo duale
		}
	}
	
	MatrixXi NewCell1DsExtrema(2, dualEdges.size()); // Inizializzo la matrice degli estremi degli spigoli duali
	vector<unsigned int> NewCell1DsId(dualEdges.size()); // Vettore per gli ID degli spigoli duali
	for (unsigned int i = 0; i < edgeToFaces.size(); ++i)
	{
		const auto& edge = edgeToFaces[i]; // Prendo la coppia di facce che condividono lo spigolo
		NewCell1DsExtrema(0, i) = edge[0]; // Assegno il primo ID della coppia come primo estremo dello spigolo duale
		NewCell1DsExtrema(1, i) = edge[1]; // Assegno il secondo ID della coppia come secondo estremo dello spigolo duale
		NewCell1DsId.push_back(i); // Aggiungo l'ID dello spigolo duale
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



	mappaVertexToDual = findClosestDualVertices(originalVertices, dualVertices); // Trovo i vertici duali più vicini per ogni vertice originale /*

	vector<unsigned int> NewCell2DsId;
	NewCell2DsId.reserve(mappaVertexToDual.size());
	vector<vector<unsigned int>> NewCell2DsVertices;
	NewCell2DsVertices.reserve(mappaVertexToDual.size()); // Riservo spazio per gli ID dei vertici duali più vicini per ogni vertice originale
	for (const auto& entry : mappaVertexToDual) { // Scorro la mappa dei vertici originali e i loro vertici duali più vicini
		unsigned int originalVertexId = entry.first; // ID del vertice originale
		const vector<unsigned int>& dualVertexIds = entry.second; // Indici dei vertici duali più vicini
		NewCell2DsId.push_back(originalVertexId); // Aggiungo l'ID del vertice originale alla lista degli ID delle facce duali
		NewCell2DsVertices.push_back(dualVertexIds); // Aggiungo gli ID dei vertici duali come un nuovo vettore di vertici per la faccia duale
	}

	vector<vector<unsigned int>> NewCell2DsEdges; 
	NewCell2DsEdges.reserve(mappaVertexToDual.size()); // Riservo spazio per gli ID degli spigoli duali per ogni faccia duale
	for (unsigned int i = 0; i < edgeToFaces.size(); ++i) //ogni i è un edge duale
	{
		const auto& edge = edgeToFaces[i]; // Prendo la coppia di facce che condividono lo spigolo
		
		for (unsigned int j = 0; j < mappaVertexToDual.size(); j++) //ogni j è una faccia duale
		{
			if (mappaVertexToDual[j].size() >= 2) { // Se ci sono almeno due vertici duali nella faccia duale
				if (find(mappaVertexToDual[j].begin(), mappaVertexToDual[j].end(), edge[0]) != mappaVertexToDual[j].end() &&
					find(mappaVertexToDual[j].begin(), mappaVertexToDual[j].end(), edge[1]) != mappaVertexToDual[j].end()) {			
					NewCell2DsEdges[j].push_back(i); // Aggiungo lo spigolo duale i nella lista di spigoli della faccia duale j
				}
			}
		}
	}

	
	//RIEMPIO CELL0
	mesh.NumCell0Ds = NumFaces; // Imposto il numero di celle 0D (vertici duali) uguale al numero di facce
	//mesh.Cell0DsId.resize(NumFaces); // Riservo spazio per gli ID dei vertici duali
	mesh.Cell0DsId = NewCell0DsId; // Assegno gli ID dei vertici duali al vettore degli ID dei vertici del poliedro duale
	//mesh.Cell0DsCoordinates.resize(3, NumFaces); // Inizializzo la matrice delle coordinate dei vertici duali
	mesh.Cell0DsCoordinates = NewCell0DsCoordinates; // Assegno le coordinate dei vertici duali alla matrice delle coordinate dei vertici del poliedro duale

	//RIEMPIO CELL1
	mesh.NumCell1Ds = edgeToFaces.size(); // Imposto il numero di celle 1D (spigoli duali) uguale al numero di spigoli duali
	//mesh.Cell1DsId.resize(edgeToFaces.size()); // Riservo spazio per gli ID degli spigoli duali
	mesh.Cell1DsId = NewCell1DsId; // Assegno gli ID degli spigoli duali al vettore degli ID degli spigoli del poliedro duale
	//mesh.Cell1DsExtrema.resize(2, edgeToFaces.size()); // Inizializzo la matrice degli estremi degli spigoli duali
	mesh.Cell1DsExtrema = NewCell1DsExtrema; // Assegno gli estremi degli spigoli duali alla matrice degli estremi degli spigoli del poliedro duale

	//RIEMPIO CELL2
	mesh.NumCell2Ds = mappaVertexToDual.size(); // Imposto il numero di celle 2D (facce duali) uguale al numero di vertici originali
	//mesh.Cell2DsId.resize(mappaVertexToDual.size());
	mesh.Cell2DsId = NewCell2DsId;
	//mesh.Cell2DsVertices.resize(mappaVertexToDual.size());
	mesh.Cell2DsVertices = NewCell2DsVertices;
	//mesh.Cell2DsEdges.resize(mappaVertexToDual.size());
	mesh.Cell2DsEdges = NewCell2DsEdges; // Assegno gli spigoli duali alla matrice degli spigoli del poliedro duale
	//mesh.Cell2DsNumVertices.resize(mappaVertexToDual.size()); // Riservo spazio per il numero di vertici di ogni faccia duale
	mesh.Cell2DsNumVertices.assign(mappaVertexToDual.size(), NewCell2DsVertices[0].size()); // Assegno il numero di vertici di ogni faccia duale (tutti hanno lo stesso numero di vertici)
	//mesh.Cell2DsNumEdges.resize(mappaVertexToDual.size()); // Riservo spazio per il numero di spigoli di ogni faccia duale
	mesh.Cell2DsNumEdges.assign(mappaVertexToDual.size(), NewCell2DsVertices[0].size()); // Assegno il numero di spigoli di ogni faccia duale (tutti hanno lo stesso numero di spigoli)

	//RIEMIPIO CELL3
	mesh.NumCell3Ds = 1; // Il poliedro duale ha una sola cella 3D (il poliedro stesso)
	mesh.Cell3DsId = 0; // L'ID della cella 3D del poliedro duale è 0
	
	//SENZA TRIANGOLAZIONE VALGONO LE RIGHE SEGUENTI
	mesh.Cell3DsNumVertices = NumFaces; // Il numero di vertici del poliedro duale è uguale al numero di facce del poliedro originale
	mesh.Cell3DsNumEdges = edgeToFaces.size(); // Il numero di spigoli del poliedro duale è uguale al numero di spigoli duali
	mesh.Cell3DsNumFaces = mappaVertexToDual.size(); // Il numero di facce del poliedro duale è uguale al numero di vertici originali
	//mesh.Cell3DsVertices.resize(NumFaces); // Riservo spazio per gli ID dei vertici del poliedro duale
	mesh.Cell3DsVertices = NewCell0DsId; // Assegno gli ID dei vertici del poliedro duale
	//mesh.Cell3DsEdges.resize(edgeToFaces.size()); // Riservo spazio per gli ID degli spigoli del poliedro duale
	mesh.Cell3DsEdges = NewCell1DsId; // Assegno gli ID degli spigoli del poliedro duale
	//mesh.Cell3DsFaces.resize(mappaVertexToDual.size()); // Riservo spazio per gli ID delle facce del poliedro duale
	mesh.Cell3DsFaces = NewCell2DsId; // Assegno gli ID delle facce del poliedro duale
    
}
}   
