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
//(la libreria Eigen dev'essere già inclusa in Duale.hpp e PolyhedraMesh.hpp)

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

vector<vector<double>> eigenMatrixToVectorVector(const MatrixXd& matrix) {
    vector<vector<double>> result;
    result.reserve(matrix.rows());  // Ottimizzazione: prealloca le righe
    
    for (int i = 0; i < matrix.rows(); ++i) {
        // Copia ogni riga della MatrixXd in un vector<double>
        vector<double> row(matrix.data() + i * matrix.cols(), 
                               matrix.data() + (i + 1) * matrix.cols());
        result.push_back(row);
    }
    
    return result;
}

//Definisco la funzione findClosestDualVertices che prende in input i due vettori di punti3D che contengono
// le coordinate dei vertici originali e dei vertici duali, e restituisce in output una mappa 
// che associa ad ogni vertice originale gli indici dei vertici duali più vicini, entro una certa tolleranza.

map<int, vector<int>> findClosestDualVertices(
    const vector<vector<double>>& originalVertices,
    const vector<vector<double>>& dualVertices,
    double tolerance = 1e-6
) {
    map<int, vector<int>> vertexToDualMap; //definisco la mappa che vorrò in output

    for (int i = 0; i < originalVertices.size(); ++i) {
        const vector<double> original = originalVertices[i]; // Prendo il vertice originale corrente (oggetto Point3D)
        double minDistance = numeric_limits<double>::max(); //inizializzo la minima distanza al valore massimo possibile, in modo tale che qualsiasi prima distanza misurata sia minore di quella inizializzata
        vector<int> closestDualIndices; // Vettore per memorizzare gli indici dei vertici duali più vicini a questo vertice originale

        // Trova la distanza minima per questo vertice originale
        for (int j = 0; j < dualVertices.size(); ++j) {
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
}



namespace PolyhedraLibrary {

bool MeshDuale(PolyhedraMesh& mesh) 
{
	//ASSOCIO OGNI FACCIA AL SUO VERTICE DUALE (RIEMPIO NewCell0DsId e NewCell0DsCoordinates)
	unsigned int NumFaces = Cell3DsNumFaces;
	if (NumFaces == 0) {
		cerr << "Non sono presenti facce nella mesh." << endl;
		return false;
	}

	//(inutile?)map<unsigned int, vector<double>> dualVertexToFace; //Mappa che associa ad ogni id di ogni faccia del poliedro iniziale le coordinate del suo baricentro (ovvero il nuovo vertice duale)
	vector<unsigned int> NewCell0DsId(NumFaces); 	
	MatrixXd NewCell0DsCoordinates(3, NumFaces); // Inizializzo la matrice delle coordinate dei vertici duali

	for (unsigned int faceID = 0; faceID < NumFaces; faceID++) {
		//gestisci errore(?)
		/*if (mesh.Cell2DsNumVertices[faceID].size() == 0) {
			cerr << "La faccia faceID " << faceID << " non ha vertici." << endl;
			return false;
		}*/
		Vector3d baricentro = Vector3d::Zero(); // Inizializzo il baricentro a zero
		for (const auto& id_vertex : mesh.Cell2DsVertices[faceID]) {
			baricentro += mesh.Cell0DsCoordinates.col(id_vertex); // Sommo le coordinate dei vertici della faccia
		}		
		baricentro /= mesh.Cell2DsNumVertices[faceID].size(); // Calcolo il baricentro dividendo per il numero di vertici della faccia
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
		if (faces.size() == 2) { // Se lo spigolo è condiviso da due facce
			dualEdges.emplace_back(facesSharingEdge[0], facesSharingEdge[1]); // Aggiungo la coppia di ID delle facce come un nuovo spigolo duale
		}
	}
	
	MatrixXi NewCell1DsExtrema(2, dualEdges.size()); // Inizializzo la matrice degli estremi degli spigoli duali
	vector<unsigned int> NewCell1DsId(dualEdges.size()); // Vettore per gli ID degli spigoli duali
	for (unsigned int i = 0; i < dualEdges.size(); ++i)
	{
		const auto& edge = dualEdges[i]; // Prendo la coppia di facce che condividono lo spigolo
		NewCell1DsExtrema(0, i) = edge.first; // Assegno il primo ID della coppia come primo estremo dello spigolo duale
		NewCell1DsExtrema(1, i) = edge.second; // Assegno il secondo ID della coppia come secondo estremo dello spigolo duale
		NewCell1DsId.push_back(i); // Aggiungo l'ID dello spigolo duale
	}


	//DEFINISCO OGNI FACCIA DUALE COME L'INSIEME DEI VERTICI DUALI PIU' VICINI AD OGNI VERTICE ORIGINARIO
	vector<vector<double>> originalVertices = eigenMatrixToVectorVector(mesh.Cell0DsCoordinates); // Converto le coordinate dei vertici originali in un vettore di vettori
	vector<vector<double>> dualVertices = eigenMatrixToVectorVector(NewCell0DsCoordinates); // Converto le coordinate dei vertici duali in un vettore di vettori

	mappaVertexToDual = findClosestDualVertices(originalVertices, dualVertices); // Trovo i vertici duali più vicini per ogni vertice originale

	vector<unsigned int> NewCell2DsId;
	NewCell2DsId.reserve(mappaVertexToDual.size());
	vector<vector<unsigned int>> NewCell2DsVertices;
	NewCell2DsVertices.reserve(mappaVertexToDual.size()); // Riservo spazio per gli ID dei vertici duali più vicini per ogni vertice originale
	for (const auto& entry : mappaVertexToDual) { // Scorro la mappa dei vertici originali e i loro vertici duali più vicini
		unsigned int originalVertexId = entry.first; // ID del vertice originale
		const vector<int>& dualVertexIds = entry.second; // Indici dei vertici duali più vicini
		NewCell2DsId.push_back(originalVertexId); // Aggiungo l'ID del vertice originale alla lista degli ID delle facce duali
		NewCell2DsVertices.push_back(dualVertexIds); // Aggiungo gli ID dei vertici duali come un nuovo vettore di vertici per la faccia duale
	}

	vector<vector<double>> NewCell2DsEdges; 
	NewCell2DsEdges.reserve(mappaVertexToDual.size()); // Riservo spazio per gli ID degli spigoli duali per ogni faccia duale
	for (unsigned int i = 0; i < dualEdges.size(); ++i) //ogni i è un edge duale
	{
		const auto& edge = dualEdges[i]; // Prendo la coppia di facce che condividono lo spigolo
		
		for (unsigned int j = 0; j < mappaVertexToDual.size(), j++) //ogni j è una faccia duale
		{
			if (mappaVertexToDual[j].size() >= 2) { // Se ci sono almeno due vertici duali nella faccia duale
				if edge.first in mappaVertexToDual[j] && edge.second in mappaVertexToDual[j]) // Se entrambi gli estremi dello spigolo sono vertici della j-esima faccia duale
				{
					NewCell2DsEdges[j].push_back(i); // Aggiungo lo spigolo duale i nella lista di spigoli della faccia duale j
				}
			}
		}
	}

	
	//RIEMPIO CELL0
	mesh.NumCell0Ds = NumFaces; // Imposto il numero di celle 0D (vertici duali) uguale al numero di facce
	mesh.Cell0DsId.resize(NumFaces); // Riservo spazio per gli ID dei vertici duali
	mesh.Cell0DsId = NewCell0DsId; // Assegno gli ID dei vertici duali al vettore degli ID dei vertici del poliedro duale
	mesh.Cell0DsCoordinates.resize(3, NumFaces); // Inizializzo la matrice delle coordinate dei vertici duali
	mesh.Cell0DsCoordinates = NewCell0DsCoordinates; // Assegno le coordinate dei vertici duali alla matrice delle coordinate dei vertici del poliedro duale

	//RIEMPIO CELL1
	mesh.NumCell1Ds = dualEdges.size(); // Imposto il numero di celle 1D (spigoli duali) uguale al numero di spigoli duali
	mesh.Cell1DsId.resize(dualEdges.size()); // Riservo spazio per gli ID degli spigoli duali
	mesh.Cell1DsId = NewCell1DsId; // Assegno gli ID degli spigoli duali al vettore degli ID degli spigoli del poliedro duale
	mesh.Cell1DsExtrema.resize(2, dualEdges.size()); // Inizializzo la matrice degli estremi degli spigoli duali
	mesh.Cell1DsExtrema = NewCell1DsExtrema; // Assegno gli estremi degli spigoli duali alla matrice degli estremi degli spigoli del poliedro duale

	//RIEMPIO CELL2
	mesh.NumCell2Ds = mappaVertexToDual.size(); // Imposto il numero di celle 2D (facce duali) uguale al numero di vertici originali
	mesh.Cell2DsId.resize(mappaVertexToDual.size());
	mesh.Cell2DsId = NewCell2DsId;
	mesh.Cell2DsVertices.resize(mappaVertexToDual.size());
	mesh.Cell2DsVertices = NewCell2DsVertices;
	mesh.Cell2DsEdges.resize(mappaVertexToDual.size());
	mesh.Cell2DsEdges = NewCell2DsEdges; // Assegno gli spigoli duali alla matrice degli spigoli del poliedro duale
	mesh.Cell2DsNumVertices.resize(mappaVertexToDual.size()); // Riservo spazio per il numero di vertici di ogni faccia duale
	mesh.Cell2DsNumVertices = NewCell2DsVertices[0].size()* Ones(mappaVertexToDual.size(), 1); // Assegno il numero di vertici di ogni faccia duale (tutti hanno lo stesso numero di vertici)
	mesh.Cell2DsNumEdges.resize(mappaVertexToDual.size()); // Riservo spazio per il numero di spigoli di ogni faccia duale
	mesh.Cell2DsNumEdges = NewCell2DsEdges[0].size() * Ones(mappaVertexToDual.size(), 1); // Assegno il numero di spigoli di ogni faccia duale (tutti hanno lo stesso numero di spigoli)

	//RIEMIPIO CELL3
	mesh.NumCell3Ds = 1; // Il poliedro duale ha una sola cella 3D (il poliedro stesso)
	mesh.Cell3DsId = 0; // L'ID della cella 3D del poliedro duale è 0
	
	//SENZA TRIANGOLAZIONE VALGONO LE RIGHE SEGUENTI
	mesh.Cell3DsNumVertices = NumFaces; // Il numero di vertici del poliedro duale è uguale al numero di facce del poliedro originale
	mesh.Cell3DsNumEdges = dualEdges.size(); // Il numero di spigoli del poliedro duale è uguale al numero di spigoli duali
	mesh.Cell3DsNumFaces = mappaVertexToDual.size(); // Il numero di facce del poliedro duale è uguale al numero di vertici originali
	mesh.Cell3DsVertices.resize(NumFaces); // Riservo spazio per gli ID dei vertici del poliedro duale
	mesh.Cell3DsVertices = NewCell0DsId; // Assegno gli ID dei vertici del poliedro duale
	mesh.Cell3DsEdges.resize(dualEdges.size()); // Riservo spazio per gli ID degli spigoli del poliedro duale
	mesh.Cell3DsEdges = NewCell1DsId; // Assegno gli ID degli spigoli del poliedro duale
	mesh.Cell3DsFaces.resize(mappaVertexToDual.size()); // Riservo spazio per gli ID delle facce del poliedro duale
	mesh.Cell3DsFaces = NewCell2DsId; // Assegno gli ID delle facce del poliedro duale


	return true; 
}