#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Eigen"
#include <string>
#include <list>
#include <map>
#include <vector>
#include <set>

namespace PolygonalLibrary
{
bool ImportMesh(PolygonalMesh& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}
// ***************************************************************************
bool ImportCell0Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell0Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    //rimuovo la prima riga
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2d coord;
		char elimin; // ';' del csv

        converter >> id >> elimin >> marker >> elimin >> mesh.Cell0DsCoordinates(0, id) >> elimin >> mesh.Cell0DsCoordinates(1, id);

        mesh.Cell0DsId.push_back(id);

        /// Memorizza i marker
        map<unsigned int, list<unsigned int>>& m = mesh.cell0Ds_markers;
        if (marker != 0)
		{
			auto [itor, bool_val] = m.try_emplace(marker);
			itor -> second.push_back(id);
		}

    }

    return true;
}

// ***************************************************************************
bool ImportCell1Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
	
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    //rimuovo la prima riga
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2i vertices;
		char elimin; //';' del csv

        converter >> id >> elimin >> marker >> elimin >> mesh.Cell1DsExtrema(0, id) >> elimin >> mesh.Cell1DsExtrema(1, id);

        mesh.Cell1DsId.push_back(id);

        /// Memorizza i marker
        map<unsigned int, list<unsigned int>>& m = mesh.cell1Ds_markers;
        if (marker != 0) 
		{
			auto [itor, bool_val] = m.try_emplace(marker);
			itor -> second.push_back(id);
		}
    }

    return true;
}

// ***************************************************************************
bool ImportCell2Ds(PolygonalMesh& mesh)
{
    ifstream file;
    file.open("./Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    //rimuovo la prima riga
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
		unsigned int marker;
		unsigned int nVert;
		unsigned int nEdg;
		char elimin;
		
		converter >> id >> elimin >> marker >> elimin >> nVert;
		
        vector<unsigned int> v_vertices;
		v_vertices.reserve(nVert);
		for(unsigned int i = 0; i < nVert; i++)
		{
			unsigned int vertices;
			converter >> elimin >> vertices;
			v_vertices.push_back(vertices);
		}
		mesh.Cell2DsVertices.push_back(v_vertices);
		
		converter >> elimin >> nEdg;
		
        vector<unsigned int> v_edges;
		v_edges.reserve(nEdg);
		for(unsigned int j = 0; j < nEdg; j++)
		{
			unsigned int edge;
			converter >> elimin >> edge;
			v_edges.push_back(edge);
		}
		
		mesh.Cell2DsEdges.push_back(v_edges);
		
		mesh.Cell2DsId.push_back(id);

		/// Memorizza i marker
        map<unsigned int, list<unsigned int>>& m = mesh.cell2Ds_markers;
        if (marker != 0) 
		{
			auto [itor, bool_val] = m.try_emplace(marker);
			itor -> second.push_back(id);
		}    
    
    }

    return true;
}

}
