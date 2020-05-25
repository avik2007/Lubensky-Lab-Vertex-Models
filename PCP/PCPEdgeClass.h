/**	\file	PCPEdgeClass.h
*	Header of PCPEdge class.
*	Contains data and function declarations associated with PCPEdge objects.
*	Depends on base class Edge.
*/

// include guard
#ifndef __PCPEDGECLASS_H_INCLUDED__
#define __PCPEDGECLASS_H_INCLUDED__

//Forward declared dependences
class PCPCell;
class PCPEdge;
class Cell;
class Edge;
class Vertex;
class Vect;

//General files
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>

//Base class files
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h"

//PCP subclass files
#include "PCPCellClass.h"
#include "PCPConstants.h"

using namespace std;


//! Class for cell objects.
class PCPEdge : public Edge{
public:
	//ADDITIONAL DATA
	vector<double> planarpolarity;

	//ADDITIONAL FUNCTIONS
	//CONSTRUCTORS AND DESTRUCTORS
	PCPEdge(){};
	virtual ~PCPEdge(){};
	//! Constructs an edge with id 'idnumb', vertices 'v1' and 'v2', and flags 'xflag', 'yflag'.
	PCPEdge(unsigned long int idnumb, int v1, int v2, bool flagx, bool flagy){
		id = idnumb;
		yflag = flagy;
		xflag = flagx;
		vertex1 = NULL;
		vertex2 = NULL;
		cell1 = NULL;
		cell2 = NULL;
		IOVertexID[0] = v1;
		IOVertexID[1] = v2;
		//Changes by Hayden for PCP: initialize polarity protein amounts to zero.
		planarpolarity.resize(4);
		fill(planarpolarity.begin(), planarpolarity.end(), 0.0);
	}

	PCPEdge(unsigned long int idnum, bool flagx, bool flagy, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2){
		id = idnum;
		yflag = flagy;
		xflag = flagx;
		vertex1 = v1;
		vertex2 = v2;
		cell1 = c1;
		cell2 = c2;
		//Changes by Hayden for PCP: initialize polarity protein amounts to zero.
		planarpolarity.resize(4);
		fill(planarpolarity.begin(), planarpolarity.end(), 0.0);
	}

	//Checks whether a edge has PCP.
	bool is_Edge_with_PCP();

	//Assign type to an edge based on the cells that share it (allow increased tension between partnered cells).
	int Define_Edge_Type();

	//Returns the magnitude of the tension on the edge (this is implemented differently in the edge subclasses) 
	virtual double GetTension_Magnitude();

	//Returns the magnitude of tension on a vertex due to a virtual edge given the concentration for a virtual edge specified in Lubensky PLOS 2012.
	double GetTension_Magnitude_Virtual(Vertex* v, PCPEdge* edge1_cellalpha, PCPEdge* edge2_cellalpha, PCPEdge* edge1_cellbeta, PCPEdge* edge2_cellbeta);

	//Returns the tension on a vertex due to a virtual edge given the concentration for a virtual edge specified in Lubensky PLOS 2012.
	Vect GetTension_ForceOnVertex_Virtual(Vertex* v, PCPEdge* edge1_cellalpha, PCPEdge* edge2_cellalpha, PCPEdge* edge1_cellbeta, PCPEdge* edge2_cellbeta);
};

#endif