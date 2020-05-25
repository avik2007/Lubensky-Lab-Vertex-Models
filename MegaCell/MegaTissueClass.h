/**	\file	MegaTissueClass.h
*	Header of MegaTissue class.
*	Contains data and function declarations associated with MegaTissue objects.
*	Depends on base class Tissue.
*/

// include guard
#ifndef __MEGATISSUECLASS_H_INCLUDED__
#define __MEGATISSUECLASS_H_INCLUDED__

//Forward-declared dependencies
class Tissue;
class Cell;
class Edge;
class Vertex;
class Vect;
class MegaCell;
class MegaVertex;

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include "VectClass.h"
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h"

using namespace std;

//! Class for PCP tissue object.
class MegaTissue : public Tissue{
public:
	//ADDITIONALDATA

	//FUNCTIONS
	//CONSTRUCTOR
	MegaTissue(){};							//!< Default constructor. 
	virtual ~MegaTissue(){};

	//! Creates a new edge and addes it to the list of edges, returns a pointer to that edge.
	virtual Edge* CreateEdge(bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2);

	//! Reads and constructs simulation from file. The file must be formatted like InputTemplate.txt. This function parses the input file then initializes all the objects and pointers.
	virtual void ReadInputData(string InputFileName);

	//! Check that there is only one inner cell and one outer cell, and that all other cell types are the standard cell type.
	virtual void HelperReadInputData_CheckTypes();

	//! Reads the data and sets all of the global variables.
	//! Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
	virtual void HelperReadInputData_SetGlobalValues(ifstream* datafile);

	//! Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids. Starts at current position of *datafile and reads until keyword 'End'.
	virtual void HelperReadInputData_CreateObjects(ifstream* datafile);

	//! Writes data into file specified by 'datafile'
	virtual void HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo);

	//! Calculates the force on every vertex and stores the result in vertex.
	virtual void ForceFunction(ofstream* errorfile);

	//! Make sure that you propogate PCP outward from origin_mega_cell.
	Edge* edge_candidate_outward(const vector<Edge*> &edge_candidates);

	//! Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration
	//! Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors
	//! returns the magnutude of Force-Forcenewedge and changes the value of Force
	virtual double HelperResolveVertices_CheckStability(Vertex* v, Edge* e1, Edge* e2, Edge* e3, Edge* e4, Cell* c1, Cell* c2, Cell* c3, Cell* c4, Vect *Force);
};

#endif