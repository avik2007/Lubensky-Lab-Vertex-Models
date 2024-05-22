#pragma once
// guard
#ifndef __AGSTIFFTISSUE_H_INCLUDED__
#define __AGSTIFFTISSUE_H_INCLUDED__

// Forward declared dependencies
class Tissue;
class Cell;
class Edge;
class Vertex;
class Vect;

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include <stddef.h>
#include <complex> 
#include "Eigen/Dense"

#include "VectClass.h"
#include "Constants.h"
#include "AGConstants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "AGCell.h"
#include "VertexClass.h"
#include "TissueClass.h"
#include "AGTissue.h"
#include "AGPerimConstants.h"

using namespace Eigen;
using namespace std;

class AGStiffTissue : public AGPerimTissue {
	/*
	* It's important to mention, this class descends from AGTissue but is not really intended to grow.
	*
	*/
public:
	//CONSTRUCTOR
	AGStiffTissue() {};
	
	//DESTRUCTOR
	virtual ~AGStiffTissue() {};

	//DATA
	bool constantK;
	double meanK;
	double Kstd;
	//ADDITIONAL FUNCTIONS
	// Updates read input function to accept new global variables 
	virtual void HelperReadInputData_SetGlobalValues(ifstream* datafile);

	//Updates read input function to create new objects (AGCells)
	virtual void HelperReadInputData_CreateObjects(ifstream* datafile);

	//! Writes data into file specified by 'datafile'
	virtual void HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo);

	//calculates the free energy of a packing 
	virtual double CalculateFreeEnergy();

	//calculates the area energy of a packing
	virtual double CalculatePressureEnergy();

	virtual bool ReadInputData(string InputFileName);

	//sets the spatially varying area modulus (KAM) for each cells
	virtual void AssignAreaModuli(bool constantK, double meanK, double Kstd);

	//! Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration
	//! Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors
	//! returns the magnutude of Force-Forcenewedge and changes the value of Force
	virtual double HelperResolveVertices_CheckStability(Vertex* v, Edge* e1, Edge* e2, Edge* e3, Edge* e4, Cell* c1, Cell* c2, Cell* c3, Cell* c4, Vect* Force);

	//! Calculates the force on every vertex and stores the result in vertex.
	virtual void ForceFunction(ofstream* errorfile);

	//! Sets the preferred area of A0 to a constant value for all the cells
	virtual void AssignPreferredAreas(double A0);
};

#endif