/**	\file	PCPCellClass.h
*	Header of PCPCell class.
*	Contains data and function declarations associated with PCPCell objects.
*	Depends on base class Cell.
*/

// include guard
#ifndef __PCPCELLCLASS_H_INCLUDED__
#define __PCPCELLCLASS_H_INCLUDED__

//Forward declared dependences
class Cell;
class Edge;
class Vertex;
class Vect;
class PCPCell;
class PCPEdge;

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h" 
#include "PCPEdgeClass.h"

//Changes by Hayden for PCP: Now I am going to use a dirichlet random number generation from gsl for
//initialize PCP concentrations on the cell edges.
#include "gsl\gsl_rng.h"			//Used of random generation of initial PCP protein concentrations
#include "gsl\gsl_randist.h"		//Used of random generation of initial PCP protein concentrations
//GSL files
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"			//GSL MATRIX LINEAR ALGEBRA
#include "gsl/gsl_randist.h"

using namespace std;


//! Class for cell objects.
class PCPCell : public Cell{
public:
	vector<vector<double>> V_matrix;				//Store V matrix (dependence on edge lengths) for the cell
	double time_of_pcp;								//This stores the time at which a cell got PCP (just a time stamp: will not be edited).
	int cell_generation;							//Store the generation of cell defining the layer of PCP cells.

	//ADDITIONAL FUNCTIONS
	//CONSTRUCTORS AND DESTRUCTORS
	PCPCell(){};
	virtual ~PCPCell(){};

	//! Checks whether a cell has PCP.
	bool is_Cell_with_PCP();

	//! This function is used to calculate the distance from the origin of PCP.
	virtual double calculate_Distance_from_PCP_Origin();

	//! This function checks if the cell is in contact with a cell with PCP.
	bool CheckNeighboringCell_forPCP();

	//! Find minimum cell generation neighboring the cell, update cell generation to min generation + 1.
	void Update_Cell_Generation();

	//! Returns a list of ordered edges (given a list of ordered vertices)
	vector<PCPEdge*> FindOrderedEdges();

	//! Return next edge CW from (*eptr)
	PCPEdge* GetNextCWEdge(PCPEdge* eptr);

	//! NOTE: This is only useful for simulations on a hexagonal field of cells.
	PCPEdge* FindUpperHorizontalEdge();

	//! Order edges from list of ordered vertices. NOTE: This is only used for simulations on a hexagonal field of cells.
	vector<PCPEdge*> OrderEdges();

	//HELPER FUNCTIONS

	//! Generate a random time for the given cell.
	void Helper_GenerateRandomTime();

	//! When a cells gets PCP, start clocks in neighboring cells.
	virtual void Helper_AddBorderingCellsPoisson();

	//! This function is used to create an initially randomized distribution of PCP proteins in a cell. Called by Evolve.
	void HelperEvolve_InitializePolarityProteins();

	//! This function is used to create a stereotyped distribution of protein for either edge polarized or vertex polarized states. Called by ReadInputData.
	void HelperReadInputData_InitializePolarityProteins();

	//! This function is used to randomized the protein distribution within a cell about the specified protein distribution (from input file). Called by ReadInputData.
	void HelperReadInputData_RandomizePolarityProteins();

	//! This function sets all protein concentrations in a cell to be identical to those of cell pointed to by cell_ptr (only applicable in case of hexagonal field of cells).
	void HelperReadInputData_SetProteinConc_fromCellPointer(PCPCell* cell_ptr);
};

#endif