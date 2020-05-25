/**	\file	PCPMegaTissueClass.h
*	Header of PCPMegaTissue class.
*	Contains data and function declarations associated with PCPMegaTissue objects.
*	Depends on base class Tissue and on subclass PCPTissue.
*/

// include guard
#ifndef __PCPMEGATISSUECLASS_H_INCLUDED__
#define __PCPMEGATISSUECLASS_H_INCLUDED__

//Forward declared dependences
class PCPCell;
class PCPEdge;
class PCPTissue;
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
#include "VectClass.h"
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h"
#include "PCPTissueClass.h"
#include "PCPEdgeClass.h"
#include "PCPConstants.h"
#include "PCPMegaCellClass.h"
#include <random>
#include <chrono>								//temporary, for measuring time and debugging evolve().
//Changes by Hayden for PCP: Now I am going to use a dirichlet random number generation from gsl for
//initialize PCP concentrations on the cell edges.
#include "gsl/gsl_rng.h"
#include "gsl/gsl_matrix.h"						//GSL MATRIX FUNCTIONALITY
#include "gsl/gsl_blas.h"						//GSL BASIC LINEAR ALGEBRA WITH MATRICES
#include "gsl/gsl_vector.h"						//GSL VECTOR FUNCTIONALITY


using namespace std;

//! Class for PCP tissue object.
class PCPMegaTissue : public PCPTissue{
public:
	//ADDITIONALDATA

	//FUNCTIONS
	//CONSTRUCTOR
	PCPMegaTissue(){};							//!< Default constructor. 
	virtual ~PCPMegaTissue(){};

	//! Initialize PCP concentrations within a line of cells in the tissue around the central MEGACELL.
	virtual void AddDifferentiatedCellColumn();

	//! Calculates the force on every vertex and stores the result in vertex: no adjustment of box size.
	virtual void ForceFunction(ofstream* errorfile);

	//! Evolve that includes special functions for mega cell simulation (for y junction formation)
	virtual void Evolve(double timeInterval, double* prot_stepSize, double* mech_stepSize, ofstream* errorfile);
	
	//HELPER FUNCTIONS

	//! Check that there is only one inner cell and one outer cell, and that all other cell types are the standard cell type.
	virtual void HelperReadInputData_CheckTypes();

	//! Reads the data and sets all of the global variables.
	//! Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
	virtual void HelperReadInputData_SetGlobalValues(ifstream* datafile);

	//! Writes data into file specified by 'datafile'
	virtual void HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo);

	//! Locate trapped_defect: in next_differentiated_cells, identify trapped cells
	virtual vector<PCPCell*> locate_trapped_defect(vector<PCPCell*> next_differentiated_cells);

	//! Locate trapped_defect: in next_differentiated_cells, identify trapped cells
	virtual vector<PCPCell*> locate_cells_for_deactivation(vector<PCPCell*> next_differentiated_cells);

	//! Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids.
	//! Starts at current position of *datafile and reads until keyword 'End'.
	virtual void HelperReadInputData_CreateObjects(ifstream* datafile);

	//! This function is used to loop through the list of undifferentiated cells that are in contact with differentiated cells at a propagation time.
	// \note Add a partial column after a couple of full columns.
	virtual void HelperEvolve_GeneratePCPProtein(bool delay, bool added_partial_column);

	//! This is a helper function for generating a checkerboard pattern in the next row: for a circular ring with a vertical stripe rather than for periodic boundary conditions.
	void HelperCheckerboard_GeneratePattern(vector<PCPCell*> next_differentiated_cells, ofstream* datafile);

	//! Order cell pointer by the angle that the cell centroid about the central mega cell.
	vector<PCPCell*> order_cells_angle(const vector<PCPCell*> &unordered_cells);

	//! Output vector of angles about the central mega cell.
	// \note The vector obeys the conventions of order_cell_angles.
	vector<double> generate_cell_angle_vector(const vector<PCPCell*> &unordered_cells);

	//! Find angle of cell about central megacell.
	double Helper_Find_Cell_Angle(PCPCell* pcp_cell_ptr);

	//! Make sure that you propogate PCP outward from origin_mega_cell.
	Edge* edge_candidate_outward(const vector<Edge*> &edge_candidates);

	//! Check fourfold vertex for stability (make sure that this function creates appropriate types of objects -> PCPEdge and PCPMegaCell)
	virtual double HelperResolveVertices_CheckStability(Vertex* v, Edge* e1, Edge* e2, Edge* e3, Edge* e4, Cell* c1, Cell* c2, Cell* c3, Cell* c4, Vect *Force);

};

#endif