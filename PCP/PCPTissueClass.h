/**	\file	PCPTissueClass.h
*	Header of PCPTissue class.
*	Contains data and function declarations associated with PCPTissue objects.
*	Depends on base class Tissue.
*/

// include guard
#ifndef __PCPTISSUECLASS_H_INCLUDED__
#define __PCPTISSUECLASS_H_INCLUDED__

//Forward declared dependences
class PCPCell;
class PCPEdge;
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
#include "PCPCellClass.h"
#include "PCPEdgeClass.h"
#include "PCPTissueClass.h"
#include "PCPConstants.h"
#include <random>
#include <chrono>								//temporary, for measuring time and debugging evolve().
//Changes by Hayden for PCP: Now I am going to use a dirichlet random number generation from gsl for
//initialize PCP concentrations on the cell edges.
#include "gsl/gsl_rng.h"
#include "gsl/gsl_matrix.h"						//GSL MATRIX FUNCTIONALITY
#include "gsl/gsl_blas.h"						//GSL BASIC LINEAR ALGEBRA WITH MATRICES
#include "gsl/gsl_vector.h"						//GSL VECTOR FUNCTIONALITY


using namespace std;

//Now, include vector (multidimensional) to store g, f, y for the past two steps (access these within the PCP evolution function). Also, store f and past state for the mechanical variables.
struct pcpparameters{
	int dimension_protein;				//!< Number of dynamical variables associated with protein concentrations.
	int EvalCountProtein;				//!< keep track of PCP evolution evaluation counts
	int dim_mechanical;					//!< keep track of number of vertices + 1 (for the box size)
};

//Form a combined struct for use in Evolve.
struct combinedparameters : parameters, pcpparameters {};

//! Class for PCP tissue object.
class PCPTissue : public Tissue{
public:
	//ADDITIONALDATA
	//keep track of the step sizes, implicit, explicit, and variables from the past two steps.
	vector<vector<double>> g_step_before;
	vector<vector<double>> f_step_before;
	vector<vector<double>> y_step_before;
	vector<vector<Vect>> f_mech_step_before;
	vector<vector<Vect>> past_mech_state;
	vector<double> past_step_sizes;

	vector<vector<double>> g_step_before_save;
	vector<vector<double>> f_step_before_save;
	vector<vector<double>> y_step_before_save;
	vector<vector<Vect>> f_mech_step_before_save;
	vector<vector<Vect>> past_mech_state_save;
	vector<double> past_step_sizes_save;
	//CONVENTION:
	//before a step is taken, [0] store info from two steps before, [1] stores info from one step before.
	//after a step is taken, [0] store info from two steps before, [1] stores info from one step before, [2] stores information from the step that just happened.
	//then, within the imex_solver_adaptive, [1]->[0], [2]->[1] such that the information is updated for the next step. 

	//FUNCTIONS
	//CONSTRUCTOR
	PCPTissue(){};							//!< Default constructor. 
	virtual ~PCPTissue(){};

	//! Creates a new edge and addes it to the list of edges, returns a pointer to that edge: initialize PCP protein concentrations to zero.
	virtual Edge* CreateEdge(Edge* subclass, bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2);

	//Changes by Hayden for PCP: Cut one line of cells and initialize them as differentiated cells. This is called within Main in order to start the process of propagation of differentiation.
	//! Initialize PCP concentrations within a line of cells in the tissue.
	virtual void AddDifferentiatedCellColumn();

	//! This function calculates the total number of edges in differentiated cells in the tissue (half of the number of dynamical protein variables)
	int Find_Number_Differentiated_Edges();

	//! This function calculates the total number of differentiated cells in the tissue (half of the number of dynamical protein variables)
	int Find_Number_Differentiated_Cells();

	//! This is a temporary function used for allowing the input of concentrations by the user.
	//! Reads and constructs simulation from file. The file must be formatted like InputTemplate.txt. This function parses the input file then initializes all the objects and pointers.
	virtual bool ReadInputData(string InputFileName, bool includeConcentrations);

	//! This function calculates the PCP free energy of a configuration (stores the free energy calculation in the location pointed to by protein_energy).
	double CalculatePCPFreeEnergy(void* params);

	//! Integrate the mechanical dynamics using the explicit implementation of the imex scheme (second-order).
	void Integrate_IMEX_Second_Order(double this_step_size, void* params, ofstream* errorfile);

	//! Integrate the mechanical dynamics using the explicit implementation of the imex scheme (third-order).
	void Integrate_IMEX_Third_Order(double this_step_size, void* params, ofstream* errorfile);

	//! Umbrella function for integrating the mechanical dynamics (immediately handles first-second-third order subtleties).
	void Integrate_IMEX_Without_Protein(double* step_size, const double mech_abs_tolerance, double* current_time, void* params, ofstream* errorfile);

	//! Function called by Integrate_IMEX_Without_Protein for adaptive step size.
	bool Integrate_IMEX_Adaptive(double* step_size, const double error_abs_tolerance, void* params, double* t, double consistency_order, ofstream* errorfile);

	//! Do evolution by integrating forward the system of ODE from t=0 to t=timeInterval.
	//! Stepsize is the inital step size of the integrator. The function will modify both stepsize and currentTime for use by subsequent call to Evolve(), useful for breaking integration into many smaller intervals. 
	virtual void Evolve(double timeInterval, double* prot_stepSize, double* mech_stepSize, ofstream* errorfile);
	
	//HELPER FUNCTIONS

	//! Helper function for checking the cell input conforms to expectations.
	virtual void HelperReadInputData_CheckTypes();

	//! Moves the vertices to their new positions and updates the approproate flags if they cross a boundary.
	virtual bool HelperEvolve_UpdatePosition(void* params);

	//! Step size control mechanical evolution
	virtual bool Double_Step_Size_Control(const vector<Vect> &y_error, const double error_tolerance_absolute, double* initial_step_size, double consistency_order);

	//Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration
	//Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors
	//returns the magnutude of Force-Forcenewedge and changes the value of Force
	//NOTE: I HAVE INCLUDED THIS FUNCTION IN CASE WE INSERT DIFFERENT CRITERIA FOR STABILITY BASED ON PCP
	double HelperResolveVertices_CheckStability(Vertex* v, Edge* e1, Edge* e2, Edge* e3, Edge* e4, Cell* c1, Cell* c2, Cell* c3, Cell* c4, Vect *Force);

	//! This function is called by HelperResolveVertices_CheckStability and by HelperResolveVertices_ResolveFourfoldVertex. It redistributes proteins to the new edge according to the rules documented in the 2012 PLOS paper.
	void HelperResolveVertices_RedistributeProteinsNewEdge(PCPEdge* newedge, PCPEdge* edge1, PCPEdge* edge2, PCPEdge* edge3, PCPEdge* edge4);
	
	//! This function is called by HelperT1RemoveEdge. It redistributes proteins from the edge to be deleted to nearby edges as described in the PLOS paper.
	void HelperT1_RedistributeProteinsDeletedEdge(PCPEdge* deledge, vector<PCPEdge*> newedgelist);
	
	//! Writes data into file specified by 'datafile'. If flag 'includeExtraInfo' is true, then additional information (such as tension, cell area, etc) is included.
	//! The extra information is non-essential in reconstructing the simulation state, and is just provided for convenience (e.g. for display in MATLAB).
	virtual void HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo);

	//! This is an umbrella function for handling the transition from first-second-third order evolution of protein concentrations.
	void HelperEvolve_IMEX(vector<double> &y_protein, double* prot_stepSize, double* mech_stepSize, const double prot_abs_tolerance, const double mech_abs_tolerance, double* current_time, void* params, int method_order, ofstream* errorfile);

	//! Initialize protein vector of log variables (for use within evolve).
	void HelperEvolve_Initialize_Protein_Vector(vector<double> &y);

	//! Update protein variables after integration of the dynamics.
	void HelperEvolve_Update_Protein_Variables(const vector<double> &y);

	//! Called by Evolve to renormalize the protein amounts in all cells. Lagrange multiplier constrains normalization, but we must correct for numerical drift.
	void HelperEvolve_RenormalizePolarityProteins();

	//! This function is used to loop through the list of undifferentiated cells that are in contact with differentiated cells at a propagation time.
	//! Allows for delays to fix defects.
	virtual void HelperEvolve_GeneratePCPProtein(bool delay);

	//! Function for allowing generation of PCP from defect list (delay).
	void HelperEvolve_GeneratePCPProtein_Defect();

	//! \note will have to adapt this such that we can pick out only one of the boundaries for differentiation.
	//! This function is used to locate cells on the boundary between differentiated and undifferentiated cells.
	vector<PCPCell*> HelperEvolve_FindBoundaryCells();
	vector<PCPCell*> HelperEvolve_FindCells_Contact_with_Boundary();

	//! Called by HelperResolveVertices_ResolveFourfoldVertex to update the pointers when the vertex needs to resolve, Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors 
	virtual Edge* HelperResolveVertices_UpdatePtrs(Vertex* vert, Edge* edge1, Edge* edge2, Edge* edge3, Edge* edge4, Cell* cell1, Cell* cell2, Cell* cell3, Cell* cell4, Vect* Force);

	//! Reads the data and sets all of the global variables.
	//! Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
	virtual void HelperReadInputData_SetGlobalValues(ifstream* datafile);

	//! Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids.
	//! Starts at current position of *datafile and reads until keyword 'End'.
	virtual void HelperReadInputData_CreateObjects(ifstream* datafile);

	//! RemoveEdge deletes an edge and merges the two vertices to create a fourfold vertex
	virtual Vertex* HelperT1_RemoveEdge(Edge* deledge, ofstream* errorfile, int* nT1);

	//! Initializes all of the pointers so that each object has pointers to it's neighbors.
	virtual void HelperReadInputData_CreatePointers();

	//! This is a temporary function used for allowing the input of concentrations by the user.
	//! Initializes all of the pointers so that each object has pointers to it's neighbors.
	void HelperReadInputData_ProteinConcentrations(ifstream* datafile);

	//! This is a helper function for generating a checkerboard pattern in the next row: for periodic boundary conditions with a vertical stripe rather than in a circular ring.
	virtual void HelperCheckerboard_GeneratePattern(vector<PCPCell*> next_differentiated_cells, ofstream* datafile);

	//! Check for any initial defect (a cell that is entirely surrounded by type intersect_type, type pcp_cell_type).
	virtual bool check_for_faulty_cell(Cell* cell_ptr);
};

//! This builds a list of cells that helps to determine fate specification in the next row.
vector<Cell*> HelperCheckerboard_BuildList_ConnectedCells(vector<PCPCell*> next_differentiated_cells);

//! Build the adjacency matrix (stores the values of contact lengths between cells in contact_matrix).
void HelperCheckerboard_BuildAdjacencyMatrix(vector<PCPCell*> next_differentiated_cells, gsl_matrix* contact_matrix);

//! Generate matrix of random spins for finding a low energy state (anti-ferromagnetic). This is the initial state for the metropolis algorithm.
vector<double> HelperCheckerboard_GenerateRandomSpins(int number_of_spins);

//! Order cell pointer by vertical position (try to identify chains of cells).
vector<PCPCell*> order_cells_position(const vector<PCPCell*> &unordered_cells);

//! This function searches through edge_candidates (edges along the differentiation boundary).
//! Output an edge from this list that faces to the right. This is used within AddDifferentiatedCellColumn (to establish front propogation to the right).
Edge* edge_candidate_to_right(const vector<Edge*> &edge_candidates);

//! Calculate the variance from the history vector.
double variance_history(const vector<double> &test_vector);

//! Find domains of spins (cell-cell type conflicts).
vector<int> find_domains_of_spin(const vector<double> &spin_vector);

//! Check if test cell is within list_of_cells
bool check_for_cell_pointer_within_list(PCPCell* test_cell, vector<PCPCell*> &list_of_cells);

#endif