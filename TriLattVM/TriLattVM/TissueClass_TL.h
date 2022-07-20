/**	\file	TissueClass.h
*	Header of Tissue class.
*	Contains data and function declarations associated with Tissue object.
*/

/*This is the header file for the vertex model. It contains the definition of the tissue class */

// include guard
#ifndef __TISSUECLASS_H_INCLUDED__
#define __TISSUECLASS_H_INCLUDED__

//Forward declared dependences
class Edge;
class Cell;
class Vertex;
class Vect;

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include "VectClass_TL.h"
#include <random>
#include <chrono>								//temporary, for measuring time and debugging evolve().

using namespace std;

struct parameters{
	std::unordered_map<long, int> vert_id2index;		//!< Look-up table (map) to retrieve vertex index given its ID. Faster than finding using for-loops.
	std::unordered_map<long, int> cell_id2index;		//!< Look-up table (map) to retrieve cell index given its ID. Faster than finding using for-loops.
};

//! Class for tissue object.
class Tissue{
public:
	//DATA

	//Create 3 big lists to store all of the edges vertices and cells in
	list<Vertex*> vertices;								//!< List of all vertex objects (not pointers!) in the simulation. This is where the real objects are stored.			
	list<Edge*>  edges;									//!< List of all edge objects (not pointers!) in the simulation. This is where the real objects are stored.	
	list<Cell*> cells;									//!< List of all cell objects (not pointers!) in the simulation. This is where the real objects are stored.	
	list<Cell*> OGcells;								//!< List of initial cells in packing

	double LXold;										//Place to store changing boxsize info during evolve
	double LYold;
	double ForceLX;
	double ForceLY;
	list<int> vertexIDs;								//lists of all of the ID lables ever used for the vertices edges and cells so new lables can be made unique
	list<int> edgeIDs;
	list<int> cellIDs;
	double A0;

	//FUNCTIONS
	//CONSTRUCTOR
	Tissue(){};							//!< Default constructor. \todo Remove!
	virtual ~Tissue(){};						//!< destructor: does not need to be called by you unless you created the object using the word "new" otherwise the compiler will take care of this for you \todo Remove! 

	//! Creates a new vertex and adds it to the list of vertices, returns a pointer to that vertex.
	virtual Vertex* CreateVertex(Vect position, vector<Edge*> edgs);

	//! Removes vertex 'vert' from the master list.
	void DestroyVertex(Vertex* vert);

	//! Creates a new edge and addes it to the list of edges, returns a pointer to that edge.
	virtual Edge* CreateEdge(bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2);

	//! Removes the vertex from the master list.
	void DestroyEdge(Edge* eptr);

	//! Creates a new cell and addes it to the list, returns a pointer to that cell.
	Cell* CreateCell(int typ, vector<Edge*> edgs);

	//! Removes the cell from the master list.
	void DestroyCell(Cell* cptr);

	//! Do mitosis on cells. Picks cells for mitosis based on certain criteria and calls HelperMitosis(Cell * parentcell) to perform mitosis on those cells.
	bool Mitosis();

	//! Looks thourgh all of the cells in the tissue and removes any that are to small, returns TRUE if at least 1 apoptosis is successful.
	virtual bool Apoptosis(ofstream* errorfile);

	//! Looks through all of the edges and find any that are to small and have them undergo T1 transitions, returns TRUE if successful.
	bool T1(ofstream* errorfile);

	//! Looks through all the fourfold vertices and splits them along the lines of force if they are unstable (or keeps them if they are stable), returns TRUE if any changes happened.
	bool ResolveVertices();

	//! Reads and constructs simulation from file. The file must be formatted like InputTemplate.txt. This function parses the input file then initializes all the objects and pointers.
	virtual void ReadInputData(string InputFileName);

	//! Similar to WriteOutputData(), but with additional informations (cell area, pressure, etc) appended to create interactive/clickable plot in MATLAB. Actual implementation is in HelperWriteOutputData_Write2File().
	void WriteOutputData4Matlab(string fileName);

	//! Similar to WriteOutputData(), but writes/appends to an existing file instead of opening a new one. Useful for output over regular interval, e.g. for making movie in MATLAB.
	void WriteOutputDataAppend(ofstream* datafile);

	//! Calculates the force on every vertex and stores the result in vertex.
	virtual void ForceFunction(ofstream* errorfile);

	//! Integrates the force function by a 4th order Runge-Kutta algorithem 
	void Integrate_RK4(double dt, ofstream* errorfile);

	//! Adaptive evolver for mechanical variables (enforcing given error tolerances). Only evolve to time t1.
	int Integrate_RK4_Adaptive(double* step_size, const double error_abs_tolerance, double* t, double t1, ofstream* errorfile);

	//! Integrates the forcefunction by one euiler step 
	void Integrate_Euiler(double dt, ofstream* errorfile);

	//! Evolves the tissue one time interval taking steps no larger than dt and updates the vertex positions
	virtual void Evolve(double interval, double* stepsize, ofstream* errorfile);

	//! Calculate Tissue Stress in both directions for evolving step size.
	virtual double CalculateTissueStressX(double cutline, ofstream* errorfile);
	virtual double CalculateTissueStressY(double cutline, ofstream* errorfile);

	//HELPER FUNCTIONS
	//Helper frunctions: these functions should not be called in the main program but exist to support the primary functions above

	//Integrates the forcefunction by one euiler step 
	void RK4Helper_EuilerStep(double dt, ofstream* errorfile);
	
	//! Reads the data and sets all of the global variables.
	//! Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
	virtual void HelperReadInputData_SetGlobalValues(ifstream* datafile);

	//! Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids. Starts at current position of *datafile and reads until keyword 'End'.
	virtual void HelperReadInputData_CreateObjects(ifstream* datafile);
	 
	//! Initializes all of the pointers so that each object has pointers to it's neighbors.
	virtual void HelperReadInputData_CreatePointers();

	//! Writes data into file specified by 'datafile'
	virtual void HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo);

	//! RemoveEdge deletes an edge and merges the two vertices to create a fourfold vertex
	virtual bool HelperT1_RemoveEdge(Edge* deledge);

	//! Removes a cell from the tissue and creates a new vertex at its center which is attached to the propper edges.
	//! Returns false if cell is not deleted due to 5-fold vertices, true otherwise.
	virtual bool HelperApoptosis_DeleteCell(Cell* delcell);

	//! Return look-up table to retrieve vertex index (i.e., position in master list), given its ID.
	//! Since looking up value from an unordered map is MUCH faster than finding by for-loops, using this function is highly-recommended if you need lots of 'find' operations.
	//! Don't forget to re-create table after every time list vertices change (e.g. after topological transition).
	//! \note Accessing using at() operator is recommended over [] operator, since [] might trigger rehash and won't return exception for 'invalid' key (see C++ documentation).
	unordered_map<long, int> CreateVertexId2IndexMap();

	//! Return look-up table to retrieve cell index (i.e., position in master list), given its ID.
	//! Since looking up value from an unordered map is MUCH faster than finding by for-loops, using this function is highly-recommended if you need lots of 'find' operations.
	//! Don't forget to re-create table after every time list vertices change (e.g. after topological transition).
	//! \note Accessing using at() operator is recommended over [] operator, since [] might trigger rehash and won't return exception for 'invalid' key (see C++ documentation).
	unordered_map<long, int> CreateCellId2IndexMap();

	//! Return look-up table to retrieve edge index (i.e., position in master list), given its ID.
	//! Since looking up value from an unordered map is MUCH faster than finding by for-loops, using this function is highly-recommended if you need lots of 'find' operations.
	//! Don't forget to re-create table after every time list vertices change (e.g. after topological transition).
	//! \note Accessing using at() operator is recommended over [] operator, since [] might trigger rehash and won't return exception for 'invalid' key (see C++ documentation).
	unordered_map<long, int> CreateEdgeId2IndexMap();

	//! Moves the vertices to their new positions and updates the approproate flags if they cross a boundary.
	virtual bool HelperEvolve_UpdatePosition();

	//! Build array of vertex positions (x and y) as well as LX and LY based on current values.
	void Store_Initial_Condition(vector<Vect> &y);

	//! Update rold with specified value.
	void Update_Old_Mechanical_Variables(vector<Vect> &y);

	//! Store the updated positions (from any of the mechanical position evolvers).
	void Store_New_Positions(vector<Vect> &y);

	//! Step size control mechanical evolution
	virtual bool Double_Step_Size_Control(const vector<Vect> &y_error, const double error_tolerance_absolute, double* initial_step_size, double consistency_order);

	//! Does the same thing as HelperEvolve_UpdatePosition but to a single vertex.
	void HelperResolve_UpdatePosition(Vertex* vptr);

	//! This is the helper function that actually do mitosis on cells.
	void HelperMitosis(Cell * parentcell);

	//! Resolves a forfold vertex into 2 threefold vertices if it is energeticaly favorable to do so, returns TRUE if the resolve is successful.
	bool HelperResolveVertices_ResolveFourfoldVertex(Vertex* v);

	//! Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration
	//! Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors
	//! returns the magnutude of Force-Forcenewedge and changes the value of Force
	virtual double HelperResolveVertices_CheckStability(Vertex* v, Edge* e1, Edge* e2, Edge* e3, Edge* e4, Cell* c1, Cell* c2, Cell* c3, Cell* c4, Vect *Force);

	//! Called by HelperResolveVertices_ResolveFourfoldVertex to update the pointers when the vertex needs to resolve, Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors 
	virtual void HelperResolveVertices_UpdatePtrs(Vertex* vert, Edge* edge1, Edge* edge2, Edge* edge3, Edge* edge4, Cell* cell1, Cell* cell2, Cell* cell3, Cell* cell4, Vect* Force);

	//! Writes the current state of the simulation to file. The file is formatted like InputTemplate.txt and thus can be readily read by ReadInputData(). Actual implementation is in HelperWriteOutputData_Write2File().
	//! \note This function will only write the essential data, enough for ReadInputData() to reconstruct the simulation state. All other information (e.g. cell area, pressure, etc) can be calculated from the data.
	//! \note For cases where this extra information is needed (e.g. for visualization with MATLAB, debugging, etc) use function WriteOutputData4Matlab(string fileName) instead.
	void WriteOutputData(string fileName);

	//Cycles through and randomly selects cells to change their A0 to zero
	void ApoptosisSelection(double p, Cell* apopcell);

	//Calculates and writes info for density fluctuations of packing 
	void CircleDensityTest_OutAppend(int circ_num, int rad_num, int t, int n, bool logdist, bool vertexCount, ofstream* densityfile);

	//Calculates and writes info for density fluctuations of packing for final time step
	void CircleDensityTest_Out(int circ_num, int rad_num, bool logdist, bool vertexCount, string fileName);

	//Randomly assigns restlengths to each edge
	void AssignRestLengths(double limit_type1, double limit_type2);

	//Randomly assigns spring constants to each edge
	void AssignSpringConstants(double limit_type1, double limit_type2);

	//Assigns A0 based on assigned restlengths (assuming all triangles)
	void CalculateA0();

	//Calculate area inside circle
	void AreaCalc(int circ_num, int rad_num, bool logdist, string fileName, string fileNameError);

	//Calculate area inside circles throughout the packing and record to a file
	void AreaCalc_Append(int circ_num, int rad_num, bool logdist, ofstream* densityfile, int t, int n);

	//Calculate area of a triangle intercepted by circle
	double GetCircleIntercept(vector<Vertex*> vert_out, vector<Vertex*> vert_in, vector<double> Rpoint, double R, bool edgeflag);

	//checks for exceptions with 1 vertex inside circle
	double CheckException(vector<Vertex*> vert_out, vector<Vertex*> vert_in, vector<double> Rpoint, double R, bool edgeflag);

	//returns cells of exception cells that have already been checked in order to avoid coutning area more tjan once
	vector<Cell*> ExceptionCellIDs(vector<Vertex*> vert_out, vector<Vertex*> vert_in, vector<double> Rpoint, double R, Cell* initialcell);

	//returns area of all exceptional cells, and weeds out cells that might be calculated repeatedly
	double ExceptionCalc(vector<Vertex*> vert_out, vector<Vertex*> vert_in, vector<double> Rpoint, double R, Cell* initialcell, vector<Cell*> checked_cells, int flag_areatype);

	//for exceptional triangles, checks to see if any verties without xflag and yflags should be moved
	vector<double> ExceptionalPBC(vector<double> vertpos, vector<double> Rpoint, double R);

	void AreaCalcTEST();

	void vertexNoise(double v_limit);

	bool RelaxedPacking(double tol);

	void IntersectionTest(int t, int n);

	//Assigns A0 based on assigned restlengths (assuming all triangles)
	void CalculateA0Random();

	void RandomTester();
};

#endif