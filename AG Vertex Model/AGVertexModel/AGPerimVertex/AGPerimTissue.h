// guard
#ifndef __AGPERIMTISSUE_H_INCLUDED__
#define __AGPERIMTISSUE_H_INCLUDED__

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

class AGPerimTissue : public AGTissue{
	/*
	* It's important to mention, this class descends from AGTissue but is not really intended to grow.
	*
	* If one would like to make the box grow, it wouldn't be crazy difficult, but you'd have to modify
	* ForceFunction and stress tensor calculation functions found in AGTissue
	*
	*/
public:
	//DATA
	//double stepsize_upperbound; // stepsize of integrator
	//double equilibriumtolerance = .00000000010; // if (energy change) / (timestep) between steps is less than this, it's equilibrium
	double equilibriumtolerance; 
	double A0; // preferred area of cell
	double P0; // preferred perimeter of cell
	list<double> energyList; // contains the tissue's energy from a fixed number of past timesteps
	list<double> timeList; // contains a fixed number of the simulation's past time steps
	//double energyListSize = 50;
	int energyListSize;
	bool equilibrated; //defines whether the tissue is considered equilibrated
	double outputPrecision = 32;
	//CONSTRUCTOR
	AGPerimTissue(){};

	//DESTRUCTOR
	virtual ~AGPerimTissue(){};

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

	//calculates the perimeter energy of a packing
	virtual double CalculatePerimeterEnergy();

	//calculates the pressure energy of a packing, assuming quadratic pressure potential
	virtual double CalculateQuadraticPressureEnergy();
	
	//calculates the pressure energy of a packing, assuming a logarithmic pressure potential
	virtual double CalculateLogPressureEnergy();
	
	//calculates the pressure energy of a packing, assuming an exponential pressure potential
	virtual double CalculateExpPressureEnergy();

	// does the numerical integration for the simulation
	virtual void Evolve(double interval, double *stepsize, ofstream* errorfile);

	//! Creates a new (AGPerim) edge and addes it to the list of edges, returns a pointer to that edge.
	virtual Edge* CreateEdge(Edge* subclass, bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2);

	// checks if tissue is equilibrated
	virtual bool IsEquilibrated();


	//! Looks through all of the edges and find any that are to small and have them undergo T1 transitions, returns TRUE if successful.
	virtual bool T1(ofstream* errorfile);

	//! Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration
	//! Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors
	//! returns the magnutude of Force-Forcenewedge and changes the value of Force
	virtual double HelperResolveVertices_CheckStability(Vertex* v, Edge* e1, Edge* e2, Edge* e3, Edge* e4, Cell* c1, Cell* c2, Cell* c3, Cell* c4, Vect *Force);

	// Do mitosis on cells. Picks cells for mitosis based on certain criteria and calls HelperMitosis(Cell * parentcell) to perform mitosis on those cells.
	virtual bool Mitosis();

	//virtual void DestroyCell(Cell* cptr);
	virtual Cell* CreateCell(int typ, vector<Edge*> edgs, string potential, double tcjstrength, double thetaan, double anstrenght);

	//! This is the helper function that actually do mitosis on cells.
	virtual void HelperMitosis(Cell* parentcell, double divisionangle);

	//Calculates the Hessian of the tissue
	virtual Eigen::MatrixXf getHessian(ofstream* errorfile);

	// 10/25/2023 - I think there is a smarter way to write this function that makes it easier to debug and faster to run
	virtual Eigen::MatrixXf calculateHessian(ofstream* errorfile);
};

double listAverage(list<double> inputList);

double getLinearRegressionSlope(list<double> xList, list<double> yList);

#endif
