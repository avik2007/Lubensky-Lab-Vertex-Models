/**	\file	PCPMegaCellClass.h
*	Header of PCPMegaCell class.
*	Contains data and function declarations associated with PCPMegaCell objects.
*	Depends on base class Cell and on subclass PCPCell.
*/

// include guard
#ifndef __PCPMEGACELLCLASS_H_INCLUDED__
#define __PCPMEGACELLCLASS_H_INCLUDED__

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

//Base class files
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h" 

//PCP subclass files
#include "PCPEdgeClass.h"
#include "PCPCellClass.h"
#include "MegaCellConstants.h"

//Changes by Hayden for PCP: Now I am going to use a dirichlet random number generation from gsl for
//initialize PCP concentrations on the cell edges.
#include "gsl\gsl_rng.h"			//Used of random generation of initial PCP protein concentrations
#include "gsl\gsl_randist.h"		//Used of random generation of initial PCP protein concentrations

using namespace std;


//! Class for cell objects.
class PCPMegaCell : public PCPCell{
public:
	//ADDITIONAL FUNCTIONS
	//CONSTRUCTORS AND DESTRUCTORS
	PCPMegaCell(){};
	virtual ~PCPMegaCell(){};

	//! Returns TRUE if the cell crosses the X boundary, FALSE otherwise.
	bool GetXFlag();

	//! Returns TRUE if the cell crosses the Y boundary, FALSE otherwise.
	bool GetYFlag();

	//! When a cells gets PCP, start clocks in neighboring cells.
	virtual void Helper_AddBorderingCellsPoisson();

	//! This function is used to calculate the distance from the origin of PCP (centroid of central megacell).
	virtual double calculate_Distance_from_PCP_Origin();

	//! Calculates area of a polygon given position of vertices in CW order.
	//! Area returned may be positive or negative. Negative area means that the vertices were ordered CCW.
	double HelperArea(vector<double> &xpos, vector<double> &ypos);

	//! Same as GetPressure() but with cell area as input.
	// \note This contains the implementation of constant inner pressure within the central mega cell and constant pressure from the outer cell.
	virtual double GetPressure(double Area){
		//Central mega cell is type center_cell_type (constant inner_pressure).
		if (type == center_cell_type)
			return inner_pressure;
		//Outer mega cell is type outer_cell_type (constant outer_pressure).
		else if (type == outer_cell_type)
			return outer_pressure;
		//Return K * (A0 - Area) for other cells.
		else
			return K * (A0 - Area);
	}
};

#endif