// define guard
#ifndef __AGPERIMCELL_H_INCLUDED__
#define __AGPERIMCELL_H_INCLUDED__

// forward class dependencies

class Cell;
class Vertex;
class Vect;
class Edge;

//includes 
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <random>

#include "AGCell.h"


using namespace std;
// AGCell = developed for anisotropically growing models.
class AGPerimCell : public AGCell{
public:
	//***************************************************************************************************************
	//********** Additional Data ************************************************************************************
	//***************************************************************************************************************
	double P0; // preferred perimeter
	
	// additional functions
	virtual double GetEnergy();
	// calculates the individual energy for a cell in the area model

	virtual double GetdPdLx(); // gets the derivative of the cell's area w.r.t LX

	virtual double GetdPdLy(); // gets the derivative of the cell's area w.r.t. LY

	virtual double Getd2PdLxdy(Vertex* v);
	// gets the derivative of cell area with respect to gamma and a vertex y-coordinate (the derivative with respect to gamma and a vertex x-coordinate would be 0)

	//virtual double GetdPdgammaForSimpleShearAlongXAxis();
	// gets the derivative of the cells area w.r.t. gamma if gamma is the shear strength for a simple shear in the x-direction

	//virtual double Getd2PdgammadyForSimpleShearAlongXAxis(Vertex* v);
	// gets the derivative of the cells area w.r.t. gamma and a vertex y coordinate if gamma is the shear strength for a simple shear in the x-direction
};
#endif