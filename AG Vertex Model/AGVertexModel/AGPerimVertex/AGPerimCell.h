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
};
#endif