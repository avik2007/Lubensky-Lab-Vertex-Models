#pragma once
// define guard
#ifndef __AGSTIFFCELL_H_INCLUDED__
#define __AGSTIFFCELL_H_INCLUDED__

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

#include "AGPerimCell.h"
#include "VectClass.h"


using namespace std;
// AGStiffCell -> allows for cells to have varying stiffnesses
class AGStiffCell : public AGPerimCell {
public:
	//AGStiffCell() {};
	//virtual ~AGStiffCell() {};
	//***************************************************************************************************************
	//********** Additional Data ************************************************************************************
	//***************************************************************************************************************
	double KAM; // (K, Area Modulus; varies per cell now)

	virtual void SetAreaModulus(double K);

	virtual double GetPressure(double Area);

	virtual double GetEnergy();
	// calculates the individual energy for a cell in the area model
};
#endif