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
#include "AGCell.h"
#include "VectClass.h"


using namespace std;
// AGStiffCell -> allows for cells to have varying stiffnesses, I'M CHANGING UNITS OF K
class AGStiffCell : public AGPerimCell {
public:
	AGStiffCell() {};
	virtual ~AGStiffCell() {};
	//***************************************************************************************************************
	//********** Additional Data ************************************************************************************
	//***************************************************************************************************************
	double KAM; // (K, Area Modulus; varies per cell now)

	double PM; //( LAMBDA, "Perimeter Modulus"; varies per cell now

	virtual double GetPressureNormalStress();

	virtual double GetTensionNormalStress();

	virtual void SetAreaModulus(double K);

	virtual void SetPreferredArea(double PA);

	virtual void SetPM(double PM);

	virtual double GetPressure(double Area);

	virtual double GetEnergy();
	// calculates the individual energy for a cell in the area model
};
#endif