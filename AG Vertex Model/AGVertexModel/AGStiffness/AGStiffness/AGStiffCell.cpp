#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <ctime>
#include <random>
#include <array>
#include <typeinfo>
#include <numeric>
#include <stddef.h>


#include "Constants.h"
#include "AGConstants.h"
#include "TissueClass.h"
#include "AGTissue.h"
#include "AGCell.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "AGPerimTissue.h"
#include "AGPerimConstants.h"
#include "AGPerimCell.h"
#include "AGPerimEdge.h"
#include "VectClass.h"
#include "AGStiffCell.h"
#include "AGStiffTissue.h"


/*
* void SetAreaModulus(double K)
*
* designed to be used by the tissue to set the area modulus of the AGStiffCells
*/
void AGStiffCell::SetAreaModulus(double K)
{
	(*this).KAM = K;
}

/*
* double GetPressure(double area)
*
* returns the pressure of AGStiffCell
*/
double AGStiffCell::GetPressure(double area)
{
	return (*this).KAM * (1 - area / (*this).A0);
}

/*
* double GetEnergy()
*
* returns the energy of an individual AGStiffCell
*/
double AGStiffCell::GetEnergy()
{
	double energy = 0.5 * ((*this).KAM * pow((*this).GetArea() - (*this).A0, 2) / (*this).A0 + LAMBDA * pow((*this).GetPerimeter() - P0, 2));
	return energy;
}