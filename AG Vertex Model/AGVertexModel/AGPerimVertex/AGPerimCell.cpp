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

/*
* virtual double GetEnergy()
*
* returns the energy of the cell
*
*/
double AGPerimCell::GetEnergy()
{
	double energy = 0.5 * ( K * pow((*this).GetArea() - (*this).A0 ,2) / (*this).A0  + LAMBDA * pow( (*this).GetPerimeter() - P0,2) );
	return energy;
}