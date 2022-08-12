/**	\file	AGPerimEdge.cpp
*	Implementation of AGPerimEdgeCell (anisotropically growing tissue edge with Manning Hamiltonian)class.
*	Contains function definitions associated with AGPerimCell objects.
*/

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            #include <cassert>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <algorithm>
#include <stddef.h>

#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "VectClass.h"
#include "Constants.h"
#include "AGConstants.h"
#include "AGCell.h"
#include "AGPerimCell.h"
#include "AGPerimConstants.h"
#include "AGPerimEdge.h"

using namespace std;

//***************************************************************************************************************
//*********** FUNCTIONS FOR CALCULATING EDGE PROPERTIES **********************************************************************
//***************************************************************************************************************

// Returns the magnitude of the tension on the edge
//Should return gamma+Lambda(cell1.P0-cell1.perimeter) + Lambda(cell2.P0- cell2.perimeter
double AGPerimEdge::GetTension_Magnitude(){
	double perimdiff1; //difference between the perimeter and P0 for cell1
	double perimdiff2; //difference between the perimeter and P0 for cell2
	double tension;

	//Get pointers to the neighboring cells 
	AGPerimCell* c1 = dynamic_cast<AGPerimCell*>(cell1); //pointer to cell1 (note that we need to cast the cell as the appropreate cell subclass) 
	AGPerimCell* c2 = dynamic_cast<AGPerimCell*>(cell2);

	//Calculate (cell.P0- cell.perim)
	perimdiff1 = (*c1).GetPerimeter()- (*c1).P0;
	perimdiff2 = (*c2).GetPerimeter()- (*c2).P0;

	//Calculate the tension 
	tension = LAMBDA*perimdiff1 + LAMBDA*perimdiff2;
	return tension;
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    