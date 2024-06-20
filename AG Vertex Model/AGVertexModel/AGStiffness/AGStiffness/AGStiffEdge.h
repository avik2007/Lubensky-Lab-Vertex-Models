// include guard
#ifndef __AGSTIFFEDGE_H_INCLUDED__
#define __AGSTIFFEDGE_H_INCLUDED__

//Forward declared dependences 
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
#include <stddef.h>
#include <unordered_map>

#include "AGPerimConstants.h"
#include "AGStiffConstants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "Constants.h"
#include "AGPerimTissue.h"
#include "AGCell.h"
#include "AGConstants.h"
#include "AGPerimCell.h"
#include "AGStiffCell.h"
#include "AGStiffTissue.h"
#include "AGPerimEdge.h"

using namespace std;

class AGStiffEdge : public AGPerimEdge {
public:
	// ADDITIONAL DATA								//!< Flag that indicates if the edge wraps around the boundary (periodic boundary condition/PBC) on the x- and y-direction.
	//CONSTRUCTORS
	AGStiffEdge() {};
	AGStiffEdge(unsigned long int idnum, bool flagx, bool flagy, Vertex* v1, Vertex* v2)
	{
		id = idnum;
		yflag = 0;
		xflag = 0;
		vertex1 = v1;
		vertex2 = v2;
	}

	//DESTRUCTORS
	virtual ~AGStiffEdge() {};

	//ADDITIONAL FUNCTIONS
	virtual double GetTension_Magnitude();

};
#endif 