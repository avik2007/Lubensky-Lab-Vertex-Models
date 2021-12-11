// include guard
#ifndef __AGPERIMEDGE_H_INCLUDED__
#define __AGPERIMEDGE_H_INCLUDED__

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
#include "EdgeClass.h"
#include "CellClass.h"
#include "Constants.h"
#include "AGPerimTissue.h"
#include "AGPerimCell.h"

using namespace std;

class AGPerimEdge : public Edge{
public:
	// ADDITIONAL DATA
	bool xflagold;
	bool yflagold;									//!< Flag that indicates if the edge wraps around the boundary (periodic boundary condition/PBC) on the x- and y-direction.
	//CONSTRUCTORS
	AGPerimEdge(){};
	AGPerimEdge(unsigned long int idnum, bool flagx, bool flagy, Vertex* v1, Vertex* v2)
	{
		id = idnum;
		yflag = 0;
		xflag = 0;
		vertex1 = v1;
		vertex2 = v2;
	}

	//DESTRUCTORS
	virtual ~AGPerimEdge(){};

	//ADDITIONAL FUNCTIONS
	virtual double GetTension_Magnitude();

};
#endif 