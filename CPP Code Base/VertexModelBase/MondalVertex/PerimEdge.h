// include guard
#ifndef __PERIMEDGE_H_INCLUDED__
#define __PERIMEDGE_H_INCLUDED__

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


#include "Constants.h"
#include "ConstantsPerim.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h" 
#include "PerimCell.h"

using namespace std;

class Edge_Perim : public Edge{
public:
	//ADDITIONAL DATA (none)
	
	//CONSTRUCTORS AND DESTRUCTORS
	Edge_Perim(){};
	virtual ~Edge_Perim(){};
	Edge_Perim(unsigned long int idnum, bool flagx, bool flagy, Vertex* v1, Vertex* v2){
		id = idnum;
		yflag = 0;
		xflag = 0;
		vertex1 = v1;
		vertex2 = v2;
	}

	//Returns the magnitude of the tension on the edge 
	double GetTension_Magnitude();
};

#endif