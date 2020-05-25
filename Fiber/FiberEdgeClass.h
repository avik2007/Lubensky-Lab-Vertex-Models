// include guard
#ifndef __FIBEREDGECLASS_H_INCLUDED__
#define __FIBEREDGECLASS_H_INCLUDED__

//Forward declared dependences
class Fiber;
class FiberCell;
class FiberEdge;
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
#include <unordered_map>
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h" 

using namespace std;


//! Class for cell objects.
class FiberEdge : public Edge{
public:
	//ADDITIONAL DATA
	vector<Fiber*> fibers; 
	double color; //color for drawing output 
	
	//ADDITIONAL FUNCTIONS

	//CLONE
	FiberEdge* Clone(){
		return new FiberEdge(*this);
	}

	//CONSTRUCTORS AND DESTRUCTORS
	FiberEdge(){};
	virtual ~FiberEdge(){};
	FiberEdge(unsigned long int idnum, bool flagx, bool flagy, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2){
		id = idnum;
		yflag = flagy;
		xflag = flagx;
		vertex1 = v1;
		vertex2 = v2;
		cell1 = c1;
		cell2 = c2;
		color = 0; 
	}

	//Returns true if the edge has a fiber
	bool HasFiber();

	//Removes the fiber from the edge
	void RemoveFiber(Fiber* f); 
};

#endif