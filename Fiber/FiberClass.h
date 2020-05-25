// include guard
#ifndef __FIBERCLASS_H_INCLUDED__
#define __FIBERCLASS_H_INCLUDED__

//Forward declared dependences
class Fiber;
class FiberCell;
class FiberEdge;
class NucleationSite;
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
#include "FiberCellClass.h"

using namespace std;


// Class for fiber objects.
class Fiber{
public:
	//DATA
	unsigned long int id;			//the id of this fiber, (the id is unique to this fiber)
	FiberCell* cell;				//the cell this fiber lives in
	FiberEdge* edge1; 
	FiberEdge* edge2; 
	Vertex* v1;						//fiber travels along edge1 out from v1 
	Vertex* v2;
	double dedge1_x;		//distance along the edge in the x direction
	double dedge2_x;		//distance along the edge in the x direction
	double birthtime;
	bool passedvertex;				//True if either end of the fiber has passed a vertex 
	
	//FUNCTIONS
	//CONSTRUCTORS AND DESTRUCTORS
	Fiber(){};
	virtual ~Fiber(){}; 
	Fiber(unsigned long int i){
		id = i;
	}

	//returns the length of the fiber
	double Length(); 

	//returns the position of the edne of the fiber (on edge one) 
	Vect r1();

	//returns the position of the edne of the fiber (on edge two) 
	Vect r2();

	//returns the edge on the other end of the fiber
	FiberEdge* GetOtherEdge(FiberEdge* e); 

	//Gets the velocity of the fiber 
	double GetVelocity();

	//Gets the appropreate death rate for the fiber 
	double GetBreakProb(); 

	//Returns true if the fiber is on its last edge (ie it will run out of cell when it swaps edges)  
	bool IsLastEdge();

	//Returns the distance of the fiber from the next vertex on Edge 1
	double DistanceToNextVertexonEdge1();

	//Returns the distance of the fiber from the next vertex on Edge 2
	double DistanceToNextVertexonEdge2();
};

#endif