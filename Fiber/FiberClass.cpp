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
#include "Constants.h"
#include "FiberConstants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "FiberCellClass.h"
#include "FiberClass.h"
#include "FiberEdgeClass.h"
#include "FiberTissueClass.h"
#include "NucleationSite.h"

using namespace std;

//Gets the length of the fiber
double Fiber::Length(){
	Vect r1 = (*this).r1(); 
	Vect r2 = (*this).r2(); 
	return sqrt((r1.x - r2.x)*(r1.x - r2.x) + (r1.y - r2.y)*(r1.y - r2.y));
}

//returns the position of the edne of the fiber (on edge one) 
Vect Fiber::r1(){
	Vect r;
	Vect u1 = (*edge1).MakeUnitVector(v1);

	// we require that the endpoint of the fiber sit on the edge
	// which means we must be adding some multiple of u1 to (*v1).r

	// the endpoint of the fiber should be moved by dedge1_x in the x_direction from v1
	double factor_for_multiplying_u1 = dedge1_x / (u1.x);
	// take the absolute value in order to avoid weird sign issues.
	factor_for_multiplying_u1 = fabs(factor_for_multiplying_u1);

	// start by setting r to the position of v1.
	r = (*v1).r;
	// then move the endpoint along the direction of the edge.
	r.x = r.x + factor_for_multiplying_u1 * u1.x;
	r.y = r.y + factor_for_multiplying_u1 * u1.y;

	return r; 
}

//returns the position of the edne of the fiber (on edge two) 
Vect Fiber::r2(){
	Vect r; 
	Vect u2 = (*edge2).MakeUnitVector(v2);

	// we require that the endpoint of the fiber sit on the edge
	// which means we must be adding some multiple of u2 to (*v2).r

	// the endpoint of the fiber should be moved by dedge1_x in the x_direction from v1
	double factor_for_multiplying_u2 = dedge2_x / (u2.x);
	// take the absolute value in order to avoid weird sign issues.
	factor_for_multiplying_u2 = fabs(factor_for_multiplying_u2);

	// start by setting r to the position of v1.
	r = (*v2).r;
	// then move the endpoint along the direction of the edge.
	r.x = r.x + factor_for_multiplying_u2 * u2.x;
	r.y = r.y + factor_for_multiplying_u2 * u2.y;
	
	return r; 
}

//returns the edge on the other end of the fiber
FiberEdge* Fiber::GetOtherEdge(FiberEdge* e){
	if ((*e).id == (*edge1).id) { return (&(*edge2)); }
	else{ return (&(*edge1)); }
}

double Fiber::GetVelocity(){
	if (!passedvertex) return INITIALVEL; 
	else return SLOWVEL; 
}


//Gets the appropreate death rate for the fiber 
double Fiber::GetBreakProb(){
	if (!passedvertex) return INITIALDEATHRATE;
	else return PASSEDDEATHRATE;
}


//Returns true if the fiber is on its last edge (ie it will run out of cell when it swaps edges)  
bool Fiber::IsLastEdge(){
	Edge *nexte1, *nexte2; 
	Vertex *nextv1, *nextv2;
	bool b=false; 

	//Get the next vertex and edge 
	nextv1 = (*edge1).GetOtherVertex(v1); 
	nextv2 = (*edge2).GetOtherVertex(v2);
	nexte1 = (*nextv1).GetOtherEdgeInCell(cell, edge1); 
	nexte2 = (*nextv2).GetOtherEdgeInCell(cell, edge2);

	//Return true if the next vertex on each edge is the same
	if ((*nextv1).id == (*nextv2).id) b = true; 

	//Return true if the next edge is the same 
	if ((*nexte1).id == (*nexte2).id) b = true;

	return b; 
}

//Returns the (signed) distance of the fiber from the next vertex (+ means you haven't crossed the next vertex yet on this edge).
double Fiber::DistanceToNextVertexonEdge1(){
	// calculate current position of vertex
	Vect r_fib = (*this).r1();
	// find next vertex that the fiber will run into along this edge
	Vertex* nextv1 = (*edge1).GetOtherVertex(v1);
	// store the position of the next vertex along this edge
	Vect nextv1_pos = (*nextv1).r;
	// store unit vector from v1 to nextv1
	Vect u1 = (*edge1).MakeUnitVector(v1);

	// calculate the distance between the current fiber endpoint and the next vertex on this edge
	Vect diff_vect_vertex_fib = nextv1_pos - r_fib;
	double distance_vertex_fib = diff_vect_vertex_fib.DotProduct(u1);

	return distance_vertex_fib;
}

//Returns the (signed) distance of the fiber from the next vertex (+ means you haven't crossed the next vertex yet on this edge).
double Fiber::DistanceToNextVertexonEdge2(){
	// calculate current position of vertex
	Vect r_fib = (*this).r2();
	// find next vertex that the fiber will run into along this edge
	Vertex* nextv2 = (*edge2).GetOtherVertex(v2);
	// store the position of the next vertex along this edge
	Vect nextv2_pos = (*nextv2).r;
	// store unit vector from v2 to nextv2
	Vect u2 = (*edge2).MakeUnitVector(v2);

	// calculate the distance between the current fiber endpoint and the next vertex on this edge
	Vect diff_vect_vertex_fib = nextv2_pos - r_fib;
	double distance_vertex_fib = diff_vect_vertex_fib.DotProduct(u2);

	return distance_vertex_fib;
}