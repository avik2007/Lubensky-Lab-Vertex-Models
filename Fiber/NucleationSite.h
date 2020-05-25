// include guard
#ifndef __NUCLEATIONSITE_H_INCLUDED__
#define __NUCLEATIONSITE_H_INCLUDED__

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


//! Class for NucleationSite objects.
//Fiber verticies are not a subclass of regular verticies they are a struct with aditional info that follow their parent vertex around
//Why? each parent vertex has 1-3 fiberverticies associated with it so we want to keep the origonal verticies 
class NucleationSite{
public:
	//DATA
	int id;
	FiberCell* cell;
	FiberEdge* edge1; 
	FiberEdge* edge2; 
	Vertex* parentvertex; 

	//CLONE
	NucleationSite* Clone(){ return new NucleationSite(*this); }

	//CONSTRUCTORS AND DESTRUCTORS
	NucleationSite(){};
	virtual ~NucleationSite(){};
	NucleationSite(unsigned long int idnumb, Vertex* parentv, FiberEdge* e1, FiberEdge* e2, FiberCell* c){
		id = idnumb;
		edge1= e1;
		edge2 = e2; 
		cell = c;
		parentvertex = parentv; 
	}

	//Finds the vector that bisects the angle formed by the edges
	// morally equivalent to GetOutwardAngle
	Vect GetOutwardDirection();

	//Finds the vector that bisects the angle formed by the edges
	// tells me the angle of bisecting vector relative to midline
	// midline along x-axis
	// close to cos(theta) == 1 if parallel to midline
	double GetOutwardAngle();

	// cos 2n of outward angle is used to characterize fiber nucleation
	double Calculate_Cos2OutwardAngle();

	//Finds the angle between edge 1 and edge 2
	// tells me opening angle
	double GetAngle();

	//Finds the position of the vertex
	Vect GetLocation(){return (*parentvertex).r;}

	//Returns the probability of birthing a fiber in 1 sec (checks if nucleating site has appropriate opening angle).
	double ProbFiberBirth(); 

	//Returns the probability of fiber death (added by Hayden)
	double ProbFiberDeath();

	double Calc_XDist_toDeathVertex();
};

#endif