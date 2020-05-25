/**	\file	MegaEdgeClass.h
*	Header of PCPEdge class.
*	Contains data and function declarations associated with MegaEdge objects.
*	Depends on base class Edge.
*/

// include guard
#ifndef __PCPEDGECLASS_H_INCLUDED__
#define __PCPEDGECLASS_H_INCLUDED__

//Forward declared dependences
class MegaCell;
class Cell;
class Edge;
class Vertex;
class Vect;

#include "EdgeClass.h"
#include "MegaCellConstants.h"

using namespace std;

//! Class for cell objects.
class MegaEdge : public Edge{
public:
	//Returns the magnitude of the tension on the edge (this is implemented differently in the edge subclasses) 
	virtual double GetTension_Magnitude();
};

#endif