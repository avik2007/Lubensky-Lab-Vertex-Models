/**	\file	VertexClass.h
*	Header of Vertex class.
*	Contains data and function declarations associated with Vertex objects.
*/

// include guard
#ifndef __VERTEXCLASS_H_INCLUDED__
#define __VERTEXCLASS_H_INCLUDED__

//Forward declared dependences
class Edge;
class Cell;

//includes
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include "VectClassClone.h"		// Class for physical vector with magnitude and direction
using namespace std;


/** Class for vertex/node objects.
*	Objects of this class describes vertices/nodes where three or more edges meet.
*/
class Vertex{
public:
	/********/
	/* DATA */
	/********/
	unsigned long int id;				//!< ID of a vertex. Used to identify if two vertices are equal.
	Vect r;								//!< Current x and y position of the vertex.
	Vect rnew;							//!< Future x and y position of the vertex. 
	Vect rold;							// Position of the vertex just before an integration step is taken 
	vector<Edge*> edges;				//!< Vector containing pointers to Edge objects meeting at a vertex.
	Vect SavedInfoForEvolve_Force;		//Place to store the force on the vertex calculated durring evolve
	Vect oldForceInfo;					//Force on Vertex just a step before



	/******************************/
	/* CONSTRUCTORS & DESCTRUCTOR */
	/******************************/
	Vertex(){}
	virtual ~Vertex(){};						//destructor: does not need to be called by you unless you created the object using the word "new" otherwise the compiler will take care of this for you
	//!	Constructs a vertex of id 'idnumb' at position 'x' and 'y'.
	Vertex(unsigned long int idnumb, double x, double y){
		id = idnumb;
		r.x = x;
		r.y = y;
	}

	//!	Constructs a vertex of id 'idnumb' at position 'x' and 'y', and with edges given in 'edg'.
	Vertex(unsigned long int idnumb, double x, double y, vector<Edge*> edg){
		id = idnumb;
		r.x = x;
		r.y = y;
		edges = edg;
	}



	/*******************************/
	/* STANDARD OPERATOR OVERLOADS */
	/*******************************/

	//! Tests if two vertices are the same.
	//! Compares ID of the two vertices and returns 'true' if they are the same, 'false' otherwise.
	bool operator==(const Vertex &) const;

	//! Tests if two vertices are NOT the same.
	//! Compares ID of the two vertices and returns 'true' if they are NOT the same, 'false' otherwise.
	bool operator!=(const Vertex &right) const{
		return !(*this == right);
	}



	/*************/
	/* FUNCTIONS */
	/*************/

	//! Returns the number of edges meeting at the vertex.
	int GetEdgeNumber();

	//! Check if the vertex has exactly 4 edges.
	bool CheckFourfold();

	//! Removes an edge from the list of edges.
	void RemoveEdge(Edge* eptr);

	//! Returns a pointer to the vertex on the other end of the specified edge.
	Vertex* GetOtherVertex(Edge* edg);

	//! Returns a vector of pointers to all vertices connected to this vertex.
	vector<Vertex*> GetAllOtherVertex();

	//! Returns a pointer to the edge between 'this' vertex and the specified vertex.
	Edge* GetEdgeBetween(Vertex* vrt);

	//! Returns a vector of pointers to neighboring cells.
	vector<Cell*> GetCells();

	//! GetCWCells returns a list of the cells around the vertex in CW order. The cells have the same starting point as the edges in GetCWEdges if if you walked in a CW circle around the vertex starting at e0 you would see: e0,c0,e1,c1,e2,c2,e3,c3,e0
	vector<Cell*> GetCWCells();

	//! Returns a list of the edges around a vertex in CW or CCW order.
	vector<Edge*> GetOrderedEdges();

	//! Returns a list of edges around a vertex in CW order
	vector<Edge*> GetCWEdges();

	//! Changes all x-flags on the edges atached to the vertex to be the oposite of what they currently are (needed for moving vertices across the boundary).
	void ChangeXEdgeFlags();

	//! Changes all y-flags on the edges atached to the vertex to be the oposite of what they currently are (needed for moving vertices across the boundary).
	void ChangeYEdgeFlags();

	//! TEMP: calculate net force due to pressure on a vertex
	Vect GetNetForce_Pressure();

	//! TEMP: calculate net force due to tension on a vertex
	Vect GetNetForce_Tension();
};


#endif	