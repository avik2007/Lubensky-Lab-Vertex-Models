/**	\file	VertexClass.cpp
*	Implementation of Vertex class.
*	Contains function definitions associated with Vertex objects.
*/

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <cassert> 
#include <algorithm>
#include <stddef.h>

#include "EdgeClass_TL.h"
#include "CellClass_TL.h"
#include "VertexClass_TL.h"

using namespace std;

// Tests if two vertices are the same.
// Compares ID of the two vertices and returns 'true' if they are the same, 'false' otherwise.
bool Vertex::operator==(const Vertex &right) const{
	//check that the id's are the same
	if (id == right.id)
		return true;
	else
		return false;
}

// Returns the number of edges meeting at the vertex.
int Vertex::GetEdgeNumber(){
	return edges.size();
}

// Check if the vertex has exactly 4 edges.
bool Vertex::CheckFourfold(){
	if (edges.size() == 4)
		return 1;
	//else if ( edges.size() >= 4) 
	//cout << "Vertex " << (*this).id << " is a five fold vertex!" << endl;
	else
		return 0;
}

// Removes an edge from the list of edges.
void Vertex::RemoveEdge(Edge* eptr){
	unsigned int i;
	for (i = 0; i<edges.size(); i++){
		if ((*edges[i]) == (*eptr)){
			edges.erase(edges.begin() + i);
		}
	}
}

// Returns a pointer to the vertex on the other end of the specified edge.
Vertex* Vertex::GetOtherVertex(Edge *edg){
	long int findid;
	int i;
	vector<long int> edgeids;
	//first check that edge is actually a neighbor of vertex
	findid = (*edg).id;
	for (i = 0; i < (*this).edges.size(); i++){
		edgeids.push_back((*edges[i]).id);
	}

	find(edgeids.begin(), edgeids.end(), findid) != edgeids.end();

	if ((*(*edg).vertex1).id == id)
		return (*edg).vertex2;
	else
		return (*edg).vertex1;


}

// Returns a vector of pointers to all vertices connected to this vertex.
vector<Vertex*> Vertex::GetAllOtherVertex(){
	int N = GetEdgeNumber();
	std::vector<Vertex*> verts(N);		//initialize vector verts of size N
	for (int i = 0; i < N; i++){
		verts[i] = (*edges[i]).GetOtherVertex(this);
	}
	return verts;
}

//GetCWCells returns a list of the cells around the vertex in CW order. The cells have the same starting point as the edges in GetCWEdges if if you walked in a CW circle around the vertex starting at e0 you would see: e0,c0,e1,c1,e2,c2,e3,c3,e0
vector<Cell*> Vertex::GetCWCells(){
	vector<Cell*> orderedcells;
	vector<Edge*> orderededges = (*this).GetCWEdges();
	Edge* startedge = orderededges[0];

	int i = 0;
	Edge* e1 = orderededges[i];			//get the first two edges
	Edge* e2 = orderededges[i + 1];

	do {
		// Find the cell between e1 and e2
		if ((*e1).cell1 == (*e2).cell1)
			orderedcells.push_back((*e1).cell1);
		else if ((*e1).cell1 == (*e2).cell2)
			orderedcells.push_back((*e1).cell1);
		else if ((*e1).cell2 == (*e2).cell1)
			orderedcells.push_back((*e1).cell2);
		else if ((*e1).cell2 == (*e2).cell2)
			orderedcells.push_back((*e1).cell2);

		//update e1 and e2
		i++;
		e1 = orderededges[i%orderededges.size()];
		e2 = orderededges[(i + 1) % orderededges.size()];

	} while ((*e1) != (*startedge));  //while e1 is not the starting edge

	return orderedcells;
}

//GetEdgeBetween returns a pointer to the edge between this vertex and the specified vertex
Edge* Vertex::GetEdgeBetween(Vertex* vrt){
	vector<Edge*> vrtedges;
	unsigned int i, j;

	vrtedges = (*vrt).edges;
	//look through all of the edges for both vertices
	for (i = 0; i<edges.size(); i++){
		for (j = 0; j<vrtedges.size(); j++){
			//when you find a match return it
			if ((*edges[i]).id == (*vrtedges[j]).id){
				return edges[i];
			}
			//if you get to the end and never find a match give an error
			else if (((*edges[i]).id != (*vrtedges[j]).id) && (i == (edges.size() - 1)) && (j == (vrtedges.size() - 1))){
				cout << "ERROR: GetEdgeBetween given an unacceptable vertex" << endl;
			}
		}
	}
}

// Returns a vector of pointers to neighboring cells.
vector<Cell*> Vertex::GetCells(){
	unsigned int i, j;
	vector<Cell*> cells;							//vector of cell pointers used as intermediate storage
	vector<Cell*> returncells;						//return a pointer to a vector

	for (i = 0; i<edges.size(); i++){					//for each neighboring edge
		cells.push_back((*edges[i]).cell1);			//add the first cell to the vector
		cells.push_back((*edges[i]).cell2);			//add the second cell to the vector
	}

	//Now remove the duplicate cells in the vector
	returncells.push_back(cells[0]);				//start by coppying over the first cell so the list to return has an element
	for (j = 0; j<cells.size(); j++){					//for each cell in cells
		if (find(returncells.begin(), returncells.end(), cells[j]) == returncells.end()) {				//if the curent edge is NOT in the list of returncells
			returncells.push_back(cells[j]);														//add it to the return list
		}
	}
	return returncells;
}

// Returns a list of the edges around a vertex in CW or CCW order.
vector<Edge*> Vertex::GetOrderedEdges(){
	vector<Edge*> returnlist;
	Edge e;
	Cell c;
	int i, j;

	returnlist.push_back(edges[0]);		//initialize list
	c = (*(*edges[0]).cell1);

	for (i = 0; i<edges.size() - 1; i++){
		e = (*returnlist.back());      //get the next edge
		if ((*e.cell1) == c)		//get the next cell
			c = (*e.cell2);
		else
			c = (*e.cell1);
		for (j = 0; j<edges.size(); j++){  //look through the list of edges and get the next one that bounds cell c as well
			if (((*(*edges[j]).cell1) == c || (*(*edges[j]).cell2) == c) && (*edges[j]) != e){
				returnlist.push_back(edges[j]);
			}
		}
	}
	return returnlist;
}

//returns a list of edges around a vertex in CW order
vector<Edge*> Vertex::GetCWEdges(){
	//Declare and initialize vertex pointers.
	Vertex* v0 = 0;
	Vertex* v1 = 0;
	vector<Edge*> ordered = Vertex::GetOrderedEdges();	//get the ordered edges
	Edge* e0 = ordered[0];			//edge0
	Edge* e1 = ordered[1];			//edge1
	v0 = (*e0).GetOtherVertex(this);		//vertex on the other side of edge0
	v1 = (*e1).GetOtherVertex(this);		//vertex on the other side of edge1
	Vect vect0 = (*e0).GetVectorDifference(this, v0);	//vector of edge0
	Vect vect1 = (*e1).GetVectorDifference(this, v1);		//vector of edge1
	int sign;

	//take the cross product of vect0 X vect1
	sign = vect0.SignCrossProduct(vect1);

	if (sign == -1){
		return ordered;		//they were alread ordered CW so return the origonal order
	}
	else{
		reverse(ordered.begin(), ordered.end());
		return ordered;
	}
}

//changes all of the flags on the edges atached to the vertex to be the oposite of what they currently are (needed for moving vertices across the boundry)
void Vertex::ChangeXEdgeFlags(){
	int i;
	for (i = 0; i<edges.size(); i++){
		if ((*edges[i]).xflag == 0)
			(*edges[i]).xflag = 1;
		else
			(*edges[i]).xflag = 0;
	}
}

// Changes all y-flags on the edges atached to the vertex to be the oposite of what they currently are (needed for moving vertices across the boundary).
void Vertex::ChangeYEdgeFlags(){
	int i;
	for (i = 0; i<edges.size(); i++){
		if ((*edges[i]).yflag == 0)
			(*edges[i]).yflag = 1;
		else
			(*edges[i]).yflag = 0;
	}
}

//! TEMP: calculate net force due to pressure on a vertex
Vect Vertex::GetNetForce_Pressure(){
	Vect FP(0.0, 0.0);
	Vect zerovec(0.0, 0.0);
	double deltaP, Pmag;
	double angle;
	Edge* eptr;

	//Declare and initialize pointers.
	Vertex* v1 = 0;
	Vertex* v2 = 0;
	Cell* c1 = 0;

	for (int edge_index = 0; edge_index != edges.size(); edge_index++){
		eptr = &(*(*this).edges[edge_index]);

		//Get the force from the pressure difference across the edge
		v1 = (*eptr).GetOtherVertex(this);
		c1 = &(*(*eptr).cell1);
		v2 = (*c1).GetNextCWVert(this, (*c1).SavedInfoForEvolve_CWVerts);

		//Gets the CW ordering of cells 1 and 2 across eptr
		if ((*v1).id == (*v2).id){ deltaP = (*(*eptr).cell2).SavedInfoForEvolve_Pressure - (*(*eptr).cell1).SavedInfoForEvolve_Pressure; }
		else{ deltaP = (*(*eptr).cell1).SavedInfoForEvolve_Pressure - (*(*eptr).cell2).SavedInfoForEvolve_Pressure; }

		Pmag = -0.5 * deltaP * (*eptr).GetLength();   		//magnitude of the force due to the pressure difference
		angle = (*eptr).GetAngle(this) + PI / 2;	 	//get the angle zXl

		Vect FP_this_edge(zerovec, angle, Pmag);

		FP = FP + FP_this_edge;
	}

	return FP;
}

//! TEMP: calculate net force due to tension on a vertex
Vect Vertex::GetNetForce_Tension(){
	Vect FT(0.0, 0.0);
	for (int edge_index = 0; edge_index != edges.size(); edge_index++){
		FT = FT + (*(edges.at(edge_index))).GetTension_ForceOnVertex(this);
	}

	return FT;
}
 
//if measuring window crosses boundary, defines where vertex positions are outide of PBC wrt original circle
void Vertex::PBCCirclePos(int scenario, double Rposx, double Rposy) {
	//case 1: vertex within original circle - position is same as original with box
	if (scenario == 1) {
		rPBCcross.x = r.x;
		rPBCcross.y = r.y;
	}
	//case 2: vertex needs to cross a single side
	else if (scenario == 2) {
		//determine which quadrant (top, bottom, left, right) circleoverlaps and move vertex over that side
		if (Rposx < 0) { //if this circle overlaps into the left side
			rPBCcross.x = r.x + LX;
			rPBCcross.y = r.y;
		}
		else if (Rposx > LX) { //if overlaps into the right side
			rPBCcross.x = r.x - LX;
			rPBCcross.y = r.y;
		}
		else if (Rposy < 0) { //if overlaps into bottom side
			rPBCcross.x = r.x;
			rPBCcross.y = r.y + LY;
		}
		else if (Rposy > LY) { //if overlaps into top side
			rPBCcross.x = r.x;
			rPBCcross.y = r.y - LY;
		}
		else {
			cout << "ERROR: Circle crossing PBC problem!!" << endl;
		}
	}
	//case 3: if circle overlaps corner
	else if (scenario == 3) {
		if ( (Rposx < 0) && (Rposy < 0) ) { //circle overlaps the top right
			rPBCcross.x = r.x + LX;
			rPBCcross.y = r.y + LY;
		}
		else if ( (Rposx > LX) && (Rposy < 0) ) { //top left
			rPBCcross.x = r.x - LX;
			rPBCcross.y = r.y + LY;
		}
		else if ( (Rposx < 0) && (Rposy > LY) ) { //bottom right
			rPBCcross.x = r.x + LX;
			rPBCcross.y = r.y - LY;
		}
		else if ( (Rposx > LX) && (Rposy > LY) ) { //bottom left
			rPBCcross.x = r.x - LX;
			rPBCcross.y = r.y - LY;
		}
		else {
			cout << "ERROR: Circle crossing PBC problem!!" << endl;
		}
	}

}