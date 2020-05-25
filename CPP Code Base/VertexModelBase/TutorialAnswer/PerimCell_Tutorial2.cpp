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
#include <typeinfo>
#include <cassert>
#include <numeric>
#include <stddef.h>
#include "Constants.h"

#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "PerimTissue.h"
#include "PerimCell.h"
#include "PerimEdge.h"

/*
//Gets the pressure of the cell as K*(1-A/A0)
double SpringCell::GetPressure(double Area){ return K*(1 - Area / A0); }

//Sets the isCompleete flag
bool SpringCell::SetIsCompleeteFlag(){
	int edgenumber;
	Edge* e;
	Edge* enext = NULL;
	Vertex* v;
	Vertex* vnext;
	int firstv;
	int index = 0;

	edgenumber = (*this).edges.size();
	//If the cell has two or less edges then is is definatly not compleete 
	if (edgenumber < 3){
		(*this).isCompleete = false;
		return false; //set the is complete flag to false and exit
	}

	//If it has some edges
	if ((*this).edges[0] == NULL){
		return false;
	}
	e = (*this).edges[0]; //get the first edge
	v = (*e).vertex1;
	firstv = (*v).id; //get the first vertex and record its id

	//If it has a bad vertex 
	if ((*v).edges.size() > 4){
		(*this).isCompleete = false;
		return false; //set the is complete flag to false and exit
	}

	//Keep walking around the cell untill you either get back to the first vertex or run out of edges
	while (true){
		vnext = (*v).GetOtherVertex(e); //Get the next vertex
		if ((*vnext).id == firstv){		//if it is the same as the starting vertex the cell is closed and we are done
			(*this).isCompleete = true;
			return true;
		}
		//if not keep going
		if (GetNextEdge(e, vnext, &index)){ //if there is a next edge reset for the next run
			enext = (*vnext).edges[index];
			v = vnext;
			e = enext;
		}
		else{ //if there is no next edge then stop
			(*this).isCompleete = false;
			return false;
		}
	}
}

bool SpringCell::SetIsBoundaryFlag(){
	bool flag = 0;
	int i;
	SpringEdge* e;
	SpringCell* c;
	vector<Cell*> cells;

	//Incompleete cells are on the boundary 
	if (!(*this).isCompleete){
		isBoundary = true; 
		return true; 
	}

	// Get this cells edges
	for (i = 0; i < (*this).edges.size(); i++){
		e = dynamic_cast<SpringEdge*>((*this).edges[i]);
		flag = (*e).isBoundary() || flag;         //if any edge is on the boundary then the cell is near the boundary
	}

	isBoundary = flag;
	return flag;
}

//sets the isCentral flag for both cells and vertices 
bool SpringCell::SetIsCentralFlag(){
	bool flag = 0;
	Vect center; 
	vector<Vertex*> verts; 

	if (!isCompleete){
		flag=false; 
	}
	else{
		center = (*this).GetCenter(); 
		if ((center.x > 0.25* LX) && (center.x < 0.75*LX)){
			flag = true;
			verts = (*this).GetVertices(); 
			for (int i = 0; i < verts.size(); i++){
				SpringVertex* v = dynamic_cast<SpringVertex*>(&(*(verts[i]))); 
				(*v).isCentral = true; 
			}
		}
		else{
			flag = false; 
		}
	}
	
	isCentral = flag;
	return flag;
}
//takes and edge and a vertex on that edge and returns true if that vertex has another edge in this cell
bool SpringCell::GetNextEdge(Edge* e, Vertex* v, int* result){
	bool flag;
	Edge* teste;
	int cid = (*this).id;

	flag = false;
	//for every edge attached to this vertex (there must be at least one because you gave it e)
	for (int i = 0; i < (*v).edges.size(); i++){
		teste = (*v).edges[i];
		if ((*e).id != (*teste).id){		//if the edge is not the same one you put in 
			if ((*teste).cell1 != NULL && (*teste).cell2 != NULL){
				if (cid == (*(*teste).cell1).id || cid == (*(*teste).cell2).id){  //and one of its cells matches this cell
					flag = true;
					(*result) = i;
				}
			}
		}
	}

	return flag;
}

//Takes a vertex and outputs the pressure force from this cell on the vertex ASSUMES CWORDERING AND PRESSURE KNOWN
Vect SpringCell::GetPressureForceOnVertex(Vertex* v){
	vector< Edge*> vedges;
	Edge *e1, *e2;
	int i;
	Vect l1, l2, la, lb;
	Vect o(0, 0);
	Vect zxl1minuml2;
	Vect l1minusl2;
	Vect force;

	//find the two edges that are on this vertex and this cell
	for (i = 0; i < edges.size(); i++){  //for each edge in this cell 
		if ((*(*edges[i]).vertex1).id == (*v).id || (*(*edges[i]).vertex2).id == (*v).id){ //if the edge contains this vertex
			vedges.push_back(edges[i]);
		}
	}
	e1 = vedges[0];
	e2 = vedges[1];

	//Get the lengths of the edges and their direction
	la = Vect(o, (*e1).GetAngle(v), (*e1).GetLength());
	lb = Vect(o, (*e2).GetAngle(v), (*e2).GetLength());

	//determine the order to subtract the edges
	if (la.SignCrossProduct(lb) > 0){
		l1 = la;
		l2 = lb;
	}
	else{
		l1 = lb;
		l2 = la;
	}

	//Get the vector part
	l1minusl2 = l1 - l2;
	zxl1minuml2 = l1minusl2.Rotate(PI / 2);

	//get the force
	force = zxl1minuml2 * SavedInfoForEvolve_Pressure / 2 * -1;
	return force;
}

//Returns 1 if the cell is convex 2 if notconvex and 3 if it is complex (ie crosses its self) 
int SpringCell::GetPolygonComplexity(){
	vector<vector<double>> xypos;
	bool convexflag;
	bool complexflag;
	bool cross;
	int i, j;
	Edge *e1, *e2;

	//if theis cell is the big outer cell return 3
	if (!(*this).isCompleete) { return 3; }

	// first check if the cell overlaps its self 
	complexflag = 0;
	for (i = 0; i < edges.size() - 1; i++){
		for (j = i + 1; j < edges.size(); j++){ //for every pair of edges 
			e1 = edges[i];
			e2 = edges[j];
			if (!(*e1).ShareVertex(e2) && EdgesIntersect(e1, e2)){      //only l0ok at edges that dont share a vertex
				return 3;
			}
		}
	}

	//If the cell is not complex then check convexity
	xypos = (*this).HelperOrderVerticesPos(SavedInfoForEvolve_CWVerts);      //first get the x,y positions (using the fact that they are already in CW order) 
	convexflag = CheckConvex(xypos);										//Determine if it is convex
	if (convexflag) return 1;
	else return 2;
}

//Returns the x or y extent of the cell (ie the total length in the x or y direction of the cell)
double SpringCell::GetXExtent(){
	vector<Vertex*> verts;					//vertices in the cell
	vector<double> xpos;					//x positions of the vertices
	int i;
	double result;

	//If the cell does not cross the X-boundary (Our cells can nexer cross the boundary!)
	verts = (*this).GetVertices();		//get the vertices
	//get the x pos of each vertex
	for (i = 0; i < verts.size(); i++){
		xpos.push_back((*verts[i]).r.x);
	}
	//find the min and max
	auto min = min_element(xpos.begin(), xpos.end());
	auto max = max_element(xpos.begin(), xpos.end());
	result = (*max) - (*min);
	return result;
}
double SpringCell::GetYExtent(){
	vector<Vertex*> verts;					//vertices in the cell
	vector<double> ypos;					//x positions of the vertices
	int i;
	double result;

	//If the cell does not cross the X-boundary (Our cells can nexer cross the boundary!)
	verts = (*this).GetVertices();		//get the vertices
	//get the y pos of each vertex
	for (i = 0; i < verts.size(); i++){
		ypos.push_back((*verts[i]).r.y);
	}
	//find the min and max
	auto min = min_element(ypos.begin(), ypos.end());
	auto max = max_element(ypos.begin(), ypos.end());
	result = (*max) - (*min);
	return result;
}

//returns true if cell is on boundary or one of its neighbors is
bool SpringCell::NearBoundary(){
	bool flag = 0;
	int i;
	SpringCell* c;
	vector<Cell*> cells;

	//if this cell is on the boundary return true
	if ( isBoundary )
	{
		return true; 
	}
	else{
		// Get this cells edges
		cells = (*this).GetNeighboringCells();
		for (i = 0; i < cells.size(); i++){
			c = dynamic_cast<SpringCell*>(cells[i]);
			flag = (*c).isBoundary || flag;         //if any edge is on the boundary then the cell is near the boundary
		}
	}

	return flag;
}

double SpringCell::GetMeanTension(){
	double sum = 0; 
	double t; 

	for (int i = 0; i < (*this).edges.size(); i++){
		t = abs((*(*this).edges[i]).GetTension_Magnitude()); 
		sum = sum + t; 
	}

	return sum / edges.size(); 
	}*/