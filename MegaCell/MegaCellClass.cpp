/**	\file	MegaCellClass.cpp
*	Implementation of MegaCell class.
*	Contains function definitions associated with MegaCell objects.
*	Depends on base class Cell.
*/

#include "Constants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "MegaCellClass.h"

using namespace std;

// Returns TRUE if the cell crosses the X boundary, FALSE otherwise.
bool MegaCell::GetXFlag(){
	int i;
	int checkflag = 0;
	Edge e;

	for (i = 0; i < edges.size(); i++){			//loop through all of the edges
		e = *edges[i];
		checkflag = checkflag + e.xflag; //keep a running sum of all of the flags
	}

	if (type == outer_cell_type) { return 1; }
	else if (checkflag == 0) { return 0; }
	else if (checkflag>0) { return 1; }
	else cout << "ERROR with GetXFlag" << endl;
}

// Returns TRUE if the cell crosses the Y boundary, FALSE otherwise.
bool MegaCell::GetYFlag(){
	int i;
	int checkflag = 0;
	Edge e;

	for (i = 0; i < edges.size(); i++){			//loop through all of the edges
		e = *edges[i];
		checkflag = checkflag + e.yflag; //keep a running sum of all of the flags
	}

	if (type == outer_cell_type) { return 1; }
	else if (checkflag == 0) { return 0; }
	else if (checkflag>0) { return 1; }
	else cout << "ERROR with GetYFlag" << endl;
}

// Calculates area of a polygon given position of vertices in CW order.
// Area returned may be positive or negative. Negative area means that the vertices were ordered CCW.
double MegaCell::HelperArea(vector<double> &xpos, vector<double> &ypos){
	int i, j;
	double area = 0;
	j = xpos.size() - 1;  // The last vertex is the 'previous' one to the first

	for (i = 0; i < xpos.size(); i++){
		area = area + ((xpos[j] + xpos[i]) * (ypos[j] - ypos[i]));
		j = i;				//j is previous vertex to i
	}

	if (type != outer_cell_type){
		return (area / 2.0);
	}
	else if (type == outer_cell_type && area < 0.0){
		return -(LX * LY - std::fabs(area / 2.0));
	}
	else if (type == outer_cell_type && area > 0.0){
		return LX * LY - std::fabs(area / 2.0);
	}
}