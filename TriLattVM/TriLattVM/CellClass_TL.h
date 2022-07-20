/**	\file	CellClass.h
*	Header of Cell class.
*	Contains data and function declarations associated with Cell objects.
*/

// include guard
#ifndef __CELLCLASS_H_INCLUDED__
#define __CELLCLASS_H_INCLUDED__

//Forward declared dependences
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
#include "Constants_TL.h"

using namespace std;


//! Class for cell objects.
class Cell{
public:
	//DATA
		unsigned long int id;							//!< ID of a cell. Used to determine if two cells are the same.
		int type;										//!< Integer specifying type of cell. 
		double A0;										//!< Preferred area of the cell.
		vector<Edge*> edges;							//!< Pointers to edges of the cell.
		vector<Vertex*> SavedInfoForEvolve_CWVerts;		//!< A place to save a list of the vertices in the cell in CW order used in evolve only!
		double SavedInfoForEvolve_Pressure;				//!< A place to save a the pressure in the cell during evolve
		vector<int> IOEdgeID;							//!< list of the ids of the edges for input and output purposes (this is nessacary for setting up the pointers when the initial data is read in from file). Used only in reading input data.

	//FUNCTIONS
		//CONSTRUCTORS AND DESTRUCTORS
		virtual ~Cell(){};						//destructor: does not need to be called by you unless you created the object using the word "new" otherwise the compiler will take care of this for you
		Cell(){};

		//! Constructs a cell with ID 'idnumb' and edge IDs 'edglist'.
		Cell(int idnumb, vector<int> edglist){id=idnumb; IOEdgeID=edglist;}

		//! Constructs a cell of type 'typeinput' with ID 'idnumb', edge IDs 'edglist', and preferred area 'A0input'.
		Cell(int idnumb, vector<int> edglist, int typeinput, double A0input){ id = idnumb; IOEdgeID = edglist; type = typeinput; A0 = A0input; }

		//! Constructs a cell of type 'typeinput' with ID 'idnumb' and pointers to edges 'edg'.
		Cell(int idnum, int typeinput, vector<Edge*> edg){
			id=idnum;
			type=typeinput;
			edges = edg;
		}


		//! Tests if two cells are the same.
		//! Compares ID of the two cells and returns 'true' if they are the same, 'false' otherwise.
		bool operator==(const Cell &) const;

		//! Tests if two cells are NOT the same.
		//! Compares ID of the two cells and returns 'true' if they are NOT the same, 'false' otherwise.
		bool operator!=(const Cell &right) const{
			return !(*this == right);
		}

		//! Returns pointers to all vertices of the cell.
		vector<Vertex*> GetVertices();

		//! Returns the number of vertices of the cell.
		int GetVertexNumber();

		//! Removes edge 'eptr' from the cell.
		void RemoveEdge( Edge* eptr ); 

		//! Returns TRUE if the cell crosses the X boundary, FALSE otherwise.
		virtual bool GetXFlag();

		//! Returns TRUE if the cell crosses the Y boundary, FALSE otherwise.
		virtual bool GetYFlag();

		//! Returns TRUE if the cell crosses any boundary boundary, FALSE otherwise.
		bool GetFlag();

		//! Returns center of the cell, respecting PBC. For cells that crosses bondaries, the function returns the center position WITHIN the box.
		//! \note This function calls GetGeoCenter() and HelperVerticesPos().
		Vect GetCenter();

		//! Updates 'rnew' of the vertices to the center of the cell (noting that the center may be outside the box for PBC), to be used for apoptosis. This calls the function GetGeoCenter().
		void MoveVertsToCenter();

		//! Returns a vector of pointers to all vertices in the cell in clockwise (CW) order.
		vector<Vertex*> GetClockWiseVertices();

		//! GetArea returns the area of the cell.
		//! Actual implementation of algorithm that calculates the area is in function HelperArea(vector<double> &xpos, vector<double> &ypos).
		double GetArea();

		//! Faster version if you already know the cw vertices (used in evolve)
		double GetArea(vector<Vertex*> cwverts);

		//! Finds the perimeter of the cell
		double GetPerimeter();

		//! Same as GetPressure() but with cell area as input.
		//! \attention This function contains the actual implementation of GetPressure. All other GetPressure functions should call this function to prevent multiple implementations. Furthermore any modification/extension to the algorithm should only be made here.
		virtual double GetPressure(double Area) { return 1; }

		//! Returns a vector of pointers to cells in contact with the given cell (used for building contact matrix)
		vector<Cell*> GetNeighboringCells();

		//! Returns a vector of contact length between neighboring cells (same ordering as GetNeighboringCells: used for building contact matrix)
		vector<double> GetNeighboringCells_ContactLength();

		//! Returns common edge between this cell and other cell.
		Edge* Find_Common_Edge(Cell* other_cell);

	//Helper functions: these functions should not be called in the main program but exist to support the primary functions above
	//some of the primary functions reuse the same code and so this code has been split into a seperate function that multiple primary functions can call
		
		//! Calculates area of a polygon given position of vertices in CW order.
		//! Area returned may be positive or negative. Negative area means that the vertices were ordered CCW.
		virtual double HelperArea(vector<double> &xpos, vector<double> &ypos);
		
		//! Returns a 2D vector of x and y position of vertices, ordered either CW or CCW. Respects PBC.
		vector<vector<double>> HelperOrderVerticesPos();

		// Faster calculation of HelperOrderVerticesPos() due to precalculated list of cwvertices 
		vector<vector<double>> HelperOrderVerticesPos(vector<Vertex*> cwverts);

		// Faster calculation of HelperOrderVerticesPos() due to precalculated list of cwvertices: starts building vertex position list with the position of *first_vertex.
		// \note This is important because with PBC, it is often useful to build a cwverts pos list that respects PBC and is on one side of the boundary (the side specified by *first_vertex).
		vector<vector<double>> HelperOrderVerticesPos(vector<Vertex*> cwverts, Vertex* first_vertex);

		//! Returns a list of pointers to the vertices in either CW or CCW order.
		vector<Vertex*> HelperOrderVertices();

		//! Returns flags of the edges between vertices given in OrderedVerts (i.e. xFlag[i] is the flag of edge between vertices OrderedVerts[i] and OrderedVerts[i+1]).
		//! Use HelperOrderVertices() or GetClockWiseVertices() to get OrderedVerts.
		//! xFlag[i] = 0 if on the same side of OrderedVerts[0], -1 if crossing the boundary from low to high, and +1 if crossing boundary from high to low.
		//! Hence, the position of vertices without PBC is given by x1[i] = x[i] + xFlag[i] * LX. This is useful in GetArea() calculation.
		vector<vector<int>> HelperOrderFlags(vector<Vertex*> OrderedVerts);
		
		// This returns order vertex flags (as HelperOrderFlags(vector<Vertex*> OrderedVerts) does), but starts the list with *first_vertex.
		// \note This is important because with PBC, it is often useful to build a list of vertex positions that respects PBC and is on one side of the boundary (the side specified by *first_vertex).
		vector<vector<int>> HelperOrderFlags(vector<Vertex*> OrderedVerts, Vertex* first_vertex);

		//! Removes an edge from the list of edges.
		void HelperRemoveEdge(Edge* eptr);

		//! Given a list of CW vertices in the cell it returns the next vertex after v0
		Vertex* GetNextCWVert(Vertex* v0, vector<Vertex*> cwvert);

		//Returns the length of the line at pos x (or y) that goes through the cell
		double GetLengthHorzCut(double ycut);
		double GetLengthVertCut(double xcut);

		//Returns the length of the line at pos x (or y) that goes through the cell
		//Faster version that starts with the positions of the verticies (where the positions have been adjusted so that the cell may lie outside the PBC box) 
		double GetLengthHorzCut(double ycut, vector<vector<double>> vertposlist);
		double GetLengthVertCut(double xcut, vector<vector<double>> vertposlist);
};

//! This is a general static function that calculates the geometrical center of a polygon, given some x positions and y positions of vertices.
//! Please note that that geometrical center may be outside of the boundary if some of these specified input position values are outside of the boundary. GetCenter() will perform additional calculations to ensure the center is within the box.
//! \note Vertex positions must be listed in order (either CW or CCW).
Vect GetGeoCenter(vector<double> xpos, vector<double> ypos);

//! Tests if the cell is convex.
//! Returns TRUE if the cell is convex, FALSE otherwise. Takes as input list of ordered vertex positions, and uses cross product along the polygon to determine convexity.
//! \note If you already have orderedvertices stored, use HelperOrderVerticesPos(const double y[], unordered_map<long, int> verts_id2index, vector<Vertex*> orderedVerts); otherwise use HelperVerticesPos().
bool CheckConvex(vector<vector<double>> vertposlist);

//! Using a list of vertex positions (vertposlist), define a ray from the geocenter to the right. 
//! Return the indices of the edge which the ray intersects. If the edge is between (vertex i, vertex i+1) in vertposlist, return i. Include all intersections.
// \note This is used for defining front propogation for PCP as well as for point in polygon.
vector<int> which_Edge_Intersection(Vect position, vector<vector<double>> vertposlist, double angle);

//! Return a vector of point where the ray defined by (cell_cut_vector_1 , cell_cut_vector_2) intersects the polygon defined by the vertex positions xypos.
// \note This is used within CalculateTissueStressX and CalculateTissueStressY to define the intersection between cells and the cutline.
vector<Vect> store_Intersections(Vect cell_cut_vector_1, Vect cell_cut_vector_2, vector<vector<double>> xypos);

//! Check if the point position is in the cell specified by vertposlist. (vertposlist list is given for speeding up calculation)
//! \note This uses the ray-casting algorithm: If a ray going outward from the point in an arbitrary direction intersects the polygon an even number of times, the point is outside the polygon. 
//! \note If a ray going outward from the point in an arbitrary direction intersects the polygon an odd number of times, the point is inside the polygon. 
bool is_Point_in_Polygon(Vect position, vector<vector<double>> vertposlist);

#endif