/**	\file	MegaCellClass.h
*	Header of MegaCell class.
*	Contains data and function declarations associated with MegaCell objects.
*	Depends on base class Cell.
*/

// include guard
#ifndef __MEGACELLCLASS_H_INCLUDED__
#define __MEGACELLCLASS_H_INCLUDED__

//Forward declared dependences
class Cell;
class Edge;
class Vertex;
class Vect;

#include "Constants.h"
#include "MegaCellConstants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h" 

//! Class for cell objects.
class MegaCell : public Cell{
public:

	//ADDITIONAL FUNCTIONS
	//CONSTRUCTORS AND DESTRUCTORS
	MegaCell(){};
	virtual ~MegaCell(){};

	//! Returns TRUE if the cell crosses the X boundary, FALSE otherwise.
	//! Outer MegaCell needs a flag
	virtual bool GetXFlag();

	//! Returns TRUE if the cell crosses the Y boundary, FALSE otherwise.
	//! Outer MegaCell needs a flag
	virtual bool GetYFlag();

	//! Calculates area of a polygon given position of vertices in CW order.
	//! Area returned may be positive or negative. Negative area means that the vertices were ordered CCW.
	//! Outer MegaCell's area must be addressed differently.
	virtual double HelperArea(vector<double> &xpos, vector<double> &ypos);

	//! Assign constant inner pressure to inner mega cell, assign constant outer pressure to outer mega cell, else assign K(A0-Area).
	virtual double GetPressure(double Area){
		//If the packing is not a voronoi packing, then impose inner and outer pressure.
		if (is_voronoi_packing == 0){
			if (type == center_cell_type)
				return inner_pressure;
			else if (type == outer_cell_type)
				return outer_pressure;
			else
				return K * (A0 - Area);
		}
		//If the packing is a voronoi packing, then impose only outer pressure.
		else{
			if (type == outer_cell_type)
				return outer_pressure;
			else
				return K * (A0 - Area);
		}
	}
};


#endif