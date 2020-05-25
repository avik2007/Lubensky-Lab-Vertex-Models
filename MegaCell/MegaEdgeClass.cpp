/**	\file	MegaEdgeClass.cpp
*	Implementation of MegaEdge class.
*	Contains function definitions associated with MegaEdge objects.
*	Depends on base class Edge.
*/

#include "Constants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "MegaCellClass.h"
#include "MegaEdgeClass.h"

using namespace std;

// \note MUST DEFINE GET TENSION MAGNITUDE WITHIN SUBCLASS IF YOU REQUIRE FANCY MAGNITUDE->CALLED INDEPENDENTLY WITHIN CALCULATETISSUESTRESS
double MegaEdge::GetTension_Magnitude(){
	//TEMP FOR MEGA CELL (MAKE INNER BOUNDARY SMOOTH)
	if ((*cell1).type == center_cell_type || (*cell2).type == center_cell_type){
		return factor_increased_tension_at_inner_boundary * GAMMA;
	}
	else{
		return GAMMA;
	}
	//
}