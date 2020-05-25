/**	\file	PCPMegaCellClass.cpp
*	Implementation of PCPMegaCell class.
*	Contains function definitions associated with PCPMegaCell objects.
*	Depends on base class Cell.
*/

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
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "PCPMegaCellClass.h"
#include "PCPEdgeClass.h"
#include "PCPMegaTissueClass.h"
#include "PCPConstants.h"

using namespace std;

// Returns TRUE if the cell crosses the X boundary, FALSE otherwise.
bool PCPMegaCell::GetXFlag(){
	int i;
	int checkflag = 0;
	Edge e;

	for (i = 0; i < edges.size(); i++){			//loop through all of the edges
		e = *edges[i];
		checkflag = checkflag + e.xflag; //keep a running sum of all of the flags
	}

	//TEMP FOR PCP WITH MEGA CELL
	if (type == outer_cell_type) { return 1; }
	//
	else if (checkflag == 0) { return 0; }
	else if (checkflag>0) { return 1; }
	else cout << "ERROR with GetXFlag" << endl;
}

// Returns TRUE if the cell crosses the Y boundary, FALSE otherwise.
bool PCPMegaCell::GetYFlag(){
	int i;
	int checkflag = 0;
	Edge e;

	for (i = 0; i < edges.size(); i++){			//loop through all of the edges
		e = *edges[i];
		checkflag = checkflag + e.yflag; //keep a running sum of all of the flags
	}

	//TEMP FOR PCP WITH MEGA CELL
	if (type == outer_cell_type) { return 1; }
	//
	else if (checkflag == 0) { return 0; }
	else if (checkflag>0) { return 1; }
	else cout << "ERROR with GetYFlag" << endl;
}

//! When a cells gets PCP, start clocks in neighboring cells (correlating the clocks by correlation_coefficient).
void PCPMegaCell::Helper_AddBorderingCellsPoisson(){
	//Find the cells neighboring the cell with the pointer cell_with_PCP.
	vector<Cell*> bordering_cells = (*this).GetNeighboringCells();
	PCPMegaCell* pcp_cell_ptr;

	//Store the number of bordering cells.
	int number_of_bordering_cells = bordering_cells.size();

	//add pcp to cell_with_PCP
	(*this).HelperEvolve_InitializePolarityProteins();

	//GAUSSIAN WITH NO CORRELATIONS:
	for (int bordering_cell_index = 0; bordering_cell_index != number_of_bordering_cells; bordering_cell_index++){
		pcp_cell_ptr = dynamic_cast<PCPMegaCell*>(bordering_cells.at(bordering_cell_index));
		bool already_ticking = ((*pcp_cell_ptr).type == ticking_type);
		bool is_inner_mega_cell = ((*pcp_cell_ptr).type == center_cell_type);
		bool is_outer_mega_cell = ((*pcp_cell_ptr).type == outer_cell_type);
		bool bordering_cell_with_pcp = (*pcp_cell_ptr).is_Cell_with_PCP();

		//If the cell does not have PCP, start the clock ticking.
		if (!already_ticking && !is_inner_mega_cell && !is_outer_mega_cell && !bordering_cell_with_pcp){
			(*pcp_cell_ptr).cell_generation = (*this).cell_generation + 1;
			(*pcp_cell_ptr).Helper_GenerateRandomTime();
		}
	}
}

//! This function is used to calculate the distance from the origin of PCP (centroid of central megacell).
double PCPMegaCell::calculate_Distance_from_PCP_Origin(){
	//Store the distance for this cell from the origin of the mega cell.
	//use origin_mega_cell (computed in force function).
	Vect origin_this_cell = (*this).GetCenter();
	double distance_from_mega_cell = (origin_mega_cell - origin_this_cell).Magnitude();

	return distance_from_mega_cell;
}

//void PCPMegaCell::Find_Five_Nearest_Neighbors(){
//	//HELPER FIND CELL ANGLE: +Pi/2, -Pi/2, +Pi (if angle specifies direction outward).
//
//	vector<double> angle_vect;
//	angle_vect.push_back(PI);
//	angle_vect.push_back(PI / 2.0);
//	angle_vect.push_back(3.0 * PI / 2.0);
//
//	//Find three neighboring cells.
//	vector<Cell*> neighboring_cells = (*this).Return_Cell_in_Same_Column(angle_vect);
//
//	//Find cell in previous column.
//	Cell* cell_in_previous_column = neighboring_cells.at(0);
//
//	vector<double> previous_column_angle_vect;
//	previous_column_angle_vect.push_back(PI / 2.0);
//	previous_column_angle_vect.push_back(3.0 * PI / 2.0);
//
//	//Find neighboring cells in previous column.
//	vector<Cell*> neighboring_cells_previous_column = (*this).Return_Cell_in_Same_Column(previous_column_angle_vect);
//
//	//Store vector of types, and use more complex type specification.
//}

// Calculates area of a polygon given position of vertices in CW order.
// Area returned may be positive or negative. Negative area means that the vertices were ordered CCW.
double PCPMegaCell::HelperArea(vector<double> &xpos, vector<double> &ypos){
	int i, j;
	double area = 0;
	j = xpos.size() - 1;	//The last vertex is the 'previous' one to the first

	for (i = 0; i < xpos.size(); i++){
		area = area + ((xpos[j] + xpos[i]) * (ypos[j] - ypos[i]));
		j = i;				//j is previous vertex to i
	}

	//TEMP FOR PCP WITH MEGA CELL
	if (type != outer_cell_type){
		return (area / 2.0);
	}
	else if (type == outer_cell_type && area < 0.0){
		return -(LX * LY - std::fabs(area / 2.0));
	}
	else if (type == outer_cell_type && area > 0.0){
		return LX * LY - std::fabs(area / 2.0);
	}
	//
}