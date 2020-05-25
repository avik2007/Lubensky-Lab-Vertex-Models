/**	\file	PCPCellClass.cpp
*	Implementation of PCPCell class.
*	Contains function definitions associated with PCPCell objects.
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
#include <cassert>					//for error checking/assertion
#include <random>
#include "Constants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "PCPCellClass.h"
#include "PCPEdgeClass.h"

using namespace std;

//! Checks whether a cell has PCP.
bool PCPCell::is_Cell_with_PCP(){
	bool differentiated_cell = false;
	//if the cell is of a type that has PCP.
	if (type == pcp_cell_type || type == checkerboard_type1 || type == checkerboard_type2){
		differentiated_cell = true;
	}
	return differentiated_cell;
}

//! This function is used to calculate the distance from the origin of PCP.
double PCPCell::calculate_Distance_from_PCP_Origin(){
	//This assumes that propagation is to the right.

	//Store the distance for this cell from the origin of propagation (the original cutline).
	Vect origin_this_cell = (*this).GetCenter();
	double origin_this_cell_x = origin_this_cell.x;
	double distance_this_cell;
	if (origin_this_cell_x >= cutline_pcp_initialization)
		distance_this_cell = origin_this_cell_x - cutline_pcp_initialization;
	else
		distance_this_cell = (LX - cutline_pcp_initialization) + origin_this_cell_x;

	return distance_this_cell;
}

//! This function checks if the cell is in contact with a cell with PCP.
bool PCPCell::CheckNeighboringCell_forPCP(){
	//Store pointers to neighboring cells.
	vector<Cell*> neighboring_cells = (*this).GetNeighboringCells();
	int number_neighboring_cells = neighboring_cells.size();
	//Store bool (is the cell connected to a cell with PCP?)
	bool connected_to_pcp_cell = false;

	for (int cell_index = 0; cell_index != number_neighboring_cells; cell_index++){
		PCPCell* pcpcellptr = dynamic_cast<PCPCell*>(neighboring_cells.at(cell_index));
		//If neighboring cell has PCP,
		if ((*pcpcellptr).is_Cell_with_PCP())
			connected_to_pcp_cell = true;
	}

	return connected_to_pcp_cell;
}

//! Find minimum cell generation neighboring the cell, update cell generation to min generation + 1.
void PCPCell::Update_Cell_Generation(){
	vector<Cell*> neighboring_cells = (*this).GetNeighboringCells();
	int number_neighboring_cells = neighboring_cells.size();
	vector<int> cell_generation_vector;
	int min_neighboring_generation;
	int this_cell_generation;

	//Loop through neighboring cells: check if at least one of the neighboring cells is of the previous cell generation.
	for (int cell_index = 0; cell_index != number_neighboring_cells; cell_index++){
		PCPCell* neighboring_pcp_cell_ptr = dynamic_cast<PCPCell*>(neighboring_cells.at(cell_index));
		this_cell_generation = (*neighboring_pcp_cell_ptr).cell_generation;
		if (this_cell_generation != -1){
			cell_generation_vector.push_back(this_cell_generation);
		}
	}

	if ((*this).type == ticking_type){
		min_neighboring_generation = *min_element(cell_generation_vector.begin(), cell_generation_vector.end());
		if (min_neighboring_generation != (*this).cell_generation - 1){
			(*this).cell_generation = min_neighboring_generation + 1;
			(*this).Helper_GenerateRandomTime();
		}
	}
}

//! Returns a list of ordered edges (given a list of ordered vertices)
vector<PCPEdge*> PCPCell::FindOrderedEdges(){
	vector<PCPEdge*> orderedEdges;
	int number_of_verts = SavedInfoForEvolve_CWVerts.size();
	orderedEdges.resize(number_of_verts);

	//Loop through vertices.
	for (int vertex_index = 0; vertex_index != number_of_verts; vertex_index++){
		//Find the next vertex index (with mod->loop to beginning).
		int next_vertex_index = (vertex_index + 1) % number_of_verts;
		//Get pointer for edge between the two vertices.
		orderedEdges.at(vertex_index) = dynamic_cast<PCPEdge*>((*SavedInfoForEvolve_CWVerts.at(vertex_index)).GetEdgeBetween(SavedInfoForEvolve_CWVerts.at(next_vertex_index)));
	}
	return orderedEdges;
}

//! Return next edge CW from (*eptr)
PCPEdge* PCPCell::GetNextCWEdge(PCPEdge* eptr){
	Vertex* this_vertex1 = (*eptr).vertex1;
	Vertex* this_vertex2 = (*eptr).vertex2;

	Vertex* vertex1_nextptr = (*this).GetNextCWVert(this_vertex1, (*this).SavedInfoForEvolve_CWVerts);
	Vertex* vertex2_nextptr = (*this).GetNextCWVert(this_vertex2, (*this).SavedInfoForEvolve_CWVerts);

	Edge* next_edge_ptr = NULL;
	PCPEdge* pcp_next_edge_ptr = NULL;

	if ((*vertex2_nextptr).id == (*this_vertex1).id){
		next_edge_ptr = (*vertex1_nextptr).GetEdgeBetween(this_vertex1);
	}
	else if ((*vertex1_nextptr).id == (*this_vertex2).id){
		next_edge_ptr = (*vertex2_nextptr).GetEdgeBetween(this_vertex2);
	}

	pcp_next_edge_ptr = dynamic_cast<PCPEdge*>(next_edge_ptr);

	return pcp_next_edge_ptr;
}

//! NOTE: This is only useful for simulations on a hexagonal field of cells.
PCPEdge* PCPCell::FindUpperHorizontalEdge(){
	vector<Vertex*> this_vertex_vector = (*this).SavedInfoForEvolve_CWVerts;
	double err_tol = 1.0e-06;
	PCPEdge* chosen_edge = NULL;
	Edge* between_edge = NULL;

	for (int vertex_index = 0; vertex_index != edges.size(); vertex_index++){
		int next_vertex_index = (vertex_index + 1) % edges.size();

		between_edge = (*this_vertex_vector[vertex_index]).GetEdgeBetween(this_vertex_vector[next_vertex_index]);

		double this_angle = (*between_edge).GetAngle(this_vertex_vector[next_vertex_index]);

		if (fabs(this_angle) <= err_tol){
			PCPEdge* pcp_between_edge = dynamic_cast<PCPEdge*>(between_edge);
			chosen_edge = pcp_between_edge;
		}
	}

	return chosen_edge;
}

//! Order edges from list of ordered vertices. NOTE: This is only used for simulations on a hexagonal field of cells.
vector<PCPEdge*> PCPCell::OrderEdges(){
	vector<Vertex*> this_vertex_vector = (*this).SavedInfoForEvolve_CWVerts;
	vector<PCPEdge*> ordered_edge_vect;
	PCPEdge* between_edge = NULL;

	PCPEdge* upper_edge = (*this).FindUpperHorizontalEdge();
	between_edge = upper_edge;

	do{
		ordered_edge_vect.push_back(between_edge);
		between_edge = (*this).GetNextCWEdge(between_edge);
	} while ((*between_edge).id != (*upper_edge).id);

	return ordered_edge_vect;
}

//HELPER FUNCTIONS

//! Generate a random time for the given cell.
void PCPCell::Helper_GenerateRandomTime(){
	double mean_value = ((*this).cell_generation) * mean_pcpwaitingtime;

	normal_distribution<double> distribution(mean_value, standard_dev_pcpwaitingtime);
	double random_time = distribution(generator);

	(*this).type = ticking_type;
	(*this).time_of_pcp = random_time;
}

//! When a cells gets PCP, start clocks in neighboring cells.
void PCPCell::Helper_AddBorderingCellsPoisson(){
	//Find the cells neighboring the cell with the pointer cell_with_PCP.
	vector<Cell*> bordering_cells = (*this).GetNeighboringCells();
	PCPCell* pcp_cell_ptr;

	//Store the number of bordering cells.
	int number_of_bordering_cells = bordering_cells.size();

	//add pcp to cell_with_PCP
	(*this).HelperEvolve_InitializePolarityProteins();

	//GAUSSIAN WITH NO CORRELATIONS:
	for (int bordering_cell_index = 0; bordering_cell_index != number_of_bordering_cells; bordering_cell_index++){
		pcp_cell_ptr = dynamic_cast<PCPCell*>(bordering_cells.at(bordering_cell_index));
		bool already_ticking = ((*pcp_cell_ptr).type == ticking_type);
		bool intersect_cell = ((*pcp_cell_ptr).type == intersect_type);
		bool bordering_cell_with_pcp = (*pcp_cell_ptr).is_Cell_with_PCP();

		//If the cell does not have PCP, start the clock ticking.
		if (!already_ticking && !intersect_cell && !bordering_cell_with_pcp){
			(*pcp_cell_ptr).cell_generation = (*this).cell_generation + 1;
			(*pcp_cell_ptr).Helper_GenerateRandomTime();
		}
	}
}

//! This function is used to create an initially randomized distribution of PCP proteins in the cell pointed to be cellptr. Called by Evolve.
void PCPCell::HelperEvolve_InitializePolarityProteins(){
	int number_of_edges = edges.size();				//How many edges does the current cells have?
	vector<double> alpha_vector;					//This will be used to bias larger total protein amounts for larger edges.
	alpha_vector.resize(number_of_edges);
	vector<double> length_vector;					//Store the length of each edge (in order to later define the concentration)
	length_vector.resize(number_of_edges);
	vector<double> proteinA_amounts;				//Vector for randomly generated amounts of protein A for each edge.
	proteinA_amounts.resize(number_of_edges);
	vector<double> proteinB_amounts;				//Vector for randomly generated amounts of protein B for each edge.
	proteinB_amounts.resize(number_of_edges);
	int edge_index = 0;
	vector<int> which_cell;							//Used to store an indicator of which element of planarpolarity to access.
	which_cell.resize(number_of_edges);

	//Used to convert pointers from edges to pointers to PCPEdges
	vector<PCPEdge*> edge_list_this_cell;
	edge_list_this_cell.resize(number_of_edges);
	for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
		edge_list_this_cell.at(edge_index) = dynamic_cast<PCPEdge*>(edges.at(edge_index));
	}
	
	//Set the type of the newly cell to be 1 (which means that it is differentiated).
	//Output an error message if the cell is already differentiated.
	if ((*this).is_Cell_with_PCP()){
		cout << "ERROR: You are reinitializing proteins in an already differentiated cell. Tell HN." << endl;
	}
	type = pcp_cell_type;

	for (edge_index = 0; edge_index != number_of_edges; edge_index++){
		length_vector[edge_index] = (*edges[edge_index]).GetLength();
		double angle_this_edge = (*edge_list_this_cell.at(edge_index)).GetAngle((*edge_list_this_cell.at(edge_index)).vertex1);
		alpha_vector[edge_index] = dirichlet_noise_param * length_vector[edge_index];
		if ((*(*edges[edge_index]).cell1).id == (*this).id){
			which_cell[edge_index] = 0;
		}
		else {
			which_cell[edge_index] = 1;
		}
	}
	//do random assignment of protein A amounts
	gsl_ran_dirichlet(polarity_random_generator, number_of_edges, &alpha_vector[0], &proteinA_amounts[0]);
	//do random assignments of protein B amounts
	gsl_ran_dirichlet(polarity_random_generator, number_of_edges, &alpha_vector[0], &proteinB_amounts[0]);
	for (edge_index = 0; edge_index != edges.size(); edge_index++){
		(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = proteinA_amounts[edge_index] * (1 + delta_protein_amount);
		(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = proteinB_amounts[edge_index] * (1 - delta_protein_amount);
	}
}

//! This function is used to create a stereotyped distribution of protein for either edge polarized or vertex polarized states. Called by ReadInputData.
void PCPCell::HelperReadInputData_InitializePolarityProteins(){
	int number_of_edges = edges.size();				//How many edges does the current cells have?
	vector<double> alpha_vector;					//This will be used to bias larger total protein amounts for larger edges.
	alpha_vector.resize(number_of_edges);
	int edge_index = 0;
	vector<int> which_cell;							//Used to store an indicator of which element of planarpolarity to access.
	which_cell.resize(number_of_edges);

	//Used to convert pointers from edges to pointers to PCPEdges
	vector<PCPEdge*> edge_list_this_cell;
	edge_list_this_cell.resize(number_of_edges);
	for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
		edge_list_this_cell.at(edge_index) = dynamic_cast<PCPEdge*>(edges.at(edge_index));
	}

	//Set the type of the newly cell to be 1 (which means that it is differentiated).
	//Output an error message if the cell is already differentiated.
	if ((*this).is_Cell_with_PCP()){
		cout << "ERROR: You are reinitializing proteins in an already differentiated cell. Tell HN." << endl;
	}
	type = pcp_cell_type;

	for (edge_index = 0; edge_index != number_of_edges; edge_index++){
		if (((*(*edges[edge_index]).cell1)).id == (*this).id){
			which_cell[edge_index] = 0;
		}
		else {
			which_cell[edge_index] = 1;
		}
	}

	(*this).SavedInfoForEvolve_CWVerts = (*this).GetClockWiseVertices();
	PCPEdge* horizontal_start_edge = (*this).FindUpperHorizontalEdge();
	PCPEdge* next_edge_ptr = (*this).GetNextCWEdge(horizontal_start_edge);
	PCPEdge* next_next_edge_ptr = (*this).GetNextCWEdge(next_edge_ptr);
	PCPEdge* next_next_next_edge_ptr = (*this).GetNextCWEdge(next_next_edge_ptr);
	PCPEdge* next_next_next_next_edge_ptr = (*this).GetNextCWEdge(next_next_next_edge_ptr);
	PCPEdge* next_next_next_next_next_edge_ptr = (*this).GetNextCWEdge(next_next_next_next_edge_ptr);

	if (simulation_type == 0){
		double high_protein_amount = 0.6;
		double low_protein_amount = 0.1;
		double medium_protein_amount = 0.3 / 4;

		for (edge_index = 0; edge_index != number_of_edges; edge_index++){
			if ((*edge_list_this_cell[edge_index]).id == (*horizontal_start_edge).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = high_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_edge_ptr).id || (*edge_list_this_cell[edge_index]).id == (*next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = medium_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = medium_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_edge_ptr).id || (*edge_list_this_cell[edge_index]).id == (*next_next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = medium_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = medium_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amount;
			}
		}
	}
	else if (simulation_type == 1){
		double high_protein_amount = 0.3;
		double low_protein_amount = 0.05;
		double medium_protein_amount = 0.15;

		for (edge_index = 0; edge_index != number_of_edges; edge_index++){
			if ((*edge_list_this_cell[edge_index]).id == (*horizontal_start_edge).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = medium_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = medium_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_edge_ptr).id || (*edge_list_this_cell[edge_index]).id == (*next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = high_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_edge_ptr).id || (*edge_list_this_cell[edge_index]).id == (*next_next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = medium_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = medium_protein_amount;
			}
		}
	}
	else if (simulation_type == 2){
		double high_protein_amount = 0.6;
		double low_protein_amount = 0.05;
		double medium_protein_amount = 0.2;

		for (edge_index = 0; edge_index != number_of_edges; edge_index++){
			if ((*edge_list_this_cell[edge_index]).id == (*horizontal_start_edge).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = high_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = medium_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = medium_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
		}
	}
	else if (simulation_type == 3){
		double high_protein_amountA = 0.2;
		double low_protein_amountA = 0.15;

		double high_protein_amountB = 0.2;
		double low_protein_amountB = 0.1;

		for (edge_index = 0; edge_index != number_of_edges; edge_index++){
			if ((*edge_list_this_cell[edge_index]).id == (*horizontal_start_edge).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = high_protein_amountA;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amountB;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amountA;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amountB;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amountA;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amountB;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = high_protein_amountA;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amountB;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amountA;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amountB;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amountA;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amountB;
			}
		}
	}
	else if (simulation_type == 4){
		double high_protein_amount = 0.3;
		double medium_protein_amount = 0.15;
		double low_protein_amount = 0.05;

		for (edge_index = 0; edge_index != number_of_edges; edge_index++){
			if ((*edge_list_this_cell[edge_index]).id == (*horizontal_start_edge).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = high_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = medium_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = medium_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = high_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = low_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = low_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = high_protein_amount;
			}
			else if ((*edge_list_this_cell[edge_index]).id == (*next_next_next_next_next_edge_ptr).id){
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = medium_protein_amount;
				(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = medium_protein_amount;
			}
		}
	}

	for (edge_index = 0; edge_index != number_of_edges; edge_index++){
		(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] * (1 - delta_protein_amount);
		(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] * (1 + delta_protein_amount);
	}
}

//! This function is used to randomized the protein distribution within a cell about the specified protein distribution (from input file). Called by ReadInputData.
void PCPCell::HelperReadInputData_RandomizePolarityProteins(){
	int number_of_edges = edges.size();				//How many edges does the current cells have?
	vector<double> alpha_params_vectorA;			//Store the length of each edge (in order to later define the concentration)
	alpha_params_vectorA.resize(number_of_edges);
	double sum_proteinA = 0.0;
	double sum_proteinA_after = 0.0;
	vector<double> alpha_params_vectorB;			//Store the length of each edge (in order to later define the concentration)
	alpha_params_vectorB.resize(number_of_edges);
	double sum_proteinB = 0.0;
	double sum_proteinB_after = 0.0;
	vector<double> proteinA_amounts;				//Vector for randomly generated amounts of protein A for each edge.
	proteinA_amounts.resize(number_of_edges);
	vector<double> proteinB_amounts;				//Vector for randomly generated amounts of protein B for each edge.
	proteinB_amounts.resize(number_of_edges);
	int edge_index = 0;
	vector<int> which_cell;							//Used to store an indicator of which element of planarpolarity to access.
	which_cell.resize(number_of_edges);

	//Used to convert pointers from edges to pointers to PCPEdges
	vector<PCPEdge*> edge_list_this_cell;
	edge_list_this_cell.resize(number_of_edges);
	for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
		edge_list_this_cell.at(edge_index) = dynamic_cast<PCPEdge*>(edges.at(edge_index));
	}

	type = pcp_cell_type;

	for (edge_index = 0; edge_index != number_of_edges; edge_index++){
		if ((*(*edges[edge_index]).cell1) == *this){
			which_cell[edge_index] = 0;
		}
		else {
			which_cell[edge_index] = 1;
		}

		sum_proteinA += (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]];
		sum_proteinB += (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1];
	}

	for (edge_index = 0; edge_index != number_of_edges; edge_index++){
		if ((*(*edges[edge_index]).cell1) == *this){
			which_cell[edge_index] = 0;
		}
		else {
			which_cell[edge_index] = 1;
		}

		alpha_params_vectorA[edge_index] = dirichlet_noise_param * (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] / sum_proteinA;
		alpha_params_vectorB[edge_index] = dirichlet_noise_param * (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] / sum_proteinB;
	}
	double max_protein_concentration = 0.0;

	do{
		max_protein_concentration = 0.0;
		//do random assignment of protein A amounts
		gsl_ran_dirichlet(polarity_random_generator, number_of_edges, &alpha_params_vectorA[0], &proteinA_amounts[0]);
		//do random assignments of protein B amounts
		gsl_ran_dirichlet(polarity_random_generator, number_of_edges, &alpha_params_vectorB[0], &proteinB_amounts[0]);
		for (edge_index = 0; edge_index != edges.size(); edge_index++){
			(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = proteinA_amounts[edge_index] * sum_proteinA;
			(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = proteinB_amounts[edge_index] * sum_proteinB;
			if (proteinA_amounts[edge_index] * sum_proteinA >= max_protein_concentration){
				max_protein_concentration = proteinA_amounts[edge_index] * sum_proteinA;
			}
			if (proteinB_amounts[edge_index] * sum_proteinB >= max_protein_concentration){
				max_protein_concentration = proteinB_amounts[edge_index] * sum_proteinB;
			}
		}
	} while (max_protein_concentration > max_concentration);
}

//! This function sets all protein concentrations in a cell to be identical to those of cell pointed to by cell_ptr (only applicable in case of hexagonal field of cells).
void PCPCell::HelperReadInputData_SetProteinConc_fromCellPointer(PCPCell* cell_ptr){
	int number_of_edges = edges.size();				//How many edges does the current cells have?
	vector<double> proteinA_amounts;				//Vector for randomly generated amounts of protein A for each edge.
	proteinA_amounts.resize(number_of_edges);
	vector<double> proteinB_amounts;				//Vector for randomly generated amounts of protein B for each edge.
	proteinB_amounts.resize(number_of_edges);
	int edge_index = 0;
	vector<int> which_cell;							//Used to store an indicator of which element of planarpolarity to access.
	which_cell.resize(number_of_edges);

	(*cell_ptr).SavedInfoForEvolve_CWVerts = (*cell_ptr).GetClockWiseVertices();
	vector<PCPEdge*> edge_list_this_cell = (*cell_ptr).OrderEdges();

	type = pcp_cell_type;

	for (edge_index = 0; edge_index != number_of_edges; edge_index++){
		if ((*(*edge_list_this_cell[edge_index]).cell1).id == (*cell_ptr).id){
			which_cell[edge_index] = 0;
		}
		else {
			which_cell[edge_index] = 1;
		}
		proteinA_amounts[edge_index] = (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]];
		proteinB_amounts[edge_index] = (*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1];
	}

	(*this).SavedInfoForEvolve_CWVerts = (*this).GetClockWiseVertices();
	edge_list_this_cell = (*this).OrderEdges();

	for (edge_index = 0; edge_index != number_of_edges; edge_index++){
		if ((*(*edge_list_this_cell[edge_index]).cell1) == *this){
			which_cell[edge_index] = 0;
		}
		else {
			which_cell[edge_index] = 1;
		}
		(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = proteinA_amounts[edge_index];
		(*edge_list_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = proteinB_amounts[edge_index];
	}
}