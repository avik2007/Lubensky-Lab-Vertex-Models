/**	\file	PCPTissueClass.cpp
*	Implementation of PCPTissue class.
*	Contains function definitions associated with PCPTissue objects.
*	Depends on base class Tissue.
*/

//General files
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
#include <iostream>
#include <cassert>					//for error checking/assertion
#include <iomanip>					//for changing output precision

//Base class files
#include "Constants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"

//PCP class files
#include "PCPCellClass.h"
#include "PCPEdgeClass.h"
#include "PCPTissueClass.h"
#include "PCPConstants.h"
#include "PCPEvolutionFunction.h"

using namespace std;

//! Creates a new edge and adds it to the list of edges, returns a pointer to that edge: initialize PCP protein concentrations to zero.
Edge* PCPTissue::CreateEdge(Edge* subclass, bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2){
	int largestID, newid;
	Edge newedge;
	//Get the id
	largestID = (*max_element(edgeIDs.begin(), edgeIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;													//Now make a new ID that is unique

	if (newid > 1000000000) { cout << "ERROR to many edges" << endl; }
	//create the vertex and add the vertex to the lists
	PCPEdge* eptr = new PCPEdge;
	(*eptr).id = newid;
	(*eptr).xflag = xflg;
	(*eptr).yflag = yflg;
	(*eptr).vertex1 = v1;
	(*eptr).vertex2 = v2;
	(*eptr).cell1 = c1;
	(*eptr).cell2 = c2;
	(edges).push_back(eptr);
	edgeIDs.push_back(newid);
	//Initialize PCP concentrations to be equal to zero.
	(*eptr).planarpolarity.resize(4);
	fill((*eptr).planarpolarity.begin(), (*eptr).planarpolarity.end(), 0.0);

	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){
		PCPEdge* pcpeit = dynamic_cast<PCPEdge*>(*eit);
		if ((*(pcpeit)).id == newid){
			return &(*(pcpeit));
			break;
		}
	}
}

//Cut line of cells (by specifying cells that intersect with a vertical line).
//Any cell in contact with this line of cells on one side is initialized with PCP.
void PCPTissue::AddDifferentiatedCellColumn(){
	//This cut will be used to cut the tissue. Any cell intersected by this cut will be of type intersect_type.
	//Then, all cells bordering a cell of type intersect_type to the right will get PCP (become type pcp_cell_type).

	//Define the value of the cutline.
	cutline_pcp_initialization = LX / 4.0;

	//CUT THE TISSUE AT cutline_pcp_initialization, find all cells intersected by the cut.
	//Note: This replicates code that is found in CalculateTissueStressX and CalculateTissueStressY
	vector<Edge*> crossingEdges;											//vector of all crossing edges
	crossingEdges.reserve(edges.size());
	unordered_map<long, Cell*>	crossingCells_id2pointer;							//map of crossing cell id to its pointer
	crossingCells_id2pointer.reserve(edges.size() / 2);
	unordered_map<long, vector<Edge*>> crossingCells_id2crossingEdges;						//map of crossing cell id to its crossing edges
	crossingCells_id2crossingEdges.reserve(edges.size() / 2);

	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){
		Vertex* v1 = (*(*eit)).vertex1;
		Vertex* v2 = (*(*eit)).vertex2;
		double x1 = (*v1).r.x;
		double x2 = (*v2).r.x;
		bool xflag = (*(*eit)).xflag;										//for PBC

		//test if edge crosses cutline_pcp_initialization (taking PBC into account)
		bool isCrossing = (!xflag && x1 < cutline_pcp_initialization && x2 > cutline_pcp_initialization) || (!xflag && x1 > cutline_pcp_initialization && x2 < cutline_pcp_initialization)		//if xflag = false and x1 & x2 are on the OPPOSITE side of xCut
			|| (xflag && x1 < cutline_pcp_initialization && x2 < cutline_pcp_initialization) || (xflag && x1 > cutline_pcp_initialization && x2 > cutline_pcp_initialization);			//if xflag = true and x1 & x2 are on the SAME side of xCut

		if (isCrossing){
			crossingEdges.push_back(&(*(*eit)));
			long cid1 = (*(*(*eit)).cell1).id;
			long cid2 = (*(*(*eit)).cell2).id;
			if (crossingCells_id2pointer.count(cid1) == 0){							//if not already there, put the pointer to the crossing cell into the map
				crossingCells_id2pointer[cid1] = (*(*eit)).cell1;
			}
			if (crossingCells_id2pointer.count(cid2) == 0){
				crossingCells_id2pointer[cid2] = (*(*eit)).cell2;
			}
		}
	}

	//Loop through crossed cells and specify them as type intersect_type (serves as a template for specifying the first cells with PCP).
	for (unordered_map<long, Cell*>::iterator mapit = crossingCells_id2pointer.begin(); mapit != crossingCells_id2pointer.end(); mapit++){
		Cell* cellptr = (*mapit).second;
		(*cellptr).type = intersect_type;
	}

	//Now, define the next round of cells for differentiation. This is done in order to ensure the propogation of PCP in only one direction.
	vector<Vertex*> boundary_vertices; 	//contains pointers to vertices along one of the boundaries (in order)
	Vertex *v0; 						//Used to make sure that when searching for vertices, the search never goes back to a previously found vertex
	int i, j, l; 						//indices in the for loops
	vector<Edge*> edgecandidates; 		//The candidate for the boundary edges that have a yflag

	list<Edge*>::iterator edgeIter;
	//SEARCH FOR EDGES THAT CROSS THE Y-BOUNDARY AND HAVE ONE CELL OF TYPE standard_cell_type and ONE CELL OF TYPE intersect_type.
	for (edgeIter = edges.begin(); edgeIter != edges.end(); edgeIter++){
		Cell* cell1ptr = (*(*edgeIter)).cell1;
		Cell* cell2ptr = (*(*edgeIter)).cell2;
		int edge_y_flag = (*(*edgeIter)).yflag;
		if (edge_y_flag == 1 && (((*cell1ptr).type == intersect_type && (*cell2ptr).type == standard_cell_type) || ((*cell1ptr).type == standard_cell_type && (*cell2ptr).type == intersect_type))){
			edgecandidates.push_back(&(*(*edgeIter))); 	//push back a pointer to the edge if it has a yflag and is at the boundary of propagation.
		}
	}

	//Look through edge candidates. Find an edge that faces to the right.
	Edge* edge_for_propagating_right = edge_candidate_to_right(edgecandidates);

	//Push back one of the vertices of the first edge: this edge faces to the right.
	boundary_vertices.push_back((*edge_for_propagating_right).vertex1);

	//v0 is initialized here: this will be used to make sure that as you go through vertices along the boundary, you never go backward.
	v0 = (*edge_for_propagating_right).vertex2;
	//l is a counter for the elements of boundary_vertices
	l = 0;
	do{
		//Loop through the the edges for boundary_vertices[l]
		for (i = 0; i != (*boundary_vertices[l]).edges.size(); i++){
			//find boundary edge (by checking for the types of cell1 and cell 2 for this edge).
			//then find vertex on other side of edge, never go back to previously found vertices.
			Cell* cell1ptr = (*(*boundary_vertices[l]).edges[i]).cell1;
			Cell* cell2ptr = (*(*boundary_vertices[l]).edges[i]).cell2;
			Vertex* other_vertex = (*(*boundary_vertices[l]).edges[i]).GetOtherVertex(boundary_vertices[l]);
			if ((((*cell1ptr).type == intersect_type && (*cell2ptr).type == standard_cell_type) || ((*cell1ptr).type == standard_cell_type && (*cell2ptr).type == intersect_type)) && other_vertex != v0){
				//push back vertex that you found, define v0 to prevent going backwards on the boundary, counter++
				boundary_vertices.push_back(other_vertex);
				v0 = boundary_vertices[l];
				l++;
				break;
			}
		}
	} while (boundary_vertices[0] != boundary_vertices[l]);		//Stop when you duplicate the first element

	//first_diff_cells is a list of cells to have PCP first.
	vector<PCPCell*> first_diff_cells;
	int number_boundary_vertices = boundary_vertices.size();
	for (i = 0; i != number_boundary_vertices; i++){
		//Find list of cells connected to all of the vertices on the boundary.
		vector<Cell*> list_of_cells = (*boundary_vertices[i]).GetCells();
		int length_of_list_of_cells = list_of_cells.size();
		//Loop through the cells attached to boundary_vertices[i]
		for (j = 0; j != length_of_list_of_cells; j++){
			//If any cell is standard_cell_type, push it back.
			if ((*list_of_cells[j]).type == standard_cell_type){
				first_diff_cells.push_back(dynamic_cast<PCPCell*>(list_of_cells[j]));
			}
		}
	}

	//sort the pointers.
	std::sort(first_diff_cells.begin(), first_diff_cells.end());

	//Only keep unique cells: delete duplicates.
	auto last = std::unique(first_diff_cells.begin(), first_diff_cells.end());
	first_diff_cells.erase(last, first_diff_cells.end());

	//! Order cell pointer by vertical position (try to identify chains of cells).
	vector<PCPCell*> ordered_first_diff_cells = order_cells_position(first_diff_cells);
	vector<PCPCell*> corrected_ordered_first_diff_cells;

	//Check if a cell is trapped (should not be differentiated).
	vector<bool> check_for_faulty_cells;
	int number_of_first_diff_cells = first_diff_cells.size();
	check_for_faulty_cells.resize(number_of_first_diff_cells);

	//Specify the first_diff_cells as type pcp_cell_type (will have PCP), just to check if any cell will be entirely surrounded by cells of type pcp_cell_type or type intersect_type.
	//If a cell is surrounded by cells of type pcp_cell_type or type intersect_type, but is supposed to get PCP, forbid it from getting PCP.
	for (i = 0; i != number_of_first_diff_cells; i++){
		(*ordered_first_diff_cells[i]).type = pcp_cell_type;
	}
	for (i = 0; i != number_of_first_diff_cells; i++){
		check_for_faulty_cells.at(i) = (*this).check_for_faulty_cell(ordered_first_diff_cells[i]);
	}

	//Loop through and only push back cells that are not trapped. Label trapped cells as intersect type.
	for (i = 0; i != number_of_first_diff_cells; i++){
		if (!check_for_faulty_cells.at(i)){
			corrected_ordered_first_diff_cells.push_back(ordered_first_diff_cells[i]);
		}
		else{
			(*ordered_first_diff_cells[i]).type = intersect_type;
		}
	}

	//Add PCP to cells.
	int number_of_corrected_first_diff_cells = corrected_ordered_first_diff_cells.size();
	for (i = 0; i != number_of_corrected_first_diff_cells; i++){
		//Change the cell type for the next cells to differentiate to type pcp_cell_type.
		(*corrected_ordered_first_diff_cells[i]).HelperEvolve_InitializePolarityProteins();
	}

	//Establish checkerboard types.
	for (i = 0; i != number_of_corrected_first_diff_cells; i++){
		if (i % 2 == 0){
			(*corrected_ordered_first_diff_cells[i]).type = checkerboard_type1;
		}
		else{
			(*corrected_ordered_first_diff_cells[i]).type = checkerboard_type2;
		}
	}
}

//! This function calculates the total number of edges in differentiated cells in the tissue (half of the number of dynamical protein variables)
int PCPTissue::Find_Number_Differentiated_Edges(){
	int cell_edge_counter = 0;
	list<Cell*>::iterator cit;
	PCPCell* pcpcit;
	for (cit = cells.begin(); cit != cells.end(); cit++){	//look through all of the cells
		pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			cell_edge_counter += (*(*cit)).edges.size();
		}
	}
	return cell_edge_counter;
}

//! This function calculates the total number of cells with PCP in the tissue
int PCPTissue::Find_Number_Differentiated_Cells(){
	int pcp_cell_counter = 0;
	list<Cell*>::iterator cit;
	PCPCell* pcpcit;
	for (cit = cells.begin(); cit != cells.end(); cit++){	//look through all of the cells
		pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			pcp_cell_counter++;
		}
	}
	return pcp_cell_counter;
}

// Reads and constructs simulation from file. The file must be formatted like InputTemplate.txt. This function parses the input file then initializes all the objects and pointers.
bool PCPTissue::ReadInputData(string InputFileName, bool includeConcentrations){
	cout << "Entered ReadInputData" << endl;
	ifstream datafile;
	datafile.open(InputFileName);											//open file
	if (!datafile.is_open()){											//check if open file successful
		cout << "ERROR cannot open file: " << InputFileName << endl;
		return 0;
	}

	cout << InputFileName << " file opened sucessfully" << endl;
	(*this).HelperReadInputData_SetGlobalValues(&datafile);
	cout << "Global values Set" << endl;
	(*this).HelperReadInputData_CreateObjects(&datafile);
	(*this).HelperReadInputData_CreatePointers();
	(*this).HelperReadInputData_CheckTypes();

	if (includeConcentrations){
		//Start back at the beginning of the file and read in initial protein concentrations from the file.
		datafile.close();
		datafile.open(InputFileName);
		(*this).HelperReadInputData_ProteinConcentrations(&datafile);
		int i = 0;
		PCPCell* first_pcpcit;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
			PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
			(*pcpcit).HelperReadInputData_RandomizePolarityProteins();

			//NOTE: use this to add noise to concentrations in such a way that each cell remains identical.
			//if (i == 0){
			//	(*pcpcit).HelperReadInputData_RandomizePolarityProteins();
			//	first_pcpcit = dynamic_cast<PCPCell*>(*cit);
			//	i++;
			//}
			//else{
			//	(*pcpcit).HelperReadInputData_SetProteinConc_fromCellPointer(first_pcpcit);
			//}
		}
	}
	else {
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
			PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
			(*pcpcit).HelperReadInputData_InitializePolarityProteins();
		}
	}

	datafile.close();
	return 1;
}

// This function calculates the PCP free energy of a configuration (stores the free energy calculation in the location pointed to by protein_energy).
// Currently used within evolve and stores value into protein_energy.
double PCPTissue::CalculatePCPFreeEnergy(void* params){
	//to find adjacent edges, use function to find list of edges based on list of vertices
	double free_energy;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;

	//This vector is important. It allows you to figure out what row or column index refers to a given cell's edges in the loop for defining U below.
	//Given some cell index, this vector will give you where the entries for this cell begin.
	vector<int> start_index_cells;
	vector<PCPEdge*> edge_list_for_this_cell;
	//Resize this vector such that it fits all cells in the tissue.
	start_index_cells.resize((*this).cells.size());
	//Initialize start_index to zero.
	int start_index = 0;
	//Use ind to refer to the current cell index.
	int cell_index = 0;
	for (list<Cell*>::iterator cit = (*this).cells.begin(); cit != (*this).cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//Resize the edge list to fit this given cell.
			edge_list_for_this_cell.resize((*(*cit)).edges.size());
			//Store a vector of pointers to the edges in the current cell.
			edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			//Store the start index for the current cell.
			start_index_cells[cell_index] = start_index;
			start_index += edge_list_for_this_cell.size();
		}
		//Increment the start_index and cell index (ind).
		cell_index++;
	}

	//Initialize protein_free_energy to 0.0
	free_energy = 0.0;
	int number_of_edges;
	//Helps to identify the index associated with the concentration on the other side of an edge.
	Cell* column_cell_pointer;
	//Helps to identify the index associated with the concentration on the other side of an edge.
	int beta_index;
	//Store protein concentrations->note: this is different from y_protein[]
	vector<double> protein_concentrations;
	protein_concentrations.resize((*(combinedparameters*)params).dimension_protein);
	int type_dimension = (*(combinedparameters*)params).dimension_protein / 2;
	int cell_id;
	//Loop through all of the cells and store the corresponding protein concentrations.
	for (list<Cell*>::iterator cit = (*this).cells.begin(); cit != (*this).cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			start_index = start_index_cells[cell_index];
			edge_list_for_this_cell.resize((*(*cit)).edges.size());
			//order edges: use casted pcp cell pointer to extract ordered edges.
			edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			number_of_edges = edge_list_for_this_cell.size();
			vector<double> edge_length_list;
			edge_length_list.resize(number_of_edges);
			//Store a list of edges lengths for this cell.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
			}
			//build the vector of protein concentrations (following the same conventions as the U matrix)
			for (int i = 0; i != number_of_edges; i++){
				if ((*edge_list_for_this_cell[i]).GetThisCellNumber(*cit) == 1){
					//Store the concentrations in a vector (i.e., protein amount/edge length).
					protein_concentrations[start_index + i] = (*edge_list_for_this_cell[i]).planarpolarity[0] / edge_length_list[i];
					protein_concentrations[start_index + i + type_dimension] = (*edge_list_for_this_cell[i]).planarpolarity[1] / edge_length_list[i];
				}
				else{
					//Store the concentrations in a vector (i.e., protein amount/edge length).
					protein_concentrations[start_index + i] = (*edge_list_for_this_cell[i]).planarpolarity[2] / edge_length_list[i];
					protein_concentrations[start_index + i + type_dimension] = (*edge_list_for_this_cell[i]).planarpolarity[3] / edge_length_list[i];
				}
			}
		}
	}

	int column_cell_id;
	int column_cell_index;
	vector<PCPEdge*> column_cell_ordered_edges;
	for (list<Cell*>::iterator cit = (*this).cells.begin(); cit != (*this).cells.end(); cit++){
		//Only add to protein free energy for differentiated cells (with non-zero protein)
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			start_index = start_index_cells[cell_index];
			edge_list_for_this_cell.resize((*(*cit)).edges.size());
			//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
			edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			number_of_edges = edge_list_for_this_cell.size();
			vector<double> edge_length_list;
			edge_length_list.resize(number_of_edges);
			//Store a list of edges lengths for this cell.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
			}
			//This is used to store the index for the cell on the other side of all edges.
			vector<int> other_cell_index;
			other_cell_index.resize(number_of_edges);
			//build the vector of protein concentrations (following the same conventions as the U matrix)
			for (int i = 0; i != number_of_edges; i++){
				if ((*edge_list_for_this_cell[i]).GetThisCellNumber(*cit) == 1){
					//Store the index for the cell on the other side of the edge.
					other_cell_index[i] = (*(combinedparameters*)params).cell_id2index.at((*(*edge_list_for_this_cell[i]).cell2).id);
				}
				else{
					//Store the index for the cell on the other side of the edge.
					other_cell_index[i] = (*(combinedparameters*)params).cell_id2index.at((*(*edge_list_for_this_cell[i]).cell1).id);
				}
			}
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				//For the edge referred to by the current column index, look to the other side of the edge.
				if ((*edge_list_for_this_cell[edge_index]).GetThisCellNumber(*cit) == 1){
					column_cell_pointer = (*edge_list_for_this_cell[edge_index]).cell2;
				}
				else {
					column_cell_pointer = (*edge_list_for_this_cell[edge_index]).cell1;
				}
				//Cast the cell on the other side of the current edge to a PCPCell pointer.
				PCPCell* pcp_column_cell_ptr = dynamic_cast<PCPCell*>(column_cell_pointer);
				column_cell_id = (*column_cell_pointer).id;
				column_cell_index = cell_id2index.at(column_cell_id);
				//For the cell on the other side of the edge, order the edges of that cell and find which edge corresponds to this edge.
				column_cell_ordered_edges.resize((*column_cell_pointer).edges.size());
				column_cell_ordered_edges = (*pcp_column_cell_ptr).FindOrderedEdges();
				//if cell on the other side of the edge is differentiated, enter this statement.
				if ((*pcp_column_cell_ptr).is_Cell_with_PCP()){
					//Find the index of that edge in the other cell's list.
					for (int i = 0; i != (*column_cell_pointer).edges.size(); i++){
						if (*column_cell_ordered_edges[i] == *edge_list_for_this_cell[edge_index]){
							beta_index = i;
						}
					}
					//probably need factor of 1/2 in front of epstwo (for consistency).->not included in PLOS 2012 though.
					free_energy += -epsone * edge_length_list[edge_index] * protein_concentrations[start_index + edge_index] * protein_concentrations[type_dimension + start_index_cells[other_cell_index[edge_index]] + beta_index] - epstwo / 2 * edge_length_list[edge_index] * (protein_concentrations[start_index + edge_index] * protein_concentrations[start_index_cells[other_cell_index[edge_index]] + beta_index] + protein_concentrations[type_dimension + start_index + edge_index] * protein_concentrations[type_dimension + start_index_cells[other_cell_index[edge_index]] + beta_index]);
				}
				free_energy += epsab * edge_length_list[edge_index] * protein_concentrations[start_index + edge_index] * protein_concentrations[start_index + edge_index + type_dimension] - epsnaught / 2 * edge_length_list[edge_index] * (pow(protein_concentrations[start_index + edge_index], 2) + pow(protein_concentrations[start_index + edge_index + type_dimension], 2));
				free_energy += te * edge_length_list[edge_index] * (protein_concentrations[start_index + edge_index] * log(protein_concentrations[start_index + edge_index]) + protein_concentrations[type_dimension + start_index + edge_index] * log(protein_concentrations[type_dimension + start_index + edge_index]));
			}
			for (int first_edge_index = 0; first_edge_index != number_of_edges; first_edge_index++){
				int second_edge_index = (first_edge_index + 1) % number_of_edges;
				free_energy += jone / 2 * (pow(protein_concentrations[start_index + first_edge_index] - protein_concentrations[start_index + second_edge_index], 2) + pow(protein_concentrations[type_dimension + start_index + first_edge_index] - protein_concentrations[type_dimension + start_index + second_edge_index], 2));
				free_energy -= jtwo / 2 * (pow(protein_concentrations[start_index + first_edge_index] - protein_concentrations[type_dimension + start_index + second_edge_index], 2) + pow(protein_concentrations[start_index + first_edge_index + type_dimension] - protein_concentrations[start_index + second_edge_index], 2));
			}
		}
	}

	return free_energy;
}

//Integrate the mechanical dynamics using the explicit implementation of the imex scheme (second-order).
void PCPTissue::Integrate_IMEX_Second_Order(double this_step_size, void* params, ofstream* errorfile){
	//Define coefficients for solving method.
	double gamma = 1. / 2.;
	double c = 1. / 8.;

	//Define ratios of step sizes.
	double omega_n_plus_1 = this_step_size / past_step_sizes[1];

	//Define coefficients for third order equations.
	double alpha0n = (2. * gamma - 1.) * pow(omega_n_plus_1, 2.) / (1. + omega_n_plus_1);
	double alpha1n = (1. - 2. * gamma) * omega_n_plus_1 - 1.;
	double alpha2n = (1. + 2. * gamma * omega_n_plus_1) / (1. + omega_n_plus_1);
	double beta0n = -gamma * omega_n_plus_1;
	double beta1n = 1. + gamma * omega_n_plus_1;
	double gamma0n = c / 2.;
	double gamma1n = 1. - gamma - (1. + 1. / omega_n_plus_1) * c / 2.;
	double gamma2n = gamma + c / 2. / omega_n_plus_1;

	//Calculate the force on every vertex 
	(*this).ForceFunction(errorfile);

	int i = 0;
	//loop through all the vertices and set rnew
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
		(*(*vit)).rnew = (past_mech_state[1][i] * (-alpha0n) + (*(*vit)).r * (-alpha1n) + f_mech_step_before[1][i] * beta0n * this_step_size / alphaM + (*(*vit)).SavedInfoForEvolve_TotalForce * beta1n * this_step_size / alphaM) * 1. / alpha2n;
		i++;
	}
	LXnew = (past_mech_state[1][i].x * (-alpha0n) + LX * (-alpha1n) + f_mech_step_before[1][i].x * beta0n * this_step_size / TISSUE_ALPHA + ForceLX * beta1n * this_step_size / TISSUE_ALPHA) * 1. / alpha2n;
	LYnew = (past_mech_state[1][i].y * (-alpha0n) + LY * (-alpha1n) + f_mech_step_before[1][i].y * beta0n * this_step_size / TISSUE_ALPHA + ForceLY * beta1n * this_step_size / TISSUE_ALPHA) * 1. / alpha2n;
}

//Integrate the mechanical dynamics using the explicit implementation of the imex scheme (third-order).
void PCPTissue::Integrate_IMEX_Third_Order(double this_step_size, void* params, ofstream* errorfile){
	double gamma = 0.5;
	double theta = -1.25;
	double c = -0.52;

	//Define ratios of step sizes.
	double omega_n_plus_2 = this_step_size / past_step_sizes[1];
	double omega_n_plus_1 = past_step_sizes[1] / past_step_sizes[0];

	//Define coefficients for third order equations.
	double alpha0n = -1. * pow(omega_n_plus_2, 2.) * pow(omega_n_plus_1, 3.) * (3. * pow(gamma, 2.) * omega_n_plus_2 + 2. * gamma * (1. - omega_n_plus_2) - 1.) / (1. + omega_n_plus_1) / (1. + omega_n_plus_1 * (1. + omega_n_plus_2));
	double alpha1n = -1. * pow(omega_n_plus_2, 2.) * (gamma * omega_n_plus_1 * omega_n_plus_2 * (2. - 3. * gamma) + (1. - 2. * gamma) * (1. + omega_n_plus_1)) / (1. + omega_n_plus_2);
	double alpha2n = -1. * (omega_n_plus_1 * (1. + gamma*omega_n_plus_2) * (1. + omega_n_plus_2 * (3. * gamma - 2.)) + omega_n_plus_2 * (2. * gamma - 1.) + 1.) / (1. + omega_n_plus_1) - theta;
	double alpha3n = (1. + 2. * gamma*omega_n_plus_2 + omega_n_plus_1 * (1. + omega_n_plus_2*gamma)*(1. + 3. * gamma*omega_n_plus_2)) / (1. + omega_n_plus_2) / (1. + omega_n_plus_1*(1. + omega_n_plus_2)) + theta;
	double beta0n = pow(omega_n_plus_1, 2.) * omega_n_plus_2 * (6. * gamma * (1. + gamma*omega_n_plus_2) + theta * (3. + 2. * omega_n_plus_2)) / 6. / (1. + omega_n_plus_1);
	double beta1n = -gamma * omega_n_plus_2 * (1. + omega_n_plus_1 * (1. + gamma*omega_n_plus_2)) - (theta * omega_n_plus_2 * (3. + omega_n_plus_1 * (3. + 2. * omega_n_plus_2))) / 6.;
	double beta2n = (1. + gamma * omega_n_plus_2) * (1. + omega_n_plus_1 * (1. + gamma * omega_n_plus_2)) / (1. + omega_n_plus_1) + theta * (1. + omega_n_plus_2 / 2. + omega_n_plus_2 * omega_n_plus_1 * (3. + 2. * omega_n_plus_2) / 6. / (1. + omega_n_plus_1));
	double gamma0n = theta * pow(omega_n_plus_1, 2.) * omega_n_plus_2 * (3. + 2. * omega_n_plus_2) / 6. / (1. + omega_n_plus_1) - c;
	double gamma1n = (c * (1. + omega_n_plus_1) * (1. + omega_n_plus_1 * (1. + omega_n_plus_2)) - pow(omega_n_plus_1, 2.) * pow(omega_n_plus_2, 2.) * gamma * (1. - gamma)) / pow(omega_n_plus_1, 2.) / (1. + omega_n_plus_2) - theta * omega_n_plus_2 * (3. + omega_n_plus_1 * (3. + 2. * omega_n_plus_2)) / 6.;
	double gamma2n = (pow(omega_n_plus_1, 2.) * omega_n_plus_2 * (1. - gamma) * (1. + gamma*omega_n_plus_2) - c * (1. + omega_n_plus_1 * (1. + omega_n_plus_2))) / pow(omega_n_plus_1, 2.) / omega_n_plus_2 + theta * (1. + omega_n_plus_2 / 2. + omega_n_plus_2 * omega_n_plus_1 * (3. + 2. * omega_n_plus_2) / 6. / (1. + omega_n_plus_1));
	double gamma3n = (pow(omega_n_plus_1, 2.) * omega_n_plus_2 * gamma * (1. + gamma*omega_n_plus_2) + c * (1. + omega_n_plus_1)) / pow(omega_n_plus_1, 2.) / omega_n_plus_2 / (1. + omega_n_plus_2);

	//Calculate the force on every vertex 
	(*this).ForceFunction(errorfile);

	int i = 0;
	//loop through all the vertices and set rnew
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
		(*(*vit)).rnew = (past_mech_state[0][i] * (-alpha0n) + past_mech_state[1][i] * (-alpha1n) + (*(*vit)).r * (-alpha2n) + f_mech_step_before[0][i] * beta0n * this_step_size / alphaM + f_mech_step_before[1][i] * beta1n * this_step_size / alphaM + (*(*vit)).SavedInfoForEvolve_TotalForce * beta2n * this_step_size / alphaM) * 1. / alpha3n;
		i++;
	}
	LXnew = (past_mech_state[0][i].x * (-alpha0n) + past_mech_state[1][i].x * (-alpha1n) + LX * (-alpha2n) + f_mech_step_before[0][i].x * beta0n * this_step_size / TISSUE_ALPHA + f_mech_step_before[1][i].x * beta1n * this_step_size / TISSUE_ALPHA + ForceLX * beta2n * this_step_size / TISSUE_ALPHA) * 1. / alpha3n;
	LYnew = (past_mech_state[0][i].y * (-alpha0n) + past_mech_state[1][i].y * (-alpha1n) + LY * (-alpha2n) + f_mech_step_before[0][i].y * beta0n * this_step_size / TISSUE_ALPHA + f_mech_step_before[1][i].y * beta1n * this_step_size / TISSUE_ALPHA + ForceLY * beta2n * this_step_size / TISSUE_ALPHA) * 1. / alpha3n;
}

//Umbrella function for integrating the mechanical dynamics (immediately handles first-second-third order subtleties).
void PCPTissue::Integrate_IMEX_Without_Protein(double* step_size, const double mech_abs_tolerance, double* current_time, void* params, ofstream* errorfile){
	int first_order_check = 0;
	int second_order_check = 0;
	int dim_mechanical = (*(combinedparameters*)params).dim_mechanical;
	bool evolver_code = false;
	double omega_n_plus1;
	double consistency_order;

	//Check the whether this step should be first, second, or third order.
	for (int vertex_index = 0; vertex_index < dim_mechanical; vertex_index++){
		if (f_mech_step_before[0][vertex_index].x == 0.0){
			first_order_check++;
		}
		if (f_mech_step_before[1][vertex_index].x == 0.0){
			second_order_check++;
		}
	}

	if ((first_order_check != 0) && (second_order_check != 0)){
		//Perform a first-order step (if you don't have sufficient past information -> from the last two steps).
		consistency_order = 1.0;
		evolver_code = (*this).Integrate_IMEX_Adaptive(step_size, mech_abs_tolerance, params, current_time, consistency_order, errorfile);
	}
	else if ((first_order_check != 0) && (second_order_check == 0)){
		//Perform a second-order step (if you only have information from the last step).
		consistency_order = 2.0;
		evolver_code = (*this).Integrate_IMEX_Adaptive(step_size, mech_abs_tolerance, params, current_time, consistency_order, errorfile);
	}
	else{
		//Perform a third-order step (if you have information from the past two steps).
		consistency_order = 3.0;
		evolver_code = (*this).Integrate_IMEX_Adaptive(step_size, mech_abs_tolerance, params, current_time, consistency_order, errorfile);
	}
}

//Function called by Integrate_IMEX_Without_Protein for adaptive step size.
bool PCPTissue::Integrate_IMEX_Adaptive(double* step_size, const double error_abs_tolerance, void* params, double* t, double consistency_order, ofstream* errorfile){
	//dimension of mechanical systems
	int dimension_mechanical = vertices.size() + 1;
	//Store the current step size (make sure the total scale factor is not out of range).
	double current_step_size = *step_size;
	//Define the size of a half step.
	double half_step = current_step_size / 2.;
	//Initialize the error check for entering while loop.
	bool step_error_check = true;
	//These vectors are used for both half and full steps below as well as error tolerance and estimated error.
	vector<Vect> initial_condition;
	vector<Vect> y_error_tolerance;
	vector<Vect> error_vector;
	vector<Vect> full_step_solution;
	vector<Vect> second_half_step_solution;
	vector<Vect> first_half_step_solution;
	initial_condition.resize(dimension_mechanical);
	y_error_tolerance.resize(dimension_mechanical);
	error_vector.resize(dimension_mechanical);
	full_step_solution.resize(dimension_mechanical);
	second_half_step_solution.resize(dimension_mechanical);
	first_half_step_solution.resize(dimension_mechanical);

	//Store a variable for storing the information from the first half step.
	combinedparameters half_step_params = (*(combinedparameters*)params);
	//The lower bound of the scaling factor for the step size (for zero-stability)
	double lower_bound_scaling_factor = 0.5818;

	//Store initial condition.
	(*this).Store_Initial_Condition(initial_condition);

	while (step_error_check == true){
		//Reinitialize half_step_params if you are reentering the loop.
		half_step_params = (*(combinedparameters*)params);
		//Reinitialize half_step_size (upon entering loop)
		half_step = current_step_size / 2.;

		if (consistency_order == 1.0){
			//Take the first half step.
			(*this).Integrate_Euiler(half_step, errorfile);
		}
		else if (consistency_order == 2.0){
			//Take the first half step.
			(*this).Integrate_IMEX_Second_Order(half_step, (void*)&half_step_params, errorfile);
		}
		else if (consistency_order == 3.0){
			//Take the first half step.
			(*this).Integrate_IMEX_Third_Order(half_step, (void*)&half_step_params, errorfile);
		}
		//Store info from the current step into temporary params.
		store_evolve_info((void*)&half_step_params, half_step);
		//Store solution from the first half step.
		(*this).Store_New_Positions(first_half_step_solution);

		bool evolve_mechanical = true;
		//Update the params for the first half step (now including the information from that step).
		update_past_information_solver((void*)&half_step_params);
		//Store values from first half step as the initial condition for the next half step.
		(*this).Update_Old_Mechanical_Variables(first_half_step_solution);

		if (consistency_order == 1.0){
			//Take the second half step.
			(*this).Integrate_Euiler(half_step, errorfile);
		}
		else if (consistency_order == 2.0){
			//Take the second half step.
			(*this).Integrate_IMEX_Second_Order(half_step, (void*)&half_step_params, errorfile);
		}
		else if (consistency_order == 3.0){
			//Take the second half step.
			(*this).Integrate_IMEX_Third_Order(half_step, (void*)&half_step_params, errorfile);
		}
		//Store solution from the second half step.
		(*this).Store_New_Positions(second_half_step_solution);

		//restore initial conditions:
		(*this).Update_Old_Mechanical_Variables(initial_condition);

		if (consistency_order == 1.0){
			//Take a full step.
			(*this).Integrate_Euiler(current_step_size, errorfile);
		}
		else if (consistency_order == 2.0){
			//Take a full step.
			(*this).Integrate_IMEX_Second_Order(current_step_size, params, errorfile);
		}
		else if (consistency_order == 3.0){
			//Take a full step.
			(*this).Integrate_IMEX_Third_Order(current_step_size, params, errorfile);
		}
		//Store info from the current step.
		store_evolve_info(params, current_step_size);
		//Store full step solution.
		(*this).Store_New_Positions(full_step_solution);

		//restore initial conditions:
		(*this).Update_Old_Mechanical_Variables(initial_condition);

		//Estimate error by subtracting two half steps from full step.
		error_vector.resize(dimension_mechanical);
		for (int index = 0; index != dimension_mechanical; index++){
			error_vector.at(index) = full_step_solution.at(index) - second_half_step_solution.at(index);
			error_vector.at(index) = error_vector.at(index).Abs();
		}

		//Check whether step size is too large. (if it is too large, this function will return step_error_check != 0)
		step_error_check = Double_Step_Size_Control(error_vector, error_abs_tolerance, &current_step_size, consistency_order);

		//if the step size was scaled down, reinitialize initial_step_size for reentering while loop.
		if (step_error_check != false){
			*step_size = current_step_size;
		}
	}

	//Note: it is necessary to calculate 2 * half_step because this is the step size used in the evolution above (current_step_size could have been scaled up after evolution).
	double used_step_size = 2 * half_step;

	//Increment the time by the utilized step size (once the step is successful).
	if (step_error_check == false){
		//It is important to increment the time by two half steps (because the variable current_step_size may have been scaled up).
		*t += used_step_size;
		//Update information
		bool evolve_mechanical = true;
		update_past_information_solver(params);
	}

	//Store the current_step_size in step size if it exceeds the step size.
	if (current_step_size > *step_size) *step_size = current_step_size;

	return step_error_check;
}

// Do evolution by integrating forward the system of ODE from t=0 to t=timeInterval.
// Stepsize is the inital step size of the integrator. The function will modify both stepsize and currentTime for use by subsequent call to Evolve(), useful for breaking integration into many smaller intervals. 
// PCP Generation, Type Specification, Normalization occur within this function.
void PCPTissue::Evolve(double interval, double* prot_stepSize, double* mech_stepSize, ofstream* errorfile){
	//find the end point of this evolve->Closest f
	double endtime = floor((CURRENTTIME + interval) / interval) * interval;
	bool flag0, flag1, flag2, flag3;
	bool reinitialize_mechanical = true;
	bool reinitialize_protein = true;
	//Changes by Hayden for PCP: protein and vertex positions are evolved separately (keep track of two times). After evolving protein concentrations, vertex positions are evolved to the corresponding time point.
	int dim_protein;												//The dimension of the mechanical system and the protein system.
	int status;
	vector<double> y_protein;										//Current log variable values.
	//Store parameters for evolution for mechanical and protein in params.
	combinedparameters params;										//parameters for protein evolution.
	//Variable used to calculate protein free energy
	double protein_free_energy;
	//Keep track of total steps taken.
	int total_steps_taken = 0;
	//Specify which order you would like to use: 1, 2, or 3.
	int method_order = 2;


	while (CURRENTTIME < endtime){

		//If topology has changed (T1, apoptosis, mitosis, or anything that changed the dimension of the system N and/or changed the ordering of vertices in list), then we need to reinitialize ODE objects
		if (reinitialize_mechanical){
			for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
				(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
			}
			params.vert_id2index = CreateVertexId2IndexMap();		//Get vertex & cell id to index map
			params.cell_id2index = CreateCellId2IndexMap();
			int dim_mechanical = vertices.size() + 1;
			params.dim_mechanical = dim_mechanical;
			//higher order methods:
			double init_val = 0.0;
			vector<vector<Vect>> vec_vect(3, vector<Vect>(params.dim_mechanical, Vect(init_val, init_val)));
			if (reinitialize_IMEX_solver){
				//Initialize the vectors for storing information from past steps.
				f_mech_step_before = vec_vect;
				past_mech_state = vec_vect;
				past_step_sizes.resize(3);
				*mech_stepSize = INITIALSTEPSIZE;
			}
			reinitialize_mechanical = false;
		} //End of initialization

		//If topology has changed (T1, apoptosis, mitosis, or anything that changed the dimension of the system N and/or changed the ordering of vertices in list), then we need to reinitialize ODE objects
		if (reinitialize_protein){
			//Changes by Hayden for PCP: Define a variable to store the number of edges in differentiated cells in the tissue.
			int M = (*this).Find_Number_Differentiated_Edges();
			params.dimension_protein = 2 * M; //Number of dynamical variables associated with protein concentrations.
			dim_protein = params.dimension_protein;
			y_protein.resize(dim_protein);

			//Initialize protein vector with current log variable values.
			(*this).HelperEvolve_Initialize_Protein_Vector(y_protein);

			//higher order methods:
			double init_val = 0.0;
			//The convention for params.y_step_before, params.f_step_before, params.g_step_before, params.past_step_sizes
			//before a step is taken, [0] store info from two steps before, [1] stores info from one step before.
			//after a step is taken, [0] store info from two steps before, [1] stores info from one step before, [2] stores information from the step that just happened.
			//then, within the imex_solver_adaptive, [1]->[0], [2]->[1] such that the information is updated for the next step.
			vector<vector<double>> vec(3, vector<double>(dim_protein, init_val));
			vector<vector<Vect>> vec_vect(3, vector<Vect>(params.dim_mechanical, Vect(init_val, init_val)));

			//TEMP FOR DEBUGGING: only do first order steps for the first two steps. (after that just use information from past steps even if those steps occurred in the last Evolve.
			if (reinitialize_IMEX_solver){
				//Initialize the vectors for storing information from past steps.
				y_step_before = vec;
				f_step_before = vec;
				g_step_before = vec;
				past_step_sizes.resize(3);
				*prot_stepSize = INITIALSTEPSIZE;
			}
			reinitialize_protein = false;		//Finally set reinitialize protein to false
			reinitialize_IMEX_solver = false;		//Set flag for topological changes between instances of evolve to false.
		}
		/*******Check if this code accomplishes my main goal: fixing the weird blowing up/setting very very small protein amounts to zero so that they don't cause numerical instability*******/

		////check protein free energy (check that it is monotonically decreasing).
		//protein_free_energy = (*this).CalculatePCPFreeEnergy((void*)&params);

		//STEP 2: Integrate 
		//Loop through the vertices and record the current position in rold 
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++) { (*(*vit)).rold = (*(*vit)).r; }
		LXold = LX;				//Record the old box size
		LYold = LY;

		//Evolve protein and vertex positions with IMEX Solver.
		(*this).HelperEvolve_IMEX(y_protein, prot_stepSize, mech_stepSize, error_tolerance_absolute_protein, error_tolerance_absolute_mechanical, &CURRENTTIME, (void*)&params, method_order, errorfile);

		//Use this for only evolving mechanical positions.
		//(*this).Integrate_IMEX_Without_Protein(mech_stepSize, error_tolerance_absolute_mechanical, &CURRENTTIME, (void*)&params, errorfile);

		//Update protein variables with new protein variables.
		(*this).HelperEvolve_Update_Protein_Variables(y_protein);

		if (evolve_mech_positions){
			//Step 3: Update
			LX = LXnew;					//Update the current box size
			LY = LYnew;
			(*this).HelperEvolve_UpdatePosition((void*)&params);  //Update the vertex positions (this automaticaly updates edge flags if a vertex has crossed the boundary)
		}

		total_steps_taken++;					//Keep track of total steps taken.

		//STEP 4: look for topological changes and set flags if any have occured 
		//No apoptosis in my simulations
		flag1 = 0;
		//Somehow let the code know if either of these occurred on differentiated edges.
		if (evolve_mech_positions){
			flag2 = (*this).T1(errorfile);
			flag3 = (*this).ResolveVertices(errorfile);
		}
		else{
			flag2 = 0;
			flag3 = 0;
		}
		if (flag1 || flag2 || flag3){
			cout << "reinitialize protein" << endl;
			reinitialize_protein = true;
			reinitialize_mechanical = true;		//initialize mechanical and protein systems separately.
			reinitialize_IMEX_solver = true;	//use external variable to check whether we need to reinitialize_IMEX_solver (because may need to reinitialize upon reentering evolve).
		}
	}

	cout << total_steps_taken << endl;

	cout << "TISSUE EVOLVED TO TIME " << CURRENTTIME << "  LX = " << LX << "  LY = " << LY << endl;

}

//HELPER FUNCTIONS

//! Helper function for checking the cell input conforms to expectations.
void PCPTissue::HelperReadInputData_CheckTypes(){
	list<Cell*>::iterator cit;
	for (cit = cells.begin(); cit != cells.end(); cit++){	//look through all of the cells
		if ((*(*cit)).type != standard_cell_type && (*(*cit)).type != pcp_cell_type){
			cout << "Not all cell types are the standard cell type. Code will exit." << endl;
			exit(EXIT_FAILURE);
		}
	}
}

// Moves the vertices to their new positions and updates the approproate flags if they cross a boundary.
bool PCPTissue::HelperEvolve_UpdatePosition(void* params){
	double newx, newy;
	bool VertexHasCrossed = 0;
	int i = 0;

	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
		newx = (*(*vit)).rnew.x;
		newy = (*(*vit)).rnew.y;
		//look at the x component and update and change flags
		if (newx > LX){
			(*(*vit)).r.x = newx - LX;
			//Update the past state of the system based on periodic boundary condition.
			past_mech_state[0][i].x = past_mech_state[0][i].x - LX;
			past_mech_state[1][i].x = past_mech_state[1][i].x - LX;
			(*(*vit)).ChangeXEdgeFlags();
			VertexHasCrossed = true;
		}
		else if (newx < 0){
			(*(*vit)).r.x = LX + newx;
			//Update the past state of the system based on periodic boundary condition.
			past_mech_state[0][i].x = past_mech_state[0][i].x + LX;
			past_mech_state[1][i].x = past_mech_state[1][i].x + LX;
			(*(*vit)).ChangeXEdgeFlags();
			VertexHasCrossed = true;
		}
		else{
			(*(*vit)).r.x = newx;
		}
		//look at the y component and update and change flags
		if (newy > LY){
			(*(*vit)).r.y = newy - LY;
			//Update the past state of the system based on periodic boundary condition.
			past_mech_state[0][i].y = past_mech_state[0][i].y - LY;
			past_mech_state[1][i].y = past_mech_state[1][i].y - LY;
			(*(*vit)).ChangeYEdgeFlags();
			VertexHasCrossed = true;
		}
		else if (newy < 0){
			(*(*vit)).r.y = LY + newy;
			//Update the past state of the system based on periodic boundary condition.
			past_mech_state[0][i].y = past_mech_state[0][i].y + LY;
			past_mech_state[1][i].y = past_mech_state[1][i].y + LY;
			(*(*vit)).ChangeYEdgeFlags();
			VertexHasCrossed = true;
		}
		else{
			(*(*vit)).r.y = newy;
		}
		i++;
	}
	return VertexHasCrossed;
}

//! Step size control mechanical evolution (using double step error estimation) -> uses parameters for scaling step size (based on IMEX zero-stability)
bool PCPTissue::Double_Step_Size_Control(const vector<Vect> &y_error, const double error_tolerance_absolute, double* initial_step_size, double consistency_order){
	//dimension of mechanical systems (vertices plus the boxsize)
	int dimension_mechanical_vect = vertices.size() + 1;
	int dimension_mechanical_double = 2 * dimension_mechanical_vect;
	//Initialize the flag that is output to show if satisying error tolerance.
	bool step_error_check = false;
	//Initialize the factor for scaling the step size.
	double scale_factor_double_step = 1.0;
	//Initialize a variable for storing the adjusted step size.
	double adjusted_step_size = *initial_step_size;
	//The order of the numerical solver (consistency_order)
	double safety_factor = .9;
	//perhaps tune this up.
	double upper_fractional_tolerance = 1.0;
	//possibly change this and the parameters
	double lower_fractional_tolerance = 0.1;
	//Check whether error exceeds tolerance.
	vector<double> excess_error;
	excess_error.resize(dimension_mechanical_double);
	//Check whether error is excessively below tolerance.
	vector<double> low_error;
	low_error.resize(dimension_mechanical_double);
	//Ratio of error to tolerance
	vector<double> error_ratio;
	error_ratio.resize(dimension_mechanical_double);
	//Store the max ratio of the error to the tolerance.
	double max_error_ratio = 0.0;
	//Check for excess error.
	double max_error_check = 0.0;

	//For stability, enforce bounds on the scale factor.
	double lower_bound_scaling_factor = 0.5818;
	double upper_bound_scaling_factor = 1.1892;


	int i = 0;
	//Subtract the tolerated from error from the observed error (for all components)
	for (int index = 0; index != dimension_mechanical_vect; index++){
		excess_error.at(2 * i) = y_error.at(index).x - upper_fractional_tolerance * error_tolerance_absolute;
		excess_error.at(2 * i + 1) = y_error.at(index).y - upper_fractional_tolerance * error_tolerance_absolute;
		//Store the ratio of error to tolerance.
		error_ratio.at(2 * i) = y_error.at(index).x / error_tolerance_absolute;
		error_ratio.at(2 * i + 1) = y_error.at(index).y / error_tolerance_absolute;
		i++;
	}

	//Store the maximum error ratio.
	max_error_ratio = *max_element(error_ratio.begin(), error_ratio.end());
	//Store the maximum error check.
	max_error_check = *max_element(excess_error.begin(), excess_error.end());

	//If any of the components violated the error tolerance as enforced, set a factor for scaling the step size.
	if (max_error_check > 0.0){
		scale_factor_double_step = safety_factor * pow(max_error_ratio, -1 / consistency_order);
	}

	i = 0;
	//Check if all of the error is below the low fractional tolerance (scale up step size)
	for (int index = 0; index != dimension_mechanical_vect; index++){
		low_error.at(2 * i) = y_error.at(index).x - lower_fractional_tolerance * error_tolerance_absolute;
		low_error.at(2 * i + 1) = y_error.at(index).y - lower_fractional_tolerance * error_tolerance_absolute;
		i++;
	}

	//If all of the errors are less than half of the tolerance, define a factor for upscaling the step size.
	max_error_check = *max_element(low_error.begin(), low_error.end());
	if (max_error_check < 0.0){
		scale_factor_double_step = safety_factor * pow(max_error_ratio, -1 / (consistency_order + 1));
	}

	//If the up scale factor is greater than the upper bound, snap the scale factor to the upper bound.
	if (scale_factor_double_step > upper_bound_scaling_factor){
		scale_factor_double_step = upper_bound_scaling_factor;
	}
	if (scale_factor_double_step < 1.00){
		step_error_check = true;
	}
	//If the scale factor is lower than the lower bound, snap the scale factor to the lower bound.
	if (scale_factor_double_step < lower_bound_scaling_factor){
		scale_factor_double_step = lower_bound_scaling_factor;
	}

	//Adjust the step size.
	adjusted_step_size = *initial_step_size * scale_factor_double_step;
	*initial_step_size = adjusted_step_size;

	//Return the flag for if the step size was too large.
	return step_error_check;
}

//It is not absolutely imperative at this point for PCPTissue to have a separate check stability, but I envision that in the future, we may want to check stability differently in the case of PCP.
//Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration. Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors retruns false if the vertex is unstable in that topology
// \note I HAVE INCLUDED THIS FUNCTION IN CASE WE INSERT DIFFERENT CRITERIA FOR STABILITY BASED ON PCP
double PCPTissue::HelperResolveVertices_CheckStability(Vertex* vert, Edge* edge1, Edge* edge2, Edge* edge3, Edge* edge4, Cell* cell1, Cell* cell2, Cell* cell3, Cell* cell4, Vect *Force){
	//First make new objects which are coppies of the objects in the linked lists we will minipulate these objecs and they will be destroyed when the function ends
	PCPCell* pcp_c1 = dynamic_cast<PCPCell*>(cell1);
	PCPCell* pcp_c2 = dynamic_cast<PCPCell*>(cell2);
	PCPCell* pcp_c3 = dynamic_cast<PCPCell*>(cell3);
	PCPCell* pcp_c4 = dynamic_cast<PCPCell*>(cell4);
	PCPCell c1 = (*pcp_c1);
	PCPCell c2 = (*pcp_c2);
	PCPCell c3 = (*pcp_c3);
	PCPCell c4 = (*pcp_c4);
	PCPEdge* pcp_edge1 = dynamic_cast<PCPEdge*>(edge1);
	PCPEdge* pcp_edge2 = dynamic_cast<PCPEdge*>(edge2);
	PCPEdge* pcp_edge3 = dynamic_cast<PCPEdge*>(edge3);
	PCPEdge* pcp_edge4 = dynamic_cast<PCPEdge*>(edge4);
	PCPEdge e1 = (*pcp_edge1);
	PCPEdge e2 = (*pcp_edge2);
	PCPEdge e3 = (*pcp_edge3);
	PCPEdge e4 = (*pcp_edge4);
	int stability;
	//Changes by Hayden for PCP.
	int largestID = (*max_element(edgeIDs.begin(), edgeIDs.end()));
	double a1, a2, a3, a4;

	//Make sure that these edges have distinct ids
	Vertex v12; //the vertex on edges 1 and 2
	Vertex v34; //the vertex on edges 3 and 4

	//Note: important change, initilizing all at once (also assigning an actual id rather than having garbage handed to it).
	//Make the new edge and vertices resulting from the split into 2 3fold vertices
	//Please note: enew has no polarity protein.
	PCPEdge enew = PCPEdge::PCPEdge(largestID + 1, 0, 0, &v12, &v34, &c4, &c2);

	//Now update the info on the coppies of the objects so they are 2 3fold vertices

	//For the new vertices
	v12.r = (*vert).r;							//put the two new verts where the old fourfold was
	v12.edges.push_back(&e1);
	v12.edges.push_back(&e2);
	v12.edges.push_back(&enew);
	v34.r = (*vert).r;							//put the two new verts where the old fourfold was
	v34.edges.push_back(&e3);
	v34.edges.push_back(&e4);
	v34.edges.push_back(&enew);

	//For the old edges
	e1.UpdateVertex(vert, &v12);
	e2.UpdateVertex(vert, &v12);
	e3.UpdateVertex(vert, &v34);
	e4.UpdateVertex(vert, &v34);

	//Now we have at the position of the fourfold vertex two threefold vertices in which everyone knows all of their neighbors correclty, we can now calculate the forces based on what the tensiions and pressures are in this configuration

	//pressure differences across cells
	a1 = c1.GetArea();
	a2 = c2.GetArea();
	a3 = c3.GetArea();
	a4 = c4.GetArea();
	double p2 = (c2.GetPressure(a2) - c1.GetPressure(a1)); 			//pressure difference acorss edge2 between cell 1 and cell 2 taken CCW
	double p3 = (c3.GetPressure(a3) - c2.GetPressure(a2));
	double p4 = (c4.GetPressure(a4) - c3.GetPressure(a3));
	double p1 = (c1.GetPressure(a1) - c4.GetPressure(a4));

	//getting z X edge (both the diredction and magnitude are correct
	Vect zXe1 = (e1).GetVectorDifference(&v12, v12.GetOtherVertex(&e1)).Rotate(PI / 2.0);
	Vect zXe2 = (e2).GetVectorDifference(&v12, v12.GetOtherVertex(&e2)).Rotate(PI / 2.0);
	Vect zXe3 = (e3).GetVectorDifference(&v34, v34.GetOtherVertex(&e3)).Rotate(PI / 2.0);
	Vect zXe4 = (e4).GetVectorDifference(&v34, v34.GetOtherVertex(&e4)).Rotate(PI / 2.0);

	//Get the forces
	Vect ForceT1 = (e1).GetTension_ForceOnVertex(&v12);
	Vect ForceT2 = (e2).GetTension_ForceOnVertex(&v12);
	Vect ForceT3 = (e3).GetTension_ForceOnVertex(&v34);
	Vect ForceT4 = (e4).GetTension_ForceOnVertex(&v34);
	Vect ForceP1 = zXe1 * (p1 * 0.5);
	Vect ForceP2 = zXe2 * (p2 * 0.5);
	Vect ForceP3 = zXe3 * (p3 * 0.5);
	Vect ForceP4 = zXe4 * (p4 * 0.5);
	*Force = ((ForceT1 + ForceT2 + ForceP1 + ForceP2) - (ForceT3 + ForceT4 + ForceP3 + ForceP4))* 0.5;

	Vect ForceNewEdge;
	//c4 -> cell1 of enew
	//c2 -> cell2 of enew
	//Check which edges belong to c4 and which belong to c2.
	//e1 and e2 are paired, and e3 and e4 are paired.

	//Calculate TENSION ON VIRTUAL EDGE.
	if (((c4 == *e1.cell1) || (c4 == *e1.cell2)) && ((c4 == *e3.cell1) || (c4 == *e3.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e3, &e2, &e4);
	}
	else if (((c4 == *e1.cell1) || (c4 == *e1.cell2)) && ((c4 == *e4.cell1) || (c4 == *e4.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e4, &e2, &e3);
	}
	else if (((c2 == *e1.cell1) || (c2 == *e1.cell2)) && ((c2 == *e3.cell1) || (c2 == *e3.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e3, &e2, &e4);
	}
	else if (((c2 == *e1.cell1) || (c2 == *e1.cell2)) && ((c2 == *e4.cell1) || (c2 == *e4.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e4, &e2, &e3);
	}

	//Check the stability
	stability = (*Force).Magnitude() - ForceNewEdge.Magnitude();

	return stability;
}

//Changes by Hayden for PCP
//! This function is called by HelperResolveVertices_CheckStability and by HelperResolveVertices_ResolveFourfoldVertex. It redistributes proteins to the new edge according to the rules documented in the 2012 PLOS paper.
void PCPTissue::HelperResolveVertices_RedistributeProteinsNewEdge(PCPEdge* newedge, PCPEdge* edge1, PCPEdge* edge2, PCPEdge* edge3, PCPEdge* edge4){
	//Please note that I am using the convention explained in the comments for HelperResolveVertices_CheckStability and _UpdatePtrs:
	//Splits vertices so edges(e1, e2) and (e3, e4) are paired and cells c2 and c4 are neighbors
	//The new edge is between c2 and c4
	//e1 is between c1 and c4, e2 is between c1 and c2, e3 is between c2 and c3, e4 is between c3 and c4.

	//Make sure that you are accessing the correct elements for edge1, edge2, edge3, edge4 protein polarity lists
	vector<int> which_cell;
	//which_cell will store for e1, e2, e3, e4 (respectively) whether cell1 or cell2 refers to the cell shared with newedge
	int number_of_edges = 4;
	which_cell.resize(number_of_edges);
	int edge_index = 0;
	Cell* c4 = (*newedge).cell1;
	Cell* c2 = (*newedge).cell2;
	double currentedge_length;


	//which_cell will store for e1, e2, e3, e4 (respectively) whether cell1 or cell2 refers to the cell shared with newedge
	if ((*edge1).GetThisCellNumber(c4) == 1){
		//storing 0 if correct cell is cell 1-> this is useful when later accessing the planarpolarityconclist
		which_cell[0] = 0;
	}
	else {
		//storing 1 if correct cell is cell 2-> this is useful when later accessing the planarpolarityconclist
		which_cell[0] = 1;
	}

	if ((*edge2).GetThisCellNumber(c2) == 1){
		//storing 0 if correct cell is cell 1-> this is useful when later accessing the planarpolarityconclist
		which_cell[1] = 0;
	}
	else {
		//storing 1 if correct cell is cell 2-> this is useful when later accessing the planarpolarityconclist
		which_cell[1] = 1;
	}

	if ((*edge3).GetThisCellNumber(c2) == 1){
		//storing 0 if correct cell is cell 1-> this is useful when later accessing the planarpolarityconclist
		which_cell[2] = 0;
	}
	else {
		//storing 1 if correct cell is cell 2-> this is useful when later accessing the planarpolarityconclist
		which_cell[2] = 1;
	}

	if ((*edge4).GetThisCellNumber(c4) == 1){
		//storing 0 if correct cell is cell 1-> this is useful when later accessing the planarpolarityconclist
		which_cell[3] = 0;
	}
	else {
		//storing 1 if correct cell is cell 2-> this is useful when later accessing the planarpolarityconclist
		which_cell[3] = 1;
	}

	//Assign polarity to the new edge.
	//(*newedge).cell1 corresponds to c4, (*newedge).cell2 corresponds to c2
	//Note the rules for which_cell==0 or which_cell==1
	//NOTE: these calculations are explicitly in terms of the total protein, not protein concentration. (CONCENTRATION IS JUST PROTEIN AMOUNT ON EDGE/EDGE LENGTH)
	(*newedge).planarpolarity[0] = NEWEDGELENGTH * (kT1 * min((*edge1).planarpolarity[2 * which_cell[0]] / (*edge1).GetLength(), (*edge4).planarpolarity[2 * which_cell[3]] / (*edge4).GetLength()));

	(*newedge).planarpolarity[1] = NEWEDGELENGTH * (kT1 * min((*edge1).planarpolarity[2 * which_cell[0] + 1] / (*edge1).GetLength(), (*edge4).planarpolarity[2 * which_cell[3] + 1] / (*edge4).GetLength()));

	(*newedge).planarpolarity[2] = NEWEDGELENGTH * (kT1 * min((*edge2).planarpolarity[2 * which_cell[1]] / (*edge2).GetLength(), (*edge3).planarpolarity[2 * which_cell[2]] / (*edge3).GetLength()));

	(*newedge).planarpolarity[3] = NEWEDGELENGTH * (kT1 * min((*edge2).planarpolarity[2 * which_cell[1] + 1] / (*edge2).GetLength(), (*edge3).planarpolarity[2 * which_cell[2] + 1] / (*edge3).GetLength()));

	//Change the polarity of the edges 1-4 (Make sure that total protein is conserve). That means the sum of delta_protein_conc*length over changed edges be zero.
	//Follow conventions from LUBENSKY PLOS 2012 for redistributing protein to the new edge.
	currentedge_length = (*edge1).GetLength();
	(*edge1).planarpolarity[2 * which_cell[0]] = (*edge1).planarpolarity[2 * which_cell[0]] - (*newedge).planarpolarity[0] / 2;
	(*edge1).planarpolarity[2 * which_cell[0] + 1] = (*edge1).planarpolarity[2 * which_cell[0] + 1] - (*newedge).planarpolarity[1] / 2;

	currentedge_length = (*edge4).GetLength();
	(*edge4).planarpolarity[2 * which_cell[3]] = (*edge4).planarpolarity[2 * which_cell[3]] - (*newedge).planarpolarity[0] / 2;
	(*edge4).planarpolarity[2 * which_cell[3] + 1] = (*edge4).planarpolarity[2 * which_cell[3] + 1] - (*newedge).planarpolarity[1] / 2;

	currentedge_length = (*edge2).GetLength();
	(*edge2).planarpolarity[2 * which_cell[1]] = (*edge2).planarpolarity[2 * which_cell[1]] - (*newedge).planarpolarity[2] / 2;
	(*edge2).planarpolarity[2 * which_cell[1] + 1] = (*edge2).planarpolarity[2 * which_cell[1] + 1] - (*newedge).planarpolarity[3] / 2;

	currentedge_length = (*edge3).GetLength();
	(*edge3).planarpolarity[2 * which_cell[2]] = (*edge3).planarpolarity[2 * which_cell[2]] - (*newedge).planarpolarity[2] / 2;
	(*edge3).planarpolarity[2 * which_cell[2] + 1] = (*edge3).planarpolarity[2 * which_cell[2] + 1] - (*newedge).planarpolarity[3] / 2;
}

//Changes by Hayden for PCP
// \note make sure that you do a thorough debugging of this code.
//! This function is called by HelperT1RemoveEdge. It redistributes proteins from the edge to be deleted to nearby edges as described in the PLOS paper.
//! Consult Figure S5B of LUBENSKY PLOS 2012 as a reference for what I am doing.
void PCPTissue::HelperT1_RedistributeProteinsDeletedEdge(PCPEdge* deledge, vector<PCPEdge*> newedgelist){
	double new_edge_length;									//Store the edge currently under consideration to receive protein from the deleted edge.
	double current_protein_amount;								//Protein amount before redistribution
	int number_of_cell_pointers = newedgelist.size() * 2;	//each edge has two cell pointers cell1 and cell2
	vector<Cell*> connected_cells;								//for each edge store pointer to cell1 and cell2.
	connected_cells.resize(number_of_cell_pointers);
	int edge_index = 0;
	//Have vector for both pairs of edges.
	vector<PCPEdge*> paired_edges1;
	vector<PCPEdge*> paired_edges2;
	//Use a vector to store for the other edges whether they contact cell1 or cell2.
	vector<int> which_cell;
	which_cell.resize(newedgelist.size());
	//Store the cells on either side of the edge to be deleted.
	Cell* selected_cell_pointer1 = (*deledge).cell1;
	Cell* selected_cell_pointer2 = (*deledge).cell2;

	/***Temp for debugging***/
	//Note the code has already checked for fivefold vertices and exited in that case. So here I assume that the newedgelist is of length 4. I check here though as a debug measure.
	assert((newedgelist.size() == 4) && "ERROR Tissue::Helper_Redistribute_Proteins - There should be four edges to be joined (in newedgelist)! Tell HN.");
	/***Temp for debugging***/

	//Store cell1 and cell2 for every new edge.
	for (edge_index = 0; edge_index != newedgelist.size(); edge_index++){
		connected_cells[2 * edge_index] = (*newedgelist[edge_index]).cell1;
		connected_cells[2 * edge_index + 1] = (*newedgelist[edge_index]).cell2;
	}

	//Check for which edges are on the same side of the edge to be deleted. (These edges will share protein that will be redistributed.)
	for (edge_index = 0; edge_index != newedgelist.size(); edge_index++){
		if (((*connected_cells[2 * edge_index]) == (*selected_cell_pointer1)) || ((*connected_cells[2 * edge_index + 1]) == (*selected_cell_pointer1))){
			paired_edges1.push_back(newedgelist[edge_index]);
		}
		else {
			paired_edges2.push_back(newedgelist[edge_index]);
		}
	}

	/***Temp for debugging***/
	//Note the code has already checked for fivefold vertices and exited in that case. So here I assume that the newedgelist is of length 4. I check here though as a debug measure.
	assert((paired_edges1.size() == 2) && "ERROR Tissue::Helper_Redistribute_Proteins - There should be two edges in paired_edges! Tell HN.");
	assert((paired_edges2.size() == 2) && "ERROR Tissue::Helper_Redistribute_Proteins - There should be two edges in paired_edges! Tell HN.");
	/***Temp for debugging***/

	//This is a common issue. Since which cell is stored in cell1 or cell2 for a given edge is arbitrary, you must make sure that the proteins that are redistributed stay within the same cell.
	for (edge_index = 0; edge_index != paired_edges1.size(); edge_index++){
		if ((*paired_edges1[edge_index]).GetThisCellNumber(selected_cell_pointer1) == 1){
			//storing 0 if correct cell is cell 1-> this is useful when later accessing the planarpolarityconclist
			which_cell[edge_index] = 0;
		}
		else {
			//storing 1 if correct cell is cell 2-> this is useful when later accessing the planarpolarityconclist
			which_cell[edge_index] = 1;
		}
	}
	//This is a common issue. Since which cell is stored in cell1 or cell2 for a given edge is arbitrary, you must make sure that the proteins that are redistributed stay within the same cell.
	for (edge_index = 0; edge_index != paired_edges1.size(); edge_index++){
		if ((*paired_edges2[edge_index]).GetThisCellNumber(selected_cell_pointer2) == 1){
			//storing 0 if correct cell is cell 1-> this is useful when later accessing the planarpolarityconclist
			which_cell[edge_index + paired_edges1.size()] = 0;
		}
		else {
			//storing 1 if correct cell is cell 2-> this is useful when later accessing the planarpolarityconclist
			which_cell[edge_index + paired_edges1.size()] = 1;
		}
	}

	//now use this to redistribute the proteins
	//remember the following conventions for planarpolarityconc (element 0 refers to protein A cell1, element 1 to protein B cell1, element 2 to protein A cell2, element 3 to protein B cell2).
	//NOTE: this calculation is explicitly in terms of total protein on an edge and not the concentration on an edge.
	//Consult Figure S5B of LUBENSKY PLOS 2012 as a reference for what I am doing.

	//Redistribute the proteins for cell1, protein type A of *deledge
	current_protein_amount = (*deledge).planarpolarity[0];
	//now loop through paired_edges1
	for (edge_index = 0; edge_index != paired_edges1.size(); edge_index++){
		//Store the length of the edge to which protein is being redistributed.
		new_edge_length = (*paired_edges1[edge_index]).GetLength();
		//Check which element of protein concentration you should access and change.
		(*paired_edges1[edge_index]).planarpolarity[2 * which_cell[edge_index]] = (*paired_edges1[edge_index]).planarpolarity[2 * which_cell[edge_index]] + current_protein_amount / 2;
	}

	//Redistribute the proteins for cell1, protein type B of *deledge
	current_protein_amount = (*deledge).planarpolarity[1];
	//now loop through paired_edges1
	for (edge_index = 0; edge_index != paired_edges1.size(); edge_index++){
		//Store the length of the edge to which protein is being redistributed.
		new_edge_length = (*paired_edges1[edge_index]).GetLength();
		(*paired_edges1[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = (*paired_edges1[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] + current_protein_amount / 2;
	}

	//Redistribute the proteins for cell2, protein type A of *deledge
	current_protein_amount = (*deledge).planarpolarity[2];
	//now loop through paired_edges2
	for (edge_index = 0; edge_index != paired_edges2.size(); edge_index++){
		new_edge_length = (*paired_edges2[edge_index]).GetLength();
		(*paired_edges2[edge_index]).planarpolarity[2 * which_cell[edge_index + paired_edges1.size()]] = (*paired_edges2[edge_index]).planarpolarity[2 * which_cell[edge_index + paired_edges1.size()]] + current_protein_amount / 2;
	}

	//Redistribute the proteins for cell2, protein type B of *deledge
	current_protein_amount = (*deledge).planarpolarity[3];
	//now loop through paired_edges2
	for (edge_index = 0; edge_index != paired_edges2.size(); edge_index++){
		new_edge_length = (*paired_edges2[edge_index]).GetLength();
		(*paired_edges2[edge_index]).planarpolarity[2 * which_cell[edge_index + paired_edges1.size()] + 1] = (*paired_edges2[edge_index]).planarpolarity[2 * which_cell[edge_index + paired_edges1.size()] + 1] + current_protein_amount / 2;
	}

	//Remove all protein from the edge to be deleted. (This is not absolutely necessary->it is about to be deleted anyway).
	for (int i = 0; i != 4; i++){
		(*deledge).planarpolarity[i] = 0.0;
	}
}

// Writes data into file specified by 'datafile'. If flag 'includeExtraInfo' is true, then additional information (such as tension, cell area, etc) is included.
// The extra information is non-essential in reconstructing the simulation state, and is just provided for convenience (e.g. for display in MATLAB).
void PCPTissue::HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo){
	time_t rawtime;					//to deal with writing the date out
	struct tm * timeinfo;
	(*datafile).precision(15);
	(*datafile).setf(ios::scientific);
	(*datafile).setf(ios::showpoint);

	//First write the info about what you just ran and the time you ran it
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;
	*datafile << "Date:\t" << asctime(timeinfo) << endl;

	//Now write the global variables you used
	*datafile << "Global:" << endl;
	*datafile << "SEED\t" << SEED << endl;
	*datafile << "LX\t" << LX << endl;
	*datafile << "LY\t" << LY << endl;
	*datafile << "MINAREA\t" << MINAREA << endl;
	*datafile << "MINLENGTH\t" << MINLENGTH << endl;
	*datafile << "NEWEDGELENGTH\t" << NEWEDGELENGTH << endl;
	*datafile << "K\t" << K << endl;
	*datafile << "GAMMA\t" << GAMMA << endl;
	*datafile << "TIMEINTERVAL\t" << TIMEINTERVAL << endl;
	*datafile << "TOTALTIME\t" << TOTALTIME << endl;
	*datafile << "CURRENTTIME\t" << CURRENTTIME << endl;
	*datafile << "MAX_STEP_SIZE\t" << max_step_size << endl;
	*datafile << "INITIALSTEPSIZE\t" << INITIALSTEPSIZE << endl;
	*datafile << "error_tolerance_absolute_mechanical\t" << error_tolerance_absolute_mechanical << endl;
	*datafile << "error_tolerance_absolute_protein\t" << error_tolerance_absolute_protein << endl;
	*datafile << "INCREMENT_FOR_RENORMALIZATION\t" << increment_for_renormalization << endl;
	*datafile << "ALPHAP\t" << alphaP << endl;
	*datafile << "ALPHAM\t" << alphaM << endl;
	*datafile << "TISSUE_ALPHA\t" << TISSUE_ALPHA << endl;
	*datafile << "TISSUE_EXT_STRESS_X\t" << TISSUE_EXT_STRESS_X << endl;
	*datafile << "TISSUE_EXT_STRESS_Y\t" << TISSUE_EXT_STRESS_Y << endl;
	*datafile << "GAMMA0\t" << GAMMA0 << endl;
	*datafile << "GAMMAC\t" << GAMMAC << endl;
	*datafile << "EPSAB\t" << epsab << endl;
	*datafile << "EPSNAUGHT\t" << epsnaught << endl;
	*datafile << "EPSONE\t" << epsone << endl;
	*datafile << "EPSTWO\t" << epstwo << endl;
	*datafile << "TE\t" << te << endl;
	*datafile << "JONE\t" << jone << endl;
	*datafile << "JTWO\t" << jtwo << endl;
	*datafile << "KT1\t" << kT1 << endl;
	*datafile << "evolve_mech_positions\t" << evolve_mech_positions << endl;
	*datafile << "simulation_type\t" << simulation_type << endl;
	*datafile << "delta_protein_amount\t" << delta_protein_amount << endl;
	*datafile << "max_concentration\t" << max_concentration << endl;
	*datafile << "dirichlet_noise_param\t" << dirichlet_noise_param << endl;
	//******ADD NEW OUTPUT HERE******//
	*datafile << "endGlobal" << endl;

	// Save list of ordered vertices for evaluating forces for output below.
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
		(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
	}

	if (!includeExtraInfo)											// if flag is false, include only essential information to capture the state of the simulation
	{
		//print the vertex info
		*datafile << "Vertices" << endl;
		*datafile << "id\tx\ty" << endl;
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){			//loop through all of the vertices in the tissue
			*datafile << (*(*vit)).id << "\t" << (*(*vit)).r.x << "\t" << (*(*vit)).r.y << endl;						//write it's info to the file
		}

		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			*datafile << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << endl;			//write its info to the file
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;
		*datafile << "id\tA0\tEdgeNumb\tEdgeList\tType" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++){				//make comma-separated list of edges
				ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
			}

			*datafile << (*(*cit)).id << "\t" << (*(*cit)).A0 << "\t" << (*(*cit)).edges.size() << "\t" << ssEdgeList.str() << "\t" << (*(*cit)).type << endl;
		}
	}

	else															// if flag is true, include only essential information to capture the state of the simulation
	{
		//print the vertex info with additional information
		*datafile << "Vertices" << endl;
		*datafile << "id\tx\ty\t_Ftension_x\t_Ftension_y\t_Fpressure_x\t_Fpressure_y\t_EdgesNumber" << endl;
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){			//loop through all of the vertices in the tissue
			Vect FT = (*(*vit)).GetNetForce_Tension();
			Vect FP = (*(*vit)).GetNetForce_Pressure();
			*datafile << (*(*vit)).id << "\t" << (*(*vit)).r.x << "\t" << (*(*vit)).r.y << "\t" << FT.x << "\t" << FT.y << "\t" << FP.x << "\t" << FP.y << "\t" << (*(*vit)).edges.size() << endl;						//write it's info to the file
		}

		//Changes by Hayden for PCP: output PCP concentrations.
		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag\t_length\t_Tension_x\t_Tension_y\t_cell1\t_cell2\t_cell1AProtein\t_cell1BProtein\t_cell2AProtein\t_cell2BProtein" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			Vect FT = (*(*eit)).GetTension_ForceOnVertex(v1);
			double L = (*(*eit)).GetLength();
			PCPEdge* pcpeit = dynamic_cast<PCPEdge*>(*eit);
			vector<double> planar_polarity_protein = (*pcpeit).planarpolarity;
			*datafile << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*pcpeit).yflag << "\t" << L << "\t" << FT.x << "\t" << FT.y
				<< "\t" << (*(*pcpeit).cell1).id << "\t" << (*(*pcpeit).cell2).id << "\t" << planar_polarity_protein[0] << "\t" << planar_polarity_protein[1] << "\t" << planar_polarity_protein[2] << "\t" << planar_polarity_protein[3] << endl;										//write its info to the file
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;
		*datafile << "id\tA0\tEdgeNumb\tType\t_xflag\t_yflag\t_Area\t_Pressure\t_Center_x\t_Center_y\tEdgeList\t_CWVertices\t_OrderedFlagsX\t_OrderedFlagsY" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++){				//make comma-separated list of edges
				ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
			}

			vector<Vertex*> CWVertices = (*(*cit)).GetClockWiseVertices();
			stringstream ssCWVertices;
			ssCWVertices.clear();
			for (i = 0; i < CWVertices.size(); i++){			//make comma-separated list of vertices ordered clockwise
				ssCWVertices << (*CWVertices[i]).id << ",";
			}

			//Changes by Hayden for PCP: Outputting this is useful for calculating shape polarization vectors in MATLAB.
			vector<vector<int>> CWVertexFlags = (*(*cit)).HelperOrderFlags(CWVertices);
			stringstream ssCWVerticesFlagsX;
			ssCWVerticesFlagsX.clear();
			stringstream ssCWVerticesFlagsY;
			ssCWVerticesFlagsY.clear();
			for (i = 0; i < CWVertices.size(); i++){			//make comma-separated list of vertices ordered clockwise
				ssCWVerticesFlagsX << CWVertexFlags[0][i] << ",";
			}
			for (i = 0; i < CWVertices.size(); i++){			//make comma-separated list of vertices ordered clockwise
				ssCWVerticesFlagsY << CWVertexFlags[1][i] << ",";
			}

			Vect centroid = (*(*cit)).GetCenter();
			*datafile << (*(*cit)).id << "\t" << (*(*cit)).A0 << "\t" << (*(*cit)).edges.size() << "\t" << (*(*cit)).type << "\t" << (*(*cit)).GetXFlag() << "\t" << (*(*cit)).GetYFlag() << "\t"
				<< (*(*cit)).GetArea() << "\t" << (*(*cit)).GetPressure( (*(*cit)).GetArea() ) << "\t" << centroid.x << "\t" << centroid.y << "\t" << ssEdgeList.str() << "\t" << ssCWVertices.str() << "\t" << ssCWVerticesFlagsX.str() << "\t" << ssCWVerticesFlagsY.str() << endl;
		}
	}


	*datafile << "End" << endl << endl;
}

//! This is an umbrella function for handling the transition from first-second-third order evolution of protein concentrations.
void PCPTissue::HelperEvolve_IMEX(vector<double> &y_protein, double* prot_stepSize, double* mech_stepSize, const double prot_abs_tolerance, const double mech_abs_tolerance, double* current_time, void* params, int method_order, ofstream* errorfile){
	//Check if the information from both past steps exist.
	int first_order_check = 0;
	int second_order_check = 0;
	int dim_protein = (*(combinedparameters*)params).dimension_protein;
	bool pcp_evolver_code = false;

	//Check the whether this step should be first, second, or third order.
	for (int edge_index = 0; edge_index < dim_protein; edge_index++){
		if (y_step_before[0][edge_index] == 0.0){
			first_order_check++;
		}
		if (y_step_before[1][edge_index] == 0.0){
			second_order_check++;
		}
	}

	//NOTE: method_order specifies which order should be used (1, 2, 3)

	//if some of the past step information if missing, and this is the first step since reinitialization
	if ((first_order_check != 0 && second_order_check != 0) || method_order == 1){
		//Perform a first-order step (if you don't have sufficient past information -> from the last two steps).
		cout << "first order";
		pcp_evolver_code = imex_solver_adaptive(y_protein, prot_stepSize, mech_stepSize, prot_abs_tolerance, mech_abs_tolerance, current_time, params, &imex_solver_PCP_first_order, errorfile);
	}
	else if ((first_order_check != 0 && second_order_check == 0) || method_order == 2){
		//Perform a second-order step (if you only have information from the last step).
		pcp_evolver_code = imex_solver_adaptive(y_protein, prot_stepSize, mech_stepSize, prot_abs_tolerance, mech_abs_tolerance, current_time, params, &imex_solver_PCP_second_order, errorfile);
	}
	else if ((first_order_check == 0) && (second_order_check == 0) && method_order == 3){
		//Perform a third-order step (if you have information from the past two steps).
		pcp_evolver_code = imex_solver_adaptive(y_protein, prot_stepSize, mech_stepSize, prot_abs_tolerance, mech_abs_tolerance, current_time, params, &imex_solver_PCP_third_order, errorfile);
	}
}

//! Initialize protein vector of log variables (for use within evolve).
void PCPTissue::HelperEvolve_Initialize_Protein_Vector(vector<double> &y_protein){
	//Changes by Hayden for PCP: store the values of polarity protein concentrations.
	//Assign type A.
	//Loop through in order to store the protein variables for protein A.
	//Store the state vector as xi=-log(li*ci)
	//Note the following convention: if the value is stored as 0.0, that means that the corresponding protein amount is 0.0. (You can't take the log of 0.0).
	int i = 0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//order edges 
			vector<PCPEdge*> edge_list_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_this_cell.size(); edge_index++){
				//Check whether the given cell corresponds to cell1 or cell2 for the edges.
				if ((*edge_list_this_cell[edge_index]).GetThisCellNumber(*cit) == 1){
					y_protein.at(i) = -log((*edge_list_this_cell[edge_index]).planarpolarity[0]);
				}
				else{
					y_protein.at(i) = -log((*edge_list_this_cell[edge_index]).planarpolarity[2]);
				}
				i++;
			}
		}
	}

	//Changes by Hayden for PCP: store the values of polarity protein concentrations.
	//Assign type B.
	//Loop through in order to store the protein variables for protein A.
	//Store the state vector as xi=-log(li*ci)
	//Note the following convention: if the value is stored as 0.0, that means that the corresponding protein amount is 0.0. (You can't take the log of 0.0).
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//order edges
			vector<PCPEdge*> edge_list_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_this_cell.size(); edge_index++){
				//Check whether the given cell corresponds to cell1 or cell2 for the edges.
				if ((*edge_list_this_cell[edge_index]).GetThisCellNumber(*cit) == 1){
					y_protein.at(i) = -log((*edge_list_this_cell[edge_index]).planarpolarity[1]);
				}
				else{
					y_protein.at(i) = -log((*edge_list_this_cell[edge_index]).planarpolarity[3]);
				}
				i++;
			}
		}
	}
}

//! Update protein variables after integration of the dynamics.
void PCPTissue::HelperEvolve_Update_Protein_Variables(const vector<double> &y_protein){
	//Changes by Hayden for PCP: store the new values of the polarity protein.
	//Assign type A.
	//Loop through in order to store the protein variables for protein A.
	//Convert the xi=-log(li*ci) to an actual protein amount.
	int i = 0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//order edges
			vector<PCPEdge*> edge_list_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_this_cell.size(); edge_index++){
				if ((*edge_list_this_cell[edge_index]).GetThisCellNumber(*cit) == 1){
					//To extract protein concentration, exp(-xi)=exp(log(li*ci))=li*ci.
					(*edge_list_this_cell[edge_index]).planarpolarity[0] = exp(-y_protein.at(i));
				}
				else{
					(*edge_list_this_cell[edge_index]).planarpolarity[2] = exp(-y_protein.at(i));
				}
				i++;
			}
		}
	}

	//Now, do the same for protein variables of type B.
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//order edges
			vector<PCPEdge*> edge_list_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_this_cell.size(); edge_index++){
				if ((*edge_list_this_cell[edge_index]).GetThisCellNumber(*cit) == 1){
					//To extract protein concentration, exp(-xi)=exp(log(li*ci))=li*ci.
					(*edge_list_this_cell[edge_index]).planarpolarity[1] = exp(-y_protein.at(i));
				}
				else{
					(*edge_list_this_cell[edge_index]).planarpolarity[3] = exp(-y_protein.at(i));
				}
				i++;
			}
		}
	}
}

//Called by Evolve to renormalize the protein amounts in all cells. Lagrange multiplier constrains normalization, but we must correct for numerical drift.
//Loops through all cells and calculates the protein with each cell (length*concentration). Enforces normalization of both types of proteins.
void PCPTissue::HelperEvolve_RenormalizePolarityProteins(){
	double sum_of_proteinA_amounts = 0.0;
	double sum_of_proteinB_amounts = 0.0;
	int edge_index;
	vector<PCPEdge*> edge_list_for_this_cell;
	vector<double> edge_proteinA_list;
	vector<double> edge_proteinB_list;
	vector<int> which_cell;

	list<Cell*>::iterator cit;
	for (cit = cells.begin(); cit != cells.end(); cit++){			//look through all of the cells
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		//if a cell is differentiated, enforce normalization of proteins of all types (otherwise, leave the cell alone)
		if ((*pcpcit).is_Cell_with_PCP()){
			sum_of_proteinA_amounts = 0.0;
			sum_of_proteinB_amounts = 0.0;
			edge_list_for_this_cell.resize((*(*cit)).edges.size());
			for (edge_index = 0; edge_index != (*(*cit)).edges.size(); edge_index++){
				edge_list_for_this_cell.at(edge_index) = dynamic_cast<PCPEdge*>((*(*cit)).edges.at(edge_index));
			}
			edge_proteinA_list.resize(edge_list_for_this_cell.size());
			edge_proteinB_list.resize(edge_list_for_this_cell.size());
			which_cell.resize(edge_list_for_this_cell.size());

			for (edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){

				//compare cell ids to see which cell is referred to by cell1 and cell2
				if ((*edge_list_for_this_cell[edge_index]).GetThisCellNumber(*cit) == 1){
					edge_proteinA_list[edge_index] = (*edge_list_for_this_cell[edge_index]).planarpolarity[0];
					edge_proteinB_list[edge_index] = (*edge_list_for_this_cell[edge_index]).planarpolarity[1];
					//Use which cell to store whether the current cellIter is pointing to cell1 or cell2 of the given edge (there is no rhyme or reason to which is cell1 and which is cell2)
					which_cell[edge_index] = 0;
				}

				//compare cell ids to see which cell is referred to by cell1 and cell2
				if ((*edge_list_for_this_cell[edge_index]).GetThisCellNumber(*cit) == 2){
					edge_proteinA_list[edge_index] = (*edge_list_for_this_cell[edge_index]).planarpolarity[2];
					edge_proteinB_list[edge_index] = (*edge_list_for_this_cell[edge_index]).planarpolarity[3];
					//Use which cell to store whether the current cellIter is pointing to cell1 or cell2 of the given edge (there is no rhyme or reason to which is cell1 and which is cell2)
					which_cell[edge_index] = 1;
				}

			}

			//Sum over total protein of both types.
			for (edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				sum_of_proteinA_amounts += edge_proteinA_list[edge_index];
				sum_of_proteinB_amounts += edge_proteinB_list[edge_index];
			}

			//now renormalize the concentrations stored in the local vector
			for (edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				edge_proteinA_list[edge_index] = edge_proteinA_list[edge_index] / sum_of_proteinA_amounts;
				edge_proteinB_list[edge_index] = edge_proteinB_list[edge_index] / sum_of_proteinB_amounts;
			}

			//now go back and assign these back to the actual values of protein concentrations held by the edges.
			for (edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				(*edge_list_for_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index]] = edge_proteinA_list[edge_index];
				(*edge_list_for_this_cell[edge_index]).planarpolarity[2 * which_cell[edge_index] + 1] = edge_proteinB_list[edge_index];
			}
		}
	}
}

//! This function is used to loop through the list of undifferentiated cells that are in contact with differentiated cells at a propagation time.
void PCPTissue::HelperEvolve_GeneratePCPProtein(bool delay){
	//Locate next cells to get PCP
	vector<PCPCell*> cells_for_differentiation = (*this).HelperEvolve_FindBoundaryCells();
	int cell_index = 0;
	//Count number of cells to get PCP
	int number_of_cells_for_differentiation = cells_for_differentiation.size();

	//Order cells y position.
	vector<PCPCell*> ordered_cells_for_differentiation = order_cells_position(cells_for_differentiation);

	//Loop through next cells to get PCP.
	for (cell_index = 0; cell_index != number_of_cells_for_differentiation; cell_index++){
		//put PCP in all cells in list right now.
		(*ordered_cells_for_differentiation[cell_index]).HelperEvolve_InitializePolarityProteins();
	}
}

//! This function is used to locate cells on the boundary between differentiated and undifferentiated cells.
vector<PCPCell*> PCPTissue::HelperEvolve_FindBoundaryCells(){
	vector<PCPCell*> boundary_cells;
	vector<PCPEdge*> edge_list_for_this_cell;
	int number_of_edges;
	int boundary_flag = 0;
	int edge_index = 0;

	//Check for cells that are at the boundary of the differentiated cells (cells of type standard_cell_type that border cells with PCP)
	list<Cell*>::iterator cit;
	for (cit = cells.begin(); cit != cells.end(); cit++){			//look through all of the cells
		boundary_flag = 0;											//Use a flag to identify the cells that are at the boundary of differentiated cells.
		number_of_edges = (*(*cit)).edges.size();					//Store the number of edges in a given cell.
		edge_list_for_this_cell.resize(number_of_edges);
		for (edge_index = 0; edge_index != (*(*cit)).edges.size(); edge_index++){
			edge_list_for_this_cell.at(edge_index) = dynamic_cast<PCPEdge*>((*(*cit)).edges.at(edge_index));
		}
		for (edge_index = 0; edge_index != number_of_edges; edge_index++){
			PCPCell* pcp_cell1 = dynamic_cast<PCPCell*>((*edge_list_for_this_cell[edge_index]).cell1);
			PCPCell* pcp_cell2 = dynamic_cast<PCPCell*>((*edge_list_for_this_cell[edge_index]).cell2);
			//If the for any of the edges in this cell, there is a cell with PCP on one side and a cell of type standard_cell_type on the other side, then increment the flag to indicate that this is the case.
			if ((*pcp_cell1).is_Cell_with_PCP() && (*pcp_cell2).type == standard_cell_type){
				boundary_flag++;
			}
			if ((*pcp_cell2).is_Cell_with_PCP() && (*pcp_cell1).type == standard_cell_type){
				boundary_flag++;
			}
		}
		//Zero signals an undifferentiated cell. Boundary flag != 0 indicates that the cell is on the boundary of differentiation.
		if ((boundary_flag != 0) && ((*(*cit)).type == standard_cell_type)){
			boundary_cells.push_back(dynamic_cast<PCPCell*>(&(*(*cit)))); //push back the pointer to the cell that is undifferentiated and at the boundary of differentiation. It will differentiate in the next round.
		}
	}

	//Return the list of cells to differentiate in this round.
	return boundary_cells;
}

vector<PCPCell*> PCPTissue::HelperEvolve_FindCells_Contact_with_Boundary(){
	vector<PCPCell*> boundary_cells;
	vector<PCPEdge*> edge_list_for_this_cell;
	int number_of_edges;
	int boundary_flag = 0;
	int edge_index = 0;
	PCPCell* pcp_cit;

	//Check for cells that are at the boundary of the differentiated cells (cells of type standard_cell_type that border cells with PCP)
	list<Cell*>::iterator cit;
	for (cit = cells.begin(); cit != cells.end(); cit++){			//look through all of the cells
		pcp_cit = dynamic_cast<PCPCell*>(*cit);

		boundary_flag = 0;											//Use a flag to identify the cells that are at the boundary of differentiated cells.
		number_of_edges = (*(*cit)).edges.size();					//Store the number of edges in a given cell.
		edge_list_for_this_cell.resize(number_of_edges);
		for (edge_index = 0; edge_index != (*(*cit)).edges.size(); edge_index++){
			edge_list_for_this_cell.at(edge_index) = dynamic_cast<PCPEdge*>((*(*cit)).edges.at(edge_index));
		}
		for (edge_index = 0; edge_index != number_of_edges; edge_index++){
			PCPCell* pcp_cell1 = dynamic_cast<PCPCell*>((*edge_list_for_this_cell[edge_index]).cell1);
			PCPCell* pcp_cell2 = dynamic_cast<PCPCell*>((*edge_list_for_this_cell[edge_index]).cell2);
			//If the for any of the edges in this cell, there is a cell with PCP on one side and a cell of type standard_cell_type on the other side, then increment the flag to indicate that this is the case.
			if ((*pcp_cell1).is_Cell_with_PCP() && (*pcp_cell2).type == standard_cell_type){
				boundary_flag++;
			}
			if ((*pcp_cell2).is_Cell_with_PCP() && (*pcp_cell1).type == standard_cell_type){
				boundary_flag++;
			}
		}
		//Zero signals an undifferentiated cell. Boundary flag != 0 indicates that the cell is on the boundary of differentiation.
		if ((boundary_flag != 0) && ((*pcp_cit).is_Cell_with_PCP())){
			boundary_cells.push_back(dynamic_cast<PCPCell*>(&(*(*cit)))); //push back the pointer to the cell that is undifferentiated and at the boundary of differentiation. It will differentiate in the next round.
		}
	}

	//Return the list of cells to differentiate in this round.
	return boundary_cells;
}

//Called by HelperResolveVertices_ResolveFourfoldVertex to update the pointers when the vertex needs to resolve, Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors 
Edge* PCPTissue::HelperResolveVertices_UpdatePtrs(Vertex* vert, Edge* edge1, Edge* edge2, Edge* edge3, Edge* edge4, Cell* cell1, Cell* cell2, Cell* cell3, Cell* cell4, Vect* Force){
	PCPEdge* pcpedge1 = dynamic_cast<PCPEdge*>(edge1);
	PCPEdge* pcpedge2 = dynamic_cast<PCPEdge*>(edge2);
	PCPEdge* pcpedge3 = dynamic_cast<PCPEdge*>(edge3);
	PCPEdge* pcpedge4 = dynamic_cast<PCPEdge*>(edge4);
	
	Vect oldvr = (*vert).r;
	Vect newv1r = Vect((*vert).r, (*Force).Angle(), NEWEDGELENGTH / 2.0);		 //find the possition of the new vertices
	Vect newv2r = Vect((*vert).r, (*Force).Angle() + PI, NEWEDGELENGTH / 2.0);

	//get the edge lists
	vector<Edge*> v12edges;
	v12edges.push_back(edge1);
	v12edges.push_back(edge2);
	vector<Edge*> v34edges;
	v34edges.push_back(edge3);
	v34edges.push_back(edge4);

	Vertex* v12ptr = (*this).CreateVertex(oldvr, v12edges);						//create two new vertices at the same position as the old vertex
	Vertex* v34ptr = (*this).CreateVertex(oldvr, v34edges);

	Edge* eptr = 0;
	eptr = (*this).CreateEdge(eptr, 0, 0, v12ptr, v34ptr, cell4, cell2);			//create new edge with vertex and cell ptrs (note is has no flags b/c it can't cross a boundary with length 0)
	PCPEdge* pcp_eptr = dynamic_cast<PCPEdge*>(eptr);

	(*v12ptr).edges.push_back(pcp_eptr);								//add the new edge to the new vertices
	(*v34ptr).edges.push_back(pcp_eptr);
	(*cell2).edges.push_back(pcp_eptr);									//add the new edge to the cells
	(*cell4).edges.push_back(pcp_eptr);
	(*edge1).UpdateVertex(vert, v12ptr);								//uptate the vertices on the edges
	(*edge2).UpdateVertex(vert, v12ptr);
	(*edge3).UpdateVertex(vert, v34ptr);
	(*edge4).UpdateVertex(vert, v34ptr);

	(*this).DestroyVertex(vert);	//destory the fourfold vertex

	//now move the vertices apart along the line of force and set the flags
	(*v12ptr).rnew = newv1r;
	(*v34ptr).rnew = newv2r;
	(*this).HelperResolve_UpdatePosition(v12ptr);
	(*this).HelperResolve_UpdatePosition(v34ptr);

	//Change by Hayden for PCP: Update protein concentration for the new edge as well as the edges joined to the new edge.
	(*this).HelperResolveVertices_RedistributeProteinsNewEdge(pcp_eptr, pcpedge1, pcpedge2, pcpedge3, pcpedge4);

	return eptr;
}

// Reads the data and sets all of the global variables
// Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
void PCPTissue::HelperReadInputData_SetGlobalValues(ifstream* datafile){

	string line, word, key, value;		//temporary variables to store lines and words 
	string type;						//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline;				//ssline is a string stream so we can comvert strings to stringstreams
	map <string, string> globals;


	//Read until you get to the word Global signifying the begining of the global variables
	while (line.compare("Global:") != 0){
		getline(*datafile, line);
		cout << line << endl;
	}

	//Now read the global variables and make a map with the keywords
	while (line.compare("endGlobal") != 0){
		getline(*datafile, line);							//get the next line
		ssline.clear();
		ssline << line;										//put it into a stream
		ssline >> key;										//the first word is the key
		ssline >> value;									// the second is the value of the key
		if (value != ""){										//if the item actually has a value
			globals[key] = value;									// stick the item and it's value in the globals map
		}
	}


	//Now assign the constants their values
	//makes sure the key is actually part of the map before giving the constant a value if the key is not defined it will have the value given at the begining of main
	// the atof is there to convert the string to a double
	if (globals.count("LX") == 1){ LX = atof(globals["LX"].c_str()); }
	else if (globals.count("LX") == 0){ cout << "WARNING: LX is not given in input file, used default value " << LX << endl; }
	else { cout << "ERROR: Global variable LX has multiple definition in input file, used default value LX = " << LX << endl; }

	if (globals.count("LY") == 1){ LY = atof(globals["LY"].c_str()); }
	else if (globals.count("LY") == 0){ cout << "WARNING: LY is not given in input file, used default value " << LY << endl; }
	else { cout << "ERROR: Global variable LX has multiple definition in input file, used default value LY = " << LY << endl; }

	if (globals.count("MINAREA") == 1){ MINAREA = atof(globals["MINAREA"].c_str()); }
	else if (globals.count("MINAREA") == 0){ cout << "WARNING: MINAREA is not given in input file, used default value " << MINAREA << endl; }
	else { cout << "ERROR: Global variable MINAREA has multiple definition in input file, used default value MINAREA = " << MINAREA << endl; }

	if (globals.count("MINLENGTH") == 1){ MINLENGTH = atof(globals["MINLENGTH"].c_str()); }
	else if (globals.count("MINLENGTH") == 0){ cout << "WARNING: LX is not given in input file, used default value " << MINLENGTH << endl; }
	else { cout << "ERROR: Global variable MINLENGTH has multiple definition in input file, used default value MINLENGTH = " << MINLENGTH << endl; }

	if (globals.count("NEWEDGELENGTH") == 1){ NEWEDGELENGTH = atof(globals["NEWEDGELENGTH"].c_str()); }
	else if (globals.count("NEWEDGELENGTH") == 0){ cout << "WARNING: NEWEDGELENGTH is not given in input file, used default value " << NEWEDGELENGTH << endl; }
	else { cout << "ERROR: Global variable NEWEDGELENGTH has multiple definition in input file, used default value NEWEDGELENGTH = " << NEWEDGELENGTH << endl; }

	if (globals.count("GAMMA") == 1){ GAMMA = atof(globals["GAMMA"].c_str()); }
	else if (globals.count("GAMMA") == 0){ cout << "WARNING: GAMMA is not given in input file, used default value " << GAMMA << endl; }
	else { cout << "ERROR: Global variable GAMMA has multiple definition in input file, used default value GAMMA = " << GAMMA << endl; }

	if (globals.count("K") == 1){ K = atof(globals["K"].c_str()); }
	else if (globals.count("K") == 0){ cout << "WARNING: K is not given in input file, used default value " << K << endl; }
	else { cout << "ERROR: Global variable K has multiple definition in input file, used default value K = " << K << endl; }

	if (globals.count("SEED") == 1){ SEED = atof(globals["SEED"].c_str()); }
	else if (globals.count("SEED") == 0){ cout << "WARNING: SEED is not given in input file, used default value " << SEED << endl; }
	else { cout << "ERROR: Global variable SEED has multiple definition in input file, used default value SEED = " << SEED << endl; }

	if (globals.count("INITIALSTEPSIZE") == 1){ INITIALSTEPSIZE = atof(globals["INITIALSTEPSIZE"].c_str()); }
	else if (globals.count("INITIALSTEPSIZE") == 0){ cout << "WARNING: INITIALSTEPSIZE is not given in input file, used default value " << INITIALSTEPSIZE << endl; }
	else { cout << "ERROR: Global variable INITIALSTEPSIZE has multiple definition in input file, used default value INITIALSTEPSIZE = " << INITIALSTEPSIZE << endl; }

	if (globals.count("TIMEINTERVAL") == 1){ TIMEINTERVAL = atof(globals["TIMEINTERVAL"].c_str()); }
	else if (globals.count("TIMEINTERVAL") == 0){ cout << "WARNING: TIMEINTERVAL is not given in input file, used default value " << TIMEINTERVAL << endl; }
	else { cout << "ERROR: Global variable TIMEINTERVAL has multiple definition in input file, used default value TIMEINTERVAL = " << TIMEINTERVAL << endl; }

	if (globals.count("TOTALTIME") == 1){ TOTALTIME = atof(globals["TOTALTIME"].c_str()); }
	else if (globals.count("TOTALTIME") == 0){ cout << "WARNING: TOTALTIME is not given in input file, used default value " << TOTALTIME << endl; }
	else { cout << "ERROR: Global variable TOTALTIME has multiple definition in input file, used default value TOTALTIME = " << TOTALTIME << endl; }

	if (globals.count("CURRENTTIME") == 1){ CURRENTTIME = atof(globals["CURRENTTIME"].c_str()); }
	else if (globals.count("CURRENTTIME") == 0){ cout << "WARNING: CURRENTTIME is not given in input file, used default value " << CURRENTTIME << endl; }
	else { cout << "ERROR: Global variable CURRENTTIME has multiple definition in input file, used default value CURRENTTIME = " << CURRENTTIME << endl; }

	if (globals.count("error_tolerance_absolute_mechanical") == 1){ error_tolerance_absolute_mechanical = atof(globals["error_tolerance_absolute_mechanical"].c_str()); }
	else if (globals.count("error_tolerance_absolute_mechanical") == 0){ cout << "WARNING: error_tolerance_absolute_mechanical is not given in input file, used default value " << error_tolerance_absolute_mechanical << endl; }
	else { cout << "ERROR: Global variable error_tolerance_absolute_mechanical has multiple definition in input file, used default value error_tolerance_absolute_mechanical = " << error_tolerance_absolute_mechanical << endl; }

	if (globals.count("error_tolerance_absolute_protein") == 1){ error_tolerance_absolute_protein = atof(globals["error_tolerance_absolute_protein"].c_str()); }
	else if (globals.count("error_tolerance_absolute_protein") == 0){ cout << "WARNING: error_tolerance_absolute_protein is not given in input file, used default value " << error_tolerance_absolute_protein << endl; }
	else { cout << "ERROR: Global variable error_tolerance_absolute_protein has multiple definition in input file, used default value error_tolerance_absolute_protein = " << error_tolerance_absolute_protein << endl; }

	if (globals.count("ALPHAM") == 1){ alphaM = atof(globals["ALPHAM"].c_str()); }
	else if (globals.count("ALPHAM") == 0){ cout << "WARNING: ALPHAM is not given in input file, used default value " << alphaM << endl; }
	else { cout << "ERROR: Global variable ALPHAM has multiple definition in input file, used default value ALPHAM = " << alphaM << endl; }

	if (globals.count("TISSUE_ALPHA") == 1){ TISSUE_ALPHA = atof(globals["TISSUE_ALPHA"].c_str()); }
	else if (globals.count("TISSUE_ALPHA") == 0){ cout << "WARNING: TISSUE_ALPHA is not given in input file, used default value " << TISSUE_ALPHA << endl; }
	else { cout << "ERROR: Global variable TISSUE_ALPHA has multiple definition in input file, used default value TISSUE_ALPHA = " << TISSUE_ALPHA << endl; }

	if (globals.count("TISSUE_EXT_STRESS_X") == 1){ TISSUE_EXT_STRESS_X = atof(globals["TISSUE_EXT_STRESS_X"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_X") == 0){ cout << "WARNING: TISSUE_EXT_STRESS_X is not given in input file, used default value " << TISSUE_EXT_STRESS_X << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_X has multiple definition in input file, used default value TISSUE_EXT_STRESS_X = " << TISSUE_EXT_STRESS_X << endl; }

	if (globals.count("TISSUE_EXT_STRESS_Y") == 1){ TISSUE_EXT_STRESS_Y = atof(globals["TISSUE_EXT_STRESS_Y"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_Y") == 0){ cout << "WARNING: TISSUE_EXT_STRESS_Y is not given in input file, used default value " << TISSUE_EXT_STRESS_Y << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_Y has multiple definition in input file, used default value TISSUE_EXT_STRESS_Y = " << TISSUE_EXT_STRESS_Y << endl; }

	if (globals.count("EPSAB") == 1){ epsab = atof(globals["EPSAB"].c_str()); }
	else if (globals.count("EPSAB") == 0){ cout << "WARNING: EPSAB is not given in input file, used default value " << epsab << endl; }
	else { cout << "ERROR: Global variable EPSAB has multiple definition in input file, used default value EPSAB = " << epsab << endl; }

	if (globals.count("EPSNAUGHT") == 1){ epsnaught = atof(globals["EPSNAUGHT"].c_str()); }
	else if (globals.count("EPSNAUGHT") == 0){ cout << "WARNING: EPSNAUGHT is not given in input file, used default value " << epsnaught << endl; }
	else { cout << "ERROR: Global variable EPSNAUGHT has multiple definition in input file, used default value EPSNAUGHT = " << epsnaught << endl; }

	if (globals.count("EPSONE") == 1){ epsone = atof(globals["EPSONE"].c_str()); }
	else if (globals.count("EPSONE") == 0){ cout << "WARNING: EPSONE is not given in input file, used default value " << epsone << endl; }
	else { cout << "ERROR: Global variable EPSONE has multiple definition in input file, used default value EPSONE = " << epsone << endl; }

	if (globals.count("EPSTWO") == 1){ epstwo = atof(globals["EPSTWO"].c_str()); }
	else if (globals.count("EPSTWO") == 0){ cout << "WARNING: EPSTWO is not given in input file, used default value " << epstwo << endl; }
	else { cout << "ERROR: Global variable EPSTWO has multiple definition in input file, used default value EPSTWO = " << epstwo << endl; }

	if (globals.count("JONE") == 1){ jone = atof(globals["JONE"].c_str()); }
	else if (globals.count("JONE") == 0){ cout << "WARNING: JONE is not given in input file, used default value " << jone << endl; }
	else { cout << "ERROR: Global variable JONE has multiple definition in input file, used default value JONE = " << jone << endl; }

	if (globals.count("JTWO") == 1){ jtwo = atof(globals["JTWO"].c_str()); }
	else if (globals.count("JTWO") == 0){ cout << "WARNING: JTWO is not given in input file, used default value " << jtwo << endl; }
	else { cout << "ERROR: Global variable JTWO has multiple definition in input file, used default value JTWO = " << jtwo << endl; }

	if (globals.count("TE") == 1){ te = atof(globals["TE"].c_str()); }
	else if (globals.count("TE") == 0){ cout << "WARNING: TE is not given in input file, used default value " << te << endl; }
	else { cout << "ERROR: Global variable TE has multiple definition in input file, used default value TE = " << te << endl; }

	if (globals.count("KT1") == 1){ kT1 = atof(globals["KT1"].c_str()); }
	else if (globals.count("KT1") == 0){ cout << "WARNING: KT1 is not given in input file, used default value " << kT1 << endl; }
	else { cout << "ERROR: Global variable KT1 has multiple definition in input file, used default value KT1 = " << kT1 << endl; }

	if (globals.count("GAMMA0") == 1){ GAMMA0 = atof(globals["GAMMA0"].c_str()); }
	else if (globals.count("GAMMA0") == 0){ cout << "WARNING: GAMMA0 is not given in input file, used default value " << GAMMA0 << endl; }
	else { cout << "ERROR: Global variable GAMMA0 has multiple definition in input file, used default value GAMMA0 = " << GAMMA0 << endl; }

	if (globals.count("GAMMAC") == 1){ GAMMAC = atof(globals["GAMMAC"].c_str()); }
	else if (globals.count("GAMMAC") == 0){ cout << "WARNING: GAMMAC is not given in input file, used default value " << GAMMAC << endl; }
	else { cout << "ERROR: Global variable GAMMAC has multiple definition in input file, used default value GAMMAC = " << GAMMAC << endl; }

	if (globals.count("ALPHAP") == 1){ alphaP = atof(globals["ALPHAP"].c_str()); }
	else if (globals.count("ALPHAP") == 0){ cout << "WARNING: ALPHAP is not given in input file, used default value " << alphaP << endl; }
	else { cout << "ERROR: Global variable ALPHAP has multiple definition in input file, used default value ALPHAP = " << alphaP << endl; }

	if (globals.count("INCREMENT_FOR_RENORMALIZATION") == 1){ increment_for_renormalization = atof(globals["INCREMENT_FOR_RENORMALIZATION"].c_str()); }
	else if (globals.count("INCREMENT_FOR_RENORMALIZATION") == 0){ cout << "WARNING: INCREMENT_FOR_RENORMALIZATION is not given in input file, used default value " << increment_for_renormalization << endl; }
	else { cout << "ERROR: Global variable INCREMENT_FOR_RENORMALIZATION has multiple definition in input file, used default value INCREMENT_FOR_RENORMALIZATION = " << increment_for_renormalization << endl; }

	if (globals.count("MAX_STEP_SIZE") == 1){ max_step_size = atof(globals["MAX_STEP_SIZE"].c_str()); }
	else if (globals.count("MAX_STEP_SIZE") == 0){ cout << "WARNING: MAX_STEP_SIZE is not given in input file, used default value " << max_step_size << endl; }
	else { cout << "ERROR: Global variable MAX_STEP_SIZE has multiple definition in input file, used default value MAX_STEP_SIZE = " << max_step_size << endl; }

	if (globals.count("evolve_mech_positions") == 1){ evolve_mech_positions = atof(globals["evolve_mech_positions"].c_str()); }
	else if (globals.count("evolve_mech_positions") == 0){ cout << "WARNING: evolve_mech_positions is not given in input file, used default value " << evolve_mech_positions << endl; }
	else { cout << "ERROR: Global variable evolve_mech_positions has multiple definition in input file, used default value evolve_mech_positions = " << evolve_mech_positions << endl; }

	if (globals.count("simulation_type") == 1){ simulation_type = atof(globals["simulation_type"].c_str()); }
	else if (globals.count("simulation_type") == 0){ cout << "WARNING: simulation_type is not given in input file, used default value " << simulation_type << endl; }
	else { cout << "ERROR: Global variable simulation_type has multiple definition in input file, used default value simulation_type = " << simulation_type << endl; }

	if (globals.count("max_concentration") == 1){ max_concentration = atof(globals["max_concentration"].c_str()); }
	else if (globals.count("max_concentration") == 0){ cout << "WARNING: max_concentration is not given in input file, used default value " << max_concentration << endl; }
	else { cout << "ERROR: Global variable max_concentration has multiple definition in input file, used default value max_concentration = " << max_concentration << endl; }

	if (globals.count("dirichlet_noise_param") == 1){ dirichlet_noise_param = atof(globals["dirichlet_noise_param"].c_str()); }
	else if (globals.count("dirichlet_noise_param") == 0){ cout << "WARNING: dirichlet_noise_param is not given in input file, used default value " << dirichlet_noise_param << endl; }
	else { cout << "ERROR: Global variable dirichlet_noise_param has multiple definition in input file, used default value dirichlet_noise_param = " << dirichlet_noise_param << endl; }
	
	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//
}

// Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids.
// Starts at current position of *datafile and reads until keyword 'End'.
void PCPTissue::HelperReadInputData_CreateObjects(ifstream* datafile){
	string line, key, value, edg;					//temporary variables to store lines and words 
	string type;							//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline, ssedges;					//ssline is a string stream so we can comvert strings to stringstreams
	vector<string> keywords;					//vector of strings to hold the keywords
	vector<int> listedges;
	int i;
	map <string, string> objectInfo;				//map to hold the keywords and their data
	int newid, v1, v2, enumb, typenumb;				//temp varables to store the object info
	double A0;
	double x, y;
	bool xf, yf;

	cout << "Making a base class tissue" << endl;

	//VERTICES
	//read the file untill you get to the word vertices
	while (line.compare("Vertices") != 0){
		getline(*datafile, line);
	}
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}
	//now recorde the data and make the vertex object untill you get to the word edges signifying the end of the data
	while (true){
		getline(*datafile, line);					//get the next line
		if ((line.compare("Edges") == 0)){ break; }
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')){				//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//make the object and as it to the linked list of vertices also record it's id in the master id list
		newid = atoi(objectInfo["id"].c_str());				//extract the data from the map
		x = atof(objectInfo["x"].c_str());
		y = atof(objectInfo["y"].c_str());
		Vertex* vptr = new Vertex;					// Create a vertes with an id and x,y pos and put it in the master list of all veritces
		(*vptr).id = newid;
		(*vptr).r.x = x;
		(*vptr).r.y = y;
		(vertices).push_back(vptr);
		vertexIDs.push_back(newid);					//and add its id to the list of ids
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*vertices.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}

	//EDGES
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}
	//now recorde the data and make the edge object untill you get to the word cells signifying the end of the data
	while (true){
		getline(*datafile, line);					//get the next line
		if (line.compare("Cells") == 0){ break; }
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')){				//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//make the object and as it to the linked list of vertices also record it's id in the master id list
		newid = atoi(objectInfo["id"].c_str());				//extract the data from the map
		v1 = atoi(objectInfo["vertex1"].c_str());
		v2 = atoi(objectInfo["vertex2"].c_str());
		xf = atoi(objectInfo["xflag"].c_str());
		yf = atoi(objectInfo["yflag"].c_str());
		PCPEdge* eptr = new PCPEdge;
		(*eptr).id = newid;
		(*eptr).IOVertexID[0] = v1;
		(*eptr).IOVertexID[1] = v2;
		(*eptr).xflag = xf;
		(*eptr).yflag = yf;
		(*eptr).cell1 = NULL;
		(*eptr).cell2 = NULL;
		//Initialize PCP concentrations to be equal to zero.
		(*eptr).planarpolarity.resize(4);
		fill((*eptr).planarpolarity.begin(), (*eptr).planarpolarity.end(), 0.0);
		(edges).push_back(eptr);					// Create a vertes with an id and x,y pos and put it in the master list of all veritces
		edgeIDs.push_back(newid);					//and add its id to the list of ids
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*edges.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}

	//CELLS
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}
	//now recorde the data and make the cell object untill you get to the word end signifying the end of the data
	while (true){
		getline(*datafile, line);					//get the next line
		if (line.compare("End") == 0){ break; }				//break if keyworkd 'End' is found
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')){				//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//extract the data from the map
		newid = atoi(objectInfo["id"].c_str());
		enumb = atoi(objectInfo["EdgeNumb"].c_str());
		if (objectInfo.count("Type")){
			typenumb = atoi(objectInfo["Type"].c_str());
		}
		else{ typenumb = standard_cell_type; }
		if (objectInfo.count("A0")){
			A0 = atof(objectInfo["A0"].c_str());
		}
		ssedges.clear();
		ssedges << objectInfo["EdgeList"];
		//now worry about the edge list
		listedges.clear();
		while (getline(ssedges, edg, ',')){
			listedges.push_back(atoi(edg.c_str()));
		}
		//make the cell object
		PCPCell* cptr = new PCPCell;
		(*cptr).id = newid;
		(*cptr).type = typenumb;
		(*cptr).A0 = A0;
		(*cptr).IOEdgeID = listedges;
		(*cptr).time_of_pcp = 0.0;
		(*cptr).cell_generation = -1;
		cells.push_back(cptr);
		cellIDs.push_back(newid);
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}
}

Vertex* PCPTissue::HelperT1_RemoveEdge(Edge* deledge, ofstream* errorfile, int* nT1){
	int newfold;							//is the new vertex 4 or 5 fold
	int i, j, k, l, m;
	Vertex *v1ptr, *v2ptr, *vnewptr;		//pointers to the two merging vertices
	Cell *c1ptr, *c2ptr;					//pointers to the two cells that will loose an edge (the other cells remain uneffected b/c their edges do not change and they do not know their vertices)
	vector<PCPEdge*> newedgelist; 
	vector<Edge*> casted_newedgelist;
	PCPEdge* casted_deledge = dynamic_cast<PCPEdge*>(deledge);

	v1ptr = (*casted_deledge).vertex1;
	v2ptr = (*casted_deledge).vertex2;
	c1ptr = (*casted_deledge).cell1;
	c2ptr = (*casted_deledge).cell2;

	//If this will create a 5 fold vertex give an error
	newfold = (*v1ptr).GetEdgeNumber() + (*v2ptr).GetEdgeNumber() - 2;
	if (newfold>4){
		cout << "ERROR: DID NOT COLLAPSE AN EDGE TO GET A 5FOLD VERTEX " << endl;
		return NULL;
	}

	//If this will create a cell with only two edges give an error (ie if the cell on either side of the deleted edge is a triangle)
	if ((*c1ptr).GetVertexNumber()<4 || (*c2ptr).GetVertexNumber()<4){
		cout << "ERROR: DID NOT COLLAPSE AN EDGE THAT IS PART OF A TRIANGLE" << endl;
		return NULL;
	}

	//If there are no errors continue

	//Record the number of t1s
	*nT1 = (*deledge).numbT1;

	//fill the new edge list
	for (i = 0; i<(*v1ptr).edges.size(); i++){					//get the edges from vertex1
		if ((*casted_deledge) != (*(*v1ptr).edges[i])){
			PCPEdge* pcpptr = dynamic_cast<PCPEdge*>((*v1ptr).edges[i]);
			newedgelist.push_back(pcpptr);
		}
	}
	for (j = 0; j<(*v2ptr).edges.size(); j++){					//get the edges from vertex 2
		if ((*casted_deledge) != (*(*v2ptr).edges[j])){
			PCPEdge* pcpptr = dynamic_cast<PCPEdge*>((*v2ptr).edges[j]);
			newedgelist.push_back(pcpptr);
		}
	}

	//Changes by Hayden: Make sure that the proteins on the edge to be deleted are redistributed to neighboring edges. 
	(*this).HelperT1_RedistributeProteinsDeletedEdge(casted_deledge, newedgelist);

	//FLAGS SET**********************************************************
	//find the new position and move the vertices there
	(*casted_deledge).MoveVertsToCenter();
	(*this).HelperResolve_UpdatePosition(v1ptr);
	(*this).HelperResolve_UpdatePosition(v2ptr);

	casted_newedgelist.resize(newedgelist.size());

	for (int i = 0; i != newedgelist.size(); i++)
		casted_newedgelist.at(i) = dynamic_cast<Edge*>(newedgelist.at(i));

	//make the new vertex
	vnewptr = (*this).CreateVertex((*casted_deledge).GetCenter(), casted_newedgelist);

	//go to the four edges and change their ptrs from the old vertices to the new fourfold vertex
	for (k = 0; k< newedgelist.size(); k++){
		if ((*(*newedgelist[k]).vertex1).id == (*v1ptr).id){
			(*newedgelist[k]).vertex1 = vnewptr;
		}
		else if ((*(*newedgelist[k]).vertex1).id == (*v2ptr).id){
			(*newedgelist[k]).vertex1 = vnewptr;
		}
		else if ((*(*newedgelist[k]).vertex2).id == (*v1ptr).id){
			(*newedgelist[k]).vertex2 = vnewptr;
		}
		else if ((*(*newedgelist[k]).vertex2).id == (*v2ptr).id){
			(*newedgelist[k]).vertex2 = vnewptr;
		}
		else { cout << "ERROR IN DELETE EDGE" << endl; }
	}

	//Remove the old vertices from the master list
	(*this).DestroyVertex(v1ptr);
	(*this).DestroyVertex(v2ptr);

	//Remove the edge from the 2 cell lists of edges
	for (l = 0; l<(*c1ptr).edges.size(); l++){
		if ((*(*c1ptr).edges[l]).id == (*casted_deledge).id) {
			(*c1ptr).edges.erase((*c1ptr).edges.begin() + l);
			break;
		}
	}
	for (m = 0; m<(*c2ptr).edges.size(); m++){
		if ((*(*c2ptr).edges[m]).id == (*casted_deledge).id) {
			(*c2ptr).edges.erase((*c2ptr).edges.begin() + m);
			break;
		}
	}

	//Remove the edge from the master list
	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){					//loop through all of the edges in the tissue
		if ((*(*eit)).id == (*casted_deledge).id){								//if the curent edge is to be removed
			DestroyEdge(*eit);
			break;
		}
	}

	//return a pointer to the vertex
	return vnewptr;
}

//! Initializes all of the pointers so that each object has pointers to it's neighbors.
void PCPTissue::HelperReadInputData_CreatePointers(){
	list<Vertex*>::iterator vit;
	cout << "Making cell <-> edge pointers!" << endl;
	
	//first make sure pointers are null (if using subclasses this might not have happened corectly) 
	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){
		(*(*eit)).cell1 = NULL;
		(*(*eit)).cell2 = NULL;
	}


	//Give the cells pointers to their edges and the edges pointers back to the cells
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){				//for all of the cells
		vector<int> eIDs;										//make a vector to store the edgeIDs
		eIDs = (*(*cit)).IOEdgeID;									//get the IDs of the initial edges

		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//look through all edges
			if (find(eIDs.begin(), eIDs.end(), (*(*eit)).id) != eIDs.end()){			//if the curent edge id is in the list of ids

				((*(*cit)).edges).push_back(&(*(*eit)));					//add edge pointer to cell

				if ((*(*eit)).cell1 == NULL){							//add cell ptr to the edge
					(*(*eit)).cell1 = &(*(*cit));
				}
				else{
					(*(*eit)).cell2 = &(*(*cit));
				}
			}
		}
	}

	cout << "Making vertex <-> edge pointers!" << endl;
	//Give the edges pointers to their vertices and the vertices pointers back to the edges
	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){				//for all of the edges
		for (vit = vertices.begin(); vit != vertices.end(); vit++){					//and all of the vertices

			if ((*(*vit)).id == (*(*eit)).IOVertexID[0]){						//if the vertex id matches set up pointers from the edge to the vertex and vice versa
				(*(*eit)).vertex1 = &(*(*vit));
				((*(*vit)).edges).push_back(&(*(*eit)));
			}

			if ((*(*vit)).id == (*(*eit)).IOVertexID[1]){
				(*(*eit)).vertex2 = &(*(*vit));
				((*(*vit)).edges).push_back(&(*(*eit)));
			}
		}
	}
}

//! This is a temporary function used for allowing the input of concentrations by the user.
//! Initializes all of the pointers so that each object has pointers to it's neighbors.
void PCPTissue::HelperReadInputData_ProteinConcentrations(ifstream* datafile){
	string line, key, value, edg;										//temporary variables to store lines and words 
	string type;												//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline, ssedges;										//ssline is a string stream so we can comvert strings to stringstreams
	vector<string> keywords;										//vector of strings to hold the keywords
	vector<int> listedges;
	int i;
	map <string, string> objectInfo;									//map to hold the keywords and their data
	vector<vector<int>> vec_init(2,vector<int>(edges.size(), 0));
	vector<vector<int>> cell12_id_list = vec_init;								//Store the ids
	int cell1id;
	int cell2id;
	double init_val = 0.0;
	vector<vector<double>> vec(2, vector<double>(edges.size(), init_val));
	vector<vector<double>> concA_list = vec;								//vector for storing order list of protein concentrations of type A for cell1, cell2 (as defined in input file).
	vector<vector<double>> concB_list = vec;								//vector for storing order list of protein concentrations of type A for cell2, cell1 (as defined in input file).
	unordered_map<long, int> edge_id2index = (*this).CreateEdgeId2IndexMap();
	int newid;

	//EDGES
	//read the file until you get to the word edges
	while (line.compare("Edges") != 0){
		getline(*datafile, line);
	}
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}

	//now recorde the data and make the edge object untill you get to the word cells signifying the end of the data
	while (true){
		getline(*datafile, line);									//get the next line
		if (line.compare("Cells") == 0){ break; }
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')){								//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//make the object and as it to the linked list of vertices also record it's id in the master id list
		newid = atoi(objectInfo["id"].c_str());								//extract the data from the map
		int e_ind = edge_id2index.at(newid);
		cell1id = atoi(objectInfo["_cell1"].c_str());
		cell2id = atoi(objectInfo["_cell2"].c_str());

		cell12_id_list[0].at(e_ind) = cell1id;
		cell12_id_list[1].at(e_ind) = cell2id;
		concA_list[0].at(e_ind) = atof(objectInfo["_cell1AProtein"].c_str());
		concA_list[1].at(e_ind) = atof(objectInfo["_cell2AProtein"].c_str());
		concB_list[0].at(e_ind) = atof(objectInfo["_cell1BProtein"].c_str());
		concB_list[1].at(e_ind) = atof(objectInfo["_cell2BProtein"].c_str());
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*edges.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}

	//Now iterate through the list of edge pointers and assign protein concentrations (making sure that the concentrations are assigned to the correct cells).
	int ind = 0;
	for (list<Edge*>::iterator it = edges.begin(); it != edges.end(); it++){
		PCPEdge* pcpeit = dynamic_cast<PCPEdge*>(*it);
		if ((*(*pcpeit).cell1).id == cell12_id_list[0].at(ind)){
			(*pcpeit).planarpolarity[0] = concA_list[0].at(ind);
			(*pcpeit).planarpolarity[1] = concB_list[0].at(ind);
			(*pcpeit).planarpolarity[2] = concA_list[1].at(ind);
			(*pcpeit).planarpolarity[3] = concB_list[1].at(ind);
		}
		else{
			(*pcpeit).planarpolarity[0] = concA_list[1].at(ind);
			(*pcpeit).planarpolarity[1] = concB_list[1].at(ind);
			(*pcpeit).planarpolarity[2] = concA_list[0].at(ind);
			(*pcpeit).planarpolarity[3] = concB_list[0].at(ind);
		}
		ind++;
	}
}

// This builds a list of cells that helps to determine fate specification in the next row.
vector<Cell*> HelperCheckerboard_BuildList_ConnectedCells(vector<PCPCell*> next_differentiated_cells){
	vector<Cell*> temp_connected_cells;																//Store vector of pointers to cells connected to the given cell.
	vector<Cell*> connected_cells;																	//Returns full vector of pointers to cells used to determine type determination.
	int number_of_cells_to_differentiate = next_differentiated_cells.size();						//Number of cells to differentiate.
	PCPCell* pcpcellptr;

	//Push back the line of cells to differentiate.
	for (int differentiate_cell_index = 0; differentiate_cell_index != number_of_cells_to_differentiate; differentiate_cell_index++){
		connected_cells.push_back(next_differentiated_cells.at(differentiate_cell_index));
	}

	for (int differentiate_cell_index = 0; differentiate_cell_index != number_of_cells_to_differentiate; differentiate_cell_index++){
		//Store pointers to cells connected to the given cell to differentiate.
		temp_connected_cells = (*next_differentiated_cells.at(differentiate_cell_index)).GetNeighboringCells();
		for (int temp_cell_index = 0; temp_cell_index != temp_connected_cells.size(); temp_cell_index++){
			pcpcellptr = dynamic_cast<PCPCell*>(temp_connected_cells.at(temp_cell_index));
			//Store cells that have PCP.
			if ((*pcpcellptr).is_Cell_with_PCP()){
				connected_cells.push_back(temp_connected_cells.at(temp_cell_index));
			}
		}
	}
	//sort the pointers.
	std::sort(connected_cells.begin(), connected_cells.end());

	//Only keep unique cells
	auto last = std::unique(connected_cells.begin(), connected_cells.end());
	connected_cells.erase(last, connected_cells.end());

	return connected_cells;
}

//! This is a helper function for generating a checkerboard pattern in the next row: for periodic boundary conditions with a vertical stripe rather than in a circular ring.
void PCPTissue::HelperCheckerboard_GeneratePattern(vector<PCPCell*> next_differentiated_cells, ofstream* datafile){
	//Order the differentiated cells based on angle: this helps later for flipping domains of spin.
	vector<PCPCell*> ordered_next_differentiated_cells = order_cells_position(next_differentiated_cells);
	int number_of_next_differentiated_cells = ordered_next_differentiated_cells.size();
	//Generate random distribution of spins along the ring.
	vector<double> next_diff_cells_types = HelperCheckerboard_GenerateRandomSpins(number_of_next_differentiated_cells);

	//Build list of all differentiated cells connected to the next_differentiated_cells (This list includes next_differentiated_cells).
	vector<Cell*> connected_cells = HelperCheckerboard_BuildList_ConnectedCells(ordered_next_differentiated_cells);
	int number_of_connected_cells = connected_cells.size();

	//Build adjacency matrix that includes contact lengths between cells (zeros on the diagonal).
	gsl_matrix* contact_matrix = gsl_matrix_calloc(number_of_connected_cells, number_of_connected_cells);
	HelperCheckerboard_BuildAdjacencyMatrix(ordered_next_differentiated_cells, contact_matrix);

	//Vector for storing the different types of the connected cells.
	gsl_vector* cell_type_vectors = gsl_vector_calloc(number_of_connected_cells);
	//J matrix times cell_type_vectors.
	gsl_vector* coupling_vector = gsl_vector_calloc(number_of_connected_cells);

	//Build id2index map for the connected cells.
	unordered_map<long, int> connectedcells_id2index;
	connectedcells_id2index.reserve(number_of_connected_cells);								//reserve to contain vertices.size() elements

	//Store initial random configuration.
	vector<double> min_energy_configuration;
	//Store calculated energy with the given configuration.
	//calculated_energy -> energy calculated for a given configuration.
	//current_energy -> store energy of current configuration.
	//min_energy -> store energy of minimum energy (found during the entire search thus far).
	//delta_energy -> store change in energy between current configuration and proposed change.
	double calculated_energy, current_energy, min_energy, delta_energy, delta_min_energy;

	//Store a history of past energy in order to see if the system has equilibrized during single-spin-flips.
	vector<double> past_energy_vector;
	//randomly chosen index for spin flip.
	int random_spin_index;
	//define uniform integer distribution
	uniform_int_distribution<int> distribution(0, RAND_MAX);

	//calculate the characteristic length (coupling energies are based on contact lengths : only use the area filled by actual cells)
	double area_of_packing = LX * LY;
	double characteristic_length = sqrt(area_of_packing / (*this).cells.size() / 3 / sqrt(3) * 2);
	double characteristic_energy_scale = characteristic_length;
	//variable for storing standard deviation of the past_energy_vector.
	double history_st_dev = characteristic_energy_scale;
	//store inverse temperature (~inverse of characteristic energy)
	double beta = 2. / characteristic_length;
	//length of history to be stored: sample time scale of approximately ~10 fluctuations from a metastable state.
	int length_of_history = (int)((1 + exp(-4. * beta * characteristic_energy_scale)) / exp(-4. * beta * characteristic_energy_scale));

	//Store map for id 2 index of the coupling matrix.
	int ind = 0;
	for (int connected_cell_index = 0; connected_cell_index != number_of_connected_cells; connected_cell_index++){
		connectedcells_id2index[(*connected_cells.at(connected_cell_index)).id] = ind;
		ind++;
	}

	//For cells that already have PCP that are connected to the next differentiated cells, store the corresponding spin (type 1 <-> -1.0, type 2 <-> 1.0)
	for (int cell_index = 0; cell_index != number_of_connected_cells; cell_index++){
		if ((*connected_cells.at(cell_index)).type == 1){
			gsl_vector_set(cell_type_vectors, cell_index, -1.0);
		}
		else if ((*connected_cells.at(cell_index)).type == 2){
			gsl_vector_set(cell_type_vectors, cell_index, 1.0);
		}
	}

	//Store the "spins" from randomly generated initial distribution.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
		gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(cell_index));
	}
	//Multiply the matrix of couplings by the spin vector (store result in coupling vector).
	gsl_blas_dgemv(CblasNoTrans, 1.0, contact_matrix, cell_type_vectors, 0.0, coupling_vector);
	//Do the dot product of coupling vector and spin vector -> store result in calculated energy.
	gsl_blas_ddot(cell_type_vectors, coupling_vector, &calculated_energy);
	//Store the initial energy as initial min energy, and initialize min_energy_configuration as initial "spin distribution"
	min_energy = calculated_energy;
	current_energy = min_energy;
	past_energy_vector.push_back(current_energy);
	min_energy_configuration = next_diff_cells_types;

	//Output initial type ("spin") distribution.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		if (min_energy_configuration.at(cell_index) == -1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 1;
		}
		else if (min_energy_configuration.at(cell_index) == 1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 2;
		}
	}
	(*this).WriteOutputDataAppend(datafile);

	//Keep track the number of attempts of spin flips (stop when you both satisfy the equilibration of spins and have at least gone length_of_history steps).
	int number_of_attempts = 0;
	do{
		//Select a random index for flipping the corresponding spin.
		random_spin_index = distribution(generator) % number_of_next_differentiated_cells;
		//Find absolute cell index -> in the convention of the coupling matrix.
		int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(random_spin_index)).id];
		//Flip that spin in the cell_type_vectors.
		gsl_vector_set(cell_type_vectors, absolute_cell_index, -1. * next_diff_cells_types.at(random_spin_index));
		//Multiply the matrix of couplings by the spin vector (store result in coupling vector).
		gsl_blas_dgemv(CblasNoTrans, 1.0, contact_matrix, cell_type_vectors, 0.0, coupling_vector);
		//Do the dot product of coupling vector and spin vector -> store result in calculated energy.
		gsl_blas_ddot(cell_type_vectors, coupling_vector, &calculated_energy);
		//Calculate the change in energy brought about by this spin flip (relative to the current configuration and the minimum energy configuration thus far).
		delta_energy = calculated_energy - current_energy;
		delta_min_energy = calculated_energy - min_energy;
		//Calculate the probability of the spin flip (based on delta_energy with Boltzmann weighting).
		double probability_flip = exp(-beta * delta_energy) / (1 + exp(-beta * delta_energy));
		//If this is currently the min energy state, store the configuration and the energy.
		if (delta_energy < 0.0){
			//Flip the corresponding spin.
			next_diff_cells_types.at(random_spin_index) = -1. * next_diff_cells_types.at(random_spin_index);
			//Store the calculated energy as the new minimum, and store the configuration.
			current_energy = calculated_energy;
		}
		//If the change in energy is positive, use a boltzmann weighting to judge if the change should proceed.
		else if (delta_energy > 0.0){
			bernoulli_distribution bern_distribution(probability_flip); 														//distribution for weighting chance of success
			//Use boltzmann weighting.
			if (bern_distribution(generator)){
				//Flip the corresponding spin.
				next_diff_cells_types.at(random_spin_index) = -1. * next_diff_cells_types.at(random_spin_index);
				//Store this energy as the current energy (but not the minimum).
				current_energy = calculated_energy;
			}
			else{
				//Change the spin in the cell type vectors back to the original state.
				gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(random_spin_index));
			}
		}
		//Increment the index for the number of attempts.
		number_of_attempts++;

		//If this configuration's energy is less than the currently stored minimum, update the minimum.
		if (delta_min_energy < 0.0){
			min_energy = calculated_energy;
			min_energy_configuration = next_diff_cells_types;
		}

		//If the number of attempts is less than the length of the history, push back the energy.
		if (number_of_attempts < length_of_history){
			past_energy_vector.push_back(current_energy);
		}
		//Else, update the energy vector on a rotating basis (dropping off oldest stored energy, adding new energy).
		else{
			vector<double> updated_energy_vector;
			for (int i = 1; i != length_of_history; i++) updated_energy_vector.push_back(past_energy_vector.at(i));
			updated_energy_vector.push_back(current_energy);
			past_energy_vector = updated_energy_vector;
			//Store the standard deviation of the history.
			history_st_dev = sqrt(variance_history(past_energy_vector));
		}
	} while (history_st_dev > characteristic_energy_scale || number_of_attempts < length_of_history + 1);

	//Output the configuration after this individual spin flip step.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		if (min_energy_configuration.at(cell_index) == -1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 1;
		}
		else if (min_energy_configuration.at(cell_index) == 1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 2;
		}
	}
	(*this).WriteOutputDataAppend(datafile);

	//Find domains of spin, attempt to flip domains: start with min energy configuration.
	next_diff_cells_types = min_energy_configuration;
	current_energy = min_energy;
	//Store this configuration in the cell_type_vectors.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
		gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(cell_index));
	}

	//Find mismatches in cell types in the chain of cells.
	vector<int> mismatched_cells = find_domains_of_spin(next_diff_cells_types);
	//Only attempt the number of domain flips consistent with the number of mismatches.
	int number_of_domain_flip_attempts = mismatched_cells.size();

	//Keep track the number of attempts of spin domain flips (stop when you reach number_of_domain_flip_attempts).
	int spin_flip_attempt = 0;
	while (spin_flip_attempt != number_of_domain_flip_attempts && mismatched_cells.size() != 0){
		//Find mismatches in cell types in the chain of cells.
		mismatched_cells = find_domains_of_spin(next_diff_cells_types);
		//Select a domain randomly by randomly selecting the corresponding mismatch.
		int random_mismatch_index = distribution(generator) % mismatched_cells.size();

		//Start at the domain mismatch.
		int cell_index = mismatched_cells.at(random_mismatch_index);
		int mod_cell_index;
		//Loop through flipping the full domain of spins until the next mismatch.
		do{
			//As you increment the cell index, use mod to wrap through the list of cells (the cells are in a ring).
			mod_cell_index = cell_index % number_of_next_differentiated_cells;
			//Find absolute cell index -> in the convention of the coupling matrix.
			int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(mod_cell_index)).id];
			//Flip the spin in the cell_type_vectors.
			gsl_vector_set(cell_type_vectors, absolute_cell_index, -1 * next_diff_cells_types.at(mod_cell_index));
			//Increment cell index.
			cell_index++;
			mod_cell_index = cell_index % number_of_next_differentiated_cells;
		} while (mod_cell_index != mismatched_cells.at((random_mismatch_index + 1) % mismatched_cells.size()));
		gsl_blas_dgemv(CblasNoTrans, 1.0, contact_matrix, cell_type_vectors, 0.0, coupling_vector);
		//Do the dot product of coupling vector and spin vector -> store result in calculated energy.
		gsl_blas_ddot(cell_type_vectors, coupling_vector, &calculated_energy);
		//Calculate the change in energy brought about by this spin flip (relative to the current configuration and the minimum energy configuration thus far).
		delta_energy = calculated_energy - current_energy;
		delta_min_energy = calculated_energy - min_energy;
		//Calculate the probability of the spin domain flip (based on delta_energy with Boltzmann weighting).
		double probability_flip = exp(-beta * delta_energy) / (1 + exp(-beta * delta_energy));
		//Increment spin_flip_attempt index.
		spin_flip_attempt++;

		//If this is currently the min energy state, store the configuration and the energy.
		if (delta_energy < 0.0){
			//Update the vector of next_diff_cells.
			for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
				int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
				next_diff_cells_types.at(cell_index) = gsl_vector_get(cell_type_vectors, absolute_cell_index);
			}
			//Store this energy as the current energy (but not the minimum).
			current_energy = calculated_energy;
		}
		//If the change in energy is positive, use a boltzmann weighting to judge if the change should proceed.
		else if (delta_energy > 0.0){
			bernoulli_distribution bern_distribution(probability_flip); 																//distribution for weighting chance of success
			//Use boltzmann weighting.
			if (bern_distribution(generator)){
				//Flip the corresponding domain of spins.
				for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
					int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
					next_diff_cells_types.at(cell_index) = gsl_vector_get(cell_type_vectors, absolute_cell_index);
				}
				//Update current energy.
				current_energy = calculated_energy;
			}
			else{
				//If not, reset the cell_type_vectors.
				cell_index = mismatched_cells.at(random_mismatch_index);
				do{
					mod_cell_index = cell_index % number_of_next_differentiated_cells;
					int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(mod_cell_index)).id];
					gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(mod_cell_index));
					cell_index++;
					mod_cell_index = cell_index % number_of_next_differentiated_cells;
				} while (mod_cell_index != mismatched_cells.at((random_mismatch_index + 1) % mismatched_cells.size()));
			}
		}

		//If this configuration's energy is less than the currently stored minimum, update the minimum.
		if (delta_min_energy < 0.0){
			min_energy = calculated_energy;
			min_energy_configuration = next_diff_cells_types;
		}
		mismatched_cells = find_domains_of_spin(next_diff_cells_types);
	}

	//Now differentiate the cells (based on minimum energy configuration).
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		if (min_energy_configuration.at(cell_index) == -1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 1;
		}
		else if (min_energy_configuration.at(cell_index) == 1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 2;
		}
	}
	(*this).WriteOutputDataAppend(datafile);

	//Free the gsl allocations.
	gsl_matrix_free(contact_matrix);
	gsl_vector_free(coupling_vector);
	gsl_vector_free(cell_type_vectors);
}
// Build the adjacency matrix (stores the values of contact lengths between cells in contact_matrix).
void HelperCheckerboard_BuildAdjacencyMatrix(vector<PCPCell*> next_differentiated_cells, gsl_matrix* contact_matrix){
	//Used to temporarily store list of cells connected to a given next_differentiated_cell.
	vector<Cell*> temp_connected_cells;
	//Used to store contact lengths for temp_connected_cells (follows same ordering conventions)
	vector<double> contact_length_connected_cells;
	//Store list of all connected cells (that will determine the types of next_differentiated_cells)
	vector<Cell*> connected_cells = HelperCheckerboard_BuildList_ConnectedCells(next_differentiated_cells);
	//Number of cells to differentiate in the next round.
	int number_of_cells_to_differentiate = next_differentiated_cells.size();
	//Map from id 2 index for contact matrix.
	unordered_map<long, int> connectedcells_id2index;
	PCPCell* pcpcellptr;

	//Store the id 2 index (for full list of connected cells)
	int number_of_connected_cells = connected_cells.size();
	connectedcells_id2index.reserve(number_of_connected_cells);							//reserve to contain vertices.size() elements

	//Define indices for contact lengths.
	int ind = 0;
	for (int connected_cell_index = 0; connected_cell_index != number_of_connected_cells; connected_cell_index++){
		connectedcells_id2index[(*connected_cells.at(connected_cell_index)).id] = ind;
		ind++;
	}

	//Loop through all of the connected cells.
	for (int connected_cell_index = 0; connected_cell_index != number_of_cells_to_differentiate; connected_cell_index++){
		//Store the index in the contact matrix that corresponds to the given cell in the next_differentiated_cells list.
		int absolute_connected_cell_index = connectedcells_id2index[(*next_differentiated_cells.at(connected_cell_index)).id];
		//Store the list of cells connected to the the given next differentiated cell.
		temp_connected_cells = (*next_differentiated_cells.at(connected_cell_index)).GetNeighboringCells();
		//Store contact lengths for temp_connected cells.
		contact_length_connected_cells = (*next_differentiated_cells.at(connected_cell_index)).GetNeighboringCells_ContactLength();
		//Loop through temp connected cells and modify contact matrix.
		for (int temp_cell_index = 0; temp_cell_index != temp_connected_cells.size(); temp_cell_index++){
			pcpcellptr = dynamic_cast<PCPCell*>(temp_connected_cells.at(temp_cell_index));
			//Make sure that you modify the contact matrix in a way that is symmetric about the diagonal.
			//Only if the given cell is differentiated or next to differentiate, modify the contact matrix.
			if ((*pcpcellptr).is_Cell_with_PCP()){
				int temp_connected_cell_index = connectedcells_id2index[(*temp_connected_cells.at(temp_cell_index)).id];
				gsl_matrix_set(contact_matrix, temp_connected_cell_index, absolute_connected_cell_index, contact_length_connected_cells.at(temp_cell_index));
				gsl_matrix_set(contact_matrix, absolute_connected_cell_index, temp_connected_cell_index, contact_length_connected_cells.at(temp_cell_index));
			}
		}
	}
}
// Generate matrix of random spins for finding a low energy state (anti-ferromagnetic). This is the initial state for the metropolis algorithm.
// number_of_spins ->number of cells in the next line for type specification.
vector<double> HelperCheckerboard_GenerateRandomSpins(int number_of_spins){
	//Store a matrix of all configuration with a given number of spin ups and total spins.
	vector<double> randomly_generated_spins;
	randomly_generated_spins.resize(number_of_spins);

	double p = 0.5;										//probability of mitosis of any individual cell
	bernoulli_distribution distribution(p); 			//distribution for weighting chance of success

	for (int spin_index = 0; spin_index != number_of_spins; spin_index++){
		if (distribution(generator)){
			randomly_generated_spins.at(spin_index) = -1.;
		}
		else{
			randomly_generated_spins.at(spin_index) = 1.;
		}
	}

	return randomly_generated_spins;
}

//! Order cell pointer by vertical position (try to identify chains of cells).
vector<PCPCell*> order_cells_position(const vector<PCPCell*> &unordered_cells){
	//Store the number of cells in the list.
	int number_of_cells = unordered_cells.size();
	//Build a vector for storing angles of cells relative to central mega cell.
	vector<double> position_vector;
	position_vector.resize(number_of_cells);
	//Build a vector for storing ordered cells (pointers).
	vector<PCPCell*> ordered_cells;
	ordered_cells.resize(number_of_cells);
	//Store a list of angle -> cell.
	unordered_map<double, PCPCell*> cell_angle2pointer;
	cell_angle2pointer.reserve(number_of_cells);							//reserve to contain vertices.size() elements

	//Store the angle for each cell in the unordered_cells list.
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		Vect center_position = (*unordered_cells.at(cell_index)).GetCenter();
		position_vector.at(cell_index) = center_position.y;
	}

	//Build map of angle->cell.
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		cell_angle2pointer[position_vector.at(cell_index)] = unordered_cells.at(cell_index);
	}

	//Sort the pointers by angle.
	std::sort(position_vector.begin(), position_vector.end());

	//Store the ordered cells (by angle).
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		ordered_cells.at(cell_index) = cell_angle2pointer[position_vector.at(cell_index)];
	}

	return ordered_cells;
}

//! This function searches through edge_candidates (edges along the differentiation boundary).
//! Output an edge from this list that faces to the right. This is used within AddDifferentiatedCellColumn (to establish front propogation to the right).
Edge* edge_candidate_to_right(const vector<Edge*> &edge_candidates){
	//Store all edges in edge_candidates that form the right boundary for PCP propogation (make sure that front propogation is consistently in one direction).
	vector<Edge*> edges_to_right;

	//Loop through edge_candidates.
	for (int edge_index = 0; edge_index != edge_candidates.size(); edge_index++){
		//Define pointer to cell 1 and cell2.
		Cell* cell1ptr = (*edge_candidates.at(edge_index)).cell1;
		Cell* cell2ptr = (*edge_candidates.at(edge_index)).cell2;
		Cell* type0_cellptr;

		//Figure out which cell is of type standard_cell_type.
		if ((*cell1ptr).type == standard_cell_type){
			type0_cellptr = cell1ptr;
		}
		else if ((*cell2ptr).type == standard_cell_type){
			type0_cellptr = cell2ptr;
		}

		//Returns a vector of pointers to all vertices in the cell in clockwise (CW) order: for the cell of type standard_cell_type.
		vector<Vertex*> cw_verts = (*type0_cellptr).GetClockWiseVertices();
		//Store corresponding ordered vert positions (respecting PBC).
		vector<vector<double>> ordered_vert_positions = (*type0_cellptr).HelperOrderVerticesPos(cw_verts);
		//Define the geocenter of the ordered_vert_positions.
		Vect position = GetGeoCenter(ordered_vert_positions[0], ordered_vert_positions[1]);
		double angle = 0.0;
		//Store a vector of indices for edges that a ray pointing to the right intersects.
		vector<int> intersection_indices = which_Edge_Intersection(position, ordered_vert_positions, angle);
		//For a convex polygon, there should only one intersection.
		assert((intersection_indices.size() == 1) && "ERROR: COULD NOT DETERMINE DIRECTION FOR FRONT PROPAGATION");
		int intersection_index = intersection_indices.at(0);
		//Store the number of vertices in this cell (type0_cell).
		int number_of_vertices = cw_verts.size();
		//Find next index in the list of vertices (wrapping to beginning of vector by mod).
		int next_index = (intersection_index + 1) % number_of_vertices;
		//Find the edge between the two vertices (the intersected edge).
		Edge* edge_between = (*cw_verts.at(intersection_index)).GetEdgeBetween(cw_verts.at(next_index));
		//Find the cell on the other side of this edge.
		Cell* other_cell = (*edge_between).GetOtherCell(type0_cellptr);
		//If the other cell is of type 0, edge_candidates.at(edge_index) can serve for propogation of PCP to the right.
		if ((*other_cell).type == standard_cell_type){
			edges_to_right.push_back(edge_candidates.at(edge_index));
		}
	}

	//Assert that you have found at least one edge candidate for propogating PCP to the right.
	assert((edges_to_right.size() != 0) && "ERROR: COULD NOT DETERMINE DIRECTION FOR FRONT PROPAGATION");

	//Return only one of the edges.
	Edge* chosen_edge_candidate = edges_to_right.at(0);

	return chosen_edge_candidate;
}

//! Check if the cell is a faulty differentiation (use right-propogation convention).
bool PCPTissue::check_for_faulty_cell(Cell* cell_ptr){
	bool is_faulty_cell = true;

	vector<int> bordering_cell_types;
	vector<Cell*> bordering_cells = (*cell_ptr).GetNeighboringCells();
	int number_of_bordering_cells = bordering_cells.size();
	bordering_cell_types.resize(number_of_bordering_cells);

	for (int cell_index = 0; cell_index != number_of_bordering_cells; cell_index++)
		bordering_cell_types.at(cell_index) = (*bordering_cells.at(cell_index)).type;

	for (int cell_index = 0; cell_index != number_of_bordering_cells; cell_index++){
		if (bordering_cell_types.at(cell_index) == standard_cell_type)
			is_faulty_cell = false;
	}

	return is_faulty_cell;
}

//Calculate the variance from the history vector.
double variance_history(const vector<double> &test_vector){
	//Store the length of the history vector.
	int length_vector = (test_vector).size();
	//Sum the history as well as the square of the history.
	double sum_history = 0.0;
	double sum_square_history = 0.0;
	for (int i = 0; i != length_vector; i++){
		sum_history += test_vector.at(i);
		sum_square_history += test_vector.at(i) * test_vector.at(i);
	}

	//Store the mean and mean of squares.
	double mean_history = sum_history / length_vector;
	double mean_square_history = sum_square_history / length_vector;
	
	//Calculate the variance.
	double variance = mean_square_history - mean_history * mean_history;

	return variance;
}

//Find domains of spins (cell-cell type conflicts).
vector<int> find_domains_of_spin(const vector<double> &spin_vector){
	//Declare vect for mismatched_cells (cell-cell type conflicts).
	vector<int> mismatched_cells;
	//Store the number of spins in spin_vector.
	int number_of_spins = spin_vector.size();

	//Try to find type-type mismatches.
	for (int spin_index = 1; spin_index < number_of_spins + 1; spin_index++){
		//Store neighboring spins.
		double first_spin = spin_vector.at(spin_index - 1);
		double second_spin = spin_vector.at(spin_index % number_of_spins);

		//If there is a spin mismatch, store the corresponding index (first spin in a new domain).
		if (first_spin * second_spin > 0.0){
			mismatched_cells.push_back(spin_index % number_of_spins);
		}
	}

	return mismatched_cells;
}

//Check if test_cell is within list_of_cells.
bool check_for_cell_pointer_within_list(PCPCell* test_cell, vector<PCPCell*> &list_of_cells){
	//Length of list for search
	int length_list_of_cells = list_of_cells.size();
	//If there are no elements, return false.
	if (length_list_of_cells == 0)
		return false;

	//Store count of matches (times that test_cell appear in list of cells).
	int match_count = 0;
	for (int cell_index = 0; cell_index != length_list_of_cells; cell_index++){
		if ((*list_of_cells.at(cell_index)).id == (*test_cell).id){
			match_count++;
		}
	}

	//If there are any matches, store true.
	bool cell_match = false;
	if (match_count != 0)
		cell_match = true;

	return cell_match;
}