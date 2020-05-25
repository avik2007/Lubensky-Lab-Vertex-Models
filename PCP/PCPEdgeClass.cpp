/**	\file	PCPEdgeClass.cpp
*	Implementation of PCPEdge class.
*	Contains function definitions associated with PCPEdge objects.
*	Depends on base class Edge.
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
#include "PCPCellClass.h"
#include "PCPMegaCellClass.h"
#include "PCPEdgeClass.h"
#include "PCPTissueClass.h"
#include "PCPMegaTissueClass.h"
#include "MegaCellConstants.h"

using namespace std;

//Checks whether a edge has PCP.
bool PCPEdge::is_Edge_with_PCP(){
	bool differentiated_edge = false;
	//Cast cell1 and cell2 pointers to pcp pointers.
	PCPCell* pcp_cell1 = dynamic_cast<PCPCell*>(cell1);
	PCPCell* pcp_cell2 = dynamic_cast<PCPCell*>(cell2);

	//if either cell1 or cell2 points to a differentiated cell, set the bool to true.
	if ((*pcp_cell1).is_Cell_with_PCP() || (*pcp_cell2).is_Cell_with_PCP()){
		differentiated_edge = true;
	}
	return differentiated_edge;
}

//Assign type to an edge based on the cells that share it (allow increased tension between cells of the same type).
int PCPEdge::Define_Edge_Type(){
	int edge_type = 0;

	if ((*this).is_Edge_with_PCP()){
		edge_type = 1;
	}

	//return 0 if the edge is between cells with no PCP protein.
	//return 1 if the edge has PCP
	return edge_type;
}

//Returns the magnitude of the tension on the edge (this is implemented differently in the edge subclasses) 
// \note MUST DEFINE GET TENSION MAGNITUDE WITHIN SUBCLASS IF YOU REQUIRE FANCY MAGNITUDE->CALLED INDEPENDENTLY WITHIN CALCULATETISSUESTRESS
double PCPEdge::GetTension_Magnitude(){
	vector<double> polarity_protein = planarpolarity;
	double sum_of_protein_concentrations = 0.0; 		//used to identify if this edges is a precursor edge.
	double precursor_tension = GAMMA;  					//define precursor tension in the case of the edge having no protein.
	int ind = 0;
	double corrected_gamma; 							//used to store the tension for the edge
	double edge_length = (*this).GetLength();  			//.r is updated.
	int edge_type = (*this).Define_Edge_Type();			//edge_type deciphers edges with PCP from those without, and also discriminates based on cell type.

	//sum over polarity protein on this edge.
	for (ind = 0; ind != polarity_protein.size(); ind++){
		sum_of_protein_concentrations += polarity_protein[ind] / edge_length;
	}
	//If the edge has no PCP protein.
	if (edge_type == 0){
		corrected_gamma = precursor_tension;
	}
	//If the edge has PCP protein, but not between cells of conflicting types.
	else if (edge_type == 1){
		corrected_gamma = GAMMA0 - GAMMAC * sum_of_protein_concentrations;
	}

	return corrected_gamma;
}

//Returns the magnitude of tension on a vertex due to a virtual edge given the concentration for a virtual edge specified in Lubensky PLOS 2012.
double PCPEdge::GetTension_Magnitude_Virtual(Vertex* v, PCPEdge* edge1_cellalpha, PCPEdge* edge2_cellalpha, PCPEdge* edge1_cellbeta, PCPEdge* edge2_cellbeta){
	bool differentiated_edge = (*this).is_Edge_with_PCP();
	int edge_type = (*this).Define_Edge_Type();
	double corrected_gamma;
	double precursor_tension = GAMMA;

	if (differentiated_edge){
		//convention:
		//cellalpha corresponds cell1
		//cellbeta corresponds cell2
		int whichcell_edge1_cellalpha;
		int whichcell_edge2_cellalpha;
		int whichcell_edge1_cellbeta;
		int whichcell_edge2_cellbeta;
		
		//Make sure to assign cells correctly about the new edge.
		if (*cell1 == *(*edge1_cellalpha).cell1){
			whichcell_edge1_cellalpha = 0;
		}
		else if (*cell1 == *(*edge1_cellalpha).cell2){
			whichcell_edge1_cellalpha = 1;
		}

		if (*cell1 == *(*edge2_cellalpha).cell1){
			whichcell_edge2_cellalpha = 0;
		}
		else if (*cell1 == *(*edge2_cellalpha).cell2){
			whichcell_edge2_cellalpha = 1;
		}

		if (*cell2 == *(*edge1_cellbeta).cell1){
			whichcell_edge1_cellbeta = 0;
		}
		else if (*cell2 == *(*edge1_cellbeta).cell2){
			whichcell_edge1_cellbeta = 1;
		}

		if (*cell2 == *(*edge2_cellbeta).cell1){
			whichcell_edge2_cellbeta = 0;
		}
		else if (*cell2 == *(*edge2_cellbeta).cell2){
			whichcell_edge2_cellbeta = 1;
		}

		double edge1_cellalpha_length = (*edge1_cellalpha).GetLength();
		double edge2_cellalpha_length = (*edge2_cellalpha).GetLength();
		double edge1_cellbeta_length = (*edge1_cellbeta).GetLength();
		double edge2_cellbeta_length = (*edge2_cellbeta).GetLength();

		//Find sum of protein concentration that would be on newly formed edge.
		double calpha_proteinA = kT1 * min((*edge1_cellalpha).planarpolarity[2 * whichcell_edge1_cellalpha] / edge1_cellalpha_length, (*edge2_cellalpha).planarpolarity[2 * whichcell_edge2_cellalpha] / edge2_cellalpha_length);
		double calpha_proteinB = kT1 * min((*edge1_cellalpha).planarpolarity[2 * whichcell_edge1_cellalpha + 1] / edge1_cellalpha_length, (*edge2_cellalpha).planarpolarity[2 * whichcell_edge2_cellalpha + 1] / edge2_cellalpha_length);
		double cbeta_proteinA = kT1 * min((*edge1_cellbeta).planarpolarity[2 * whichcell_edge1_cellbeta] / edge1_cellbeta_length, (*edge2_cellbeta).planarpolarity[2 * whichcell_edge2_cellbeta] / edge2_cellbeta_length);
		double cbeta_proteinB = kT1 * min((*edge1_cellbeta).planarpolarity[2 * whichcell_edge1_cellbeta + 1] / edge1_cellbeta_length, (*edge2_cellbeta).planarpolarity[2 * whichcell_edge2_cellbeta + 1] / edge2_cellbeta_length);
		double sum_of_protein_concentrations = calpha_proteinA + calpha_proteinB + cbeta_proteinA + cbeta_proteinB;
		//If the edge has PCP protein,
		if (edge_type == 1){
			corrected_gamma = GAMMA0 - GAMMAC * sum_of_protein_concentrations;
		}
	}
	else{
		if (edge_type == 0){
			corrected_gamma = precursor_tension;
		}
	}

	return corrected_gamma;
}

//Returns the tension on a vertex due to a virtual edge given the concentration for a virtual edge specified in Lubensky PLOS 2012.
Vect PCPEdge::GetTension_ForceOnVertex_Virtual(Vertex* v, PCPEdge* edge1_cellalpha, PCPEdge* edge2_cellalpha, PCPEdge* edge1_cellbeta, PCPEdge* edge2_cellbeta){
	double angle = (*this).GetAngle(v);
	Vect origin(0.0, 0.0);																										//define the origin
	double mag = (*this).GetTension_Magnitude_Virtual(v, edge1_cellalpha, edge2_cellalpha, edge1_cellbeta, edge2_cellbeta);		//Get the magnitude of the tension
	return Vect(origin, angle, mag);																							//same magnitude for all tension
}