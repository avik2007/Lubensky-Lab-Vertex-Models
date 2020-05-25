//! \file PCPEvolutionFunction.h
//! This file contains declaration and implementation of function(s) that evaluate time evolution of amount of PCP on edges (i.e., the RHS of the ODE).
//! This function implements the GSL root finding functions for the numerical ODE solver adapted from ASCHER, RUUTH, WETTON

// include guard
#ifndef __PCPEVOLUTIONFUNCTION_H_INCLUDED__
#define __PCPEVOLUTIONFUNCTION_H_INCLUDED__

//Forward declared dependences
#include <cmath>
#include <list>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <random>
#include <omp.h>

//Base class functions
#include "VectClass.h"
#include "TissueClass.h"
#include "EdgeClass.h"
#include "VertexClass.h"
#include "Constants.h"
#include "CellClass.h"

//PCP sub-class files
#include "PCPTissueClass.h"
#include "PCPCellClass.h"
#include "PCPEdgeClass.h"
#include "PCPConstants.h"
#include "PCPMegaTissueClass.h"
#include "PCPMegaCellClass.h"
#include "PCPMegaPoissonTissueClass.h"
#include "PCPPoissonTissueClass.h"

//gsl header files
#include "gsl/gsl_errno.h"			//Needed for GSL_SUCCESS return code
#include "gsl/gsl_matrix.h"			//GSL MATRIX FUNCTIONS
#include "gsl/gsl_vector.h"
#include "gsl/gsl_blas.h"			//GSL MATRIX LINEAR ALGEBRA
#include "gsl/gsl_multiroots.h"		//GSL MULTIDIMENSIONAL ROOTFINDING

extern TISSUE_TYPE tissue;				//Declare 'extern' so it can be used by the global-scope function Odeiv2System_RHS(). Definition of tissue is in main.cpp.

//rparams stores the necessary information for the two numerical solvers
struct rparams
{
	//rhs of dynamics (not including the next step terms)
	vector<double> rhs_evaluation;
	//the initial state of the system.
	vector<double> y_in;
	//edge lengths (in order corresponding to y_in)
	vector<double> length_edges;
	//Coefficients for the variable step algorithm for the next step terms.
	double alphastepforward;
	double gammastepforward;
	//The step size.
	double delta_t;
	//The V_matrix that couples edges within the same cell.
	vector<vector<double>> V_matrix;
};

//This function is used to implement the third order imex scheme cell by cell.
//Consult IMPLICIT-EXPLICIT METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by ASCHER, RUUTH, WETTON  (page 807)
//Consult VARIABLE STEP-SIZE IMPLICIT-EXPLICIT LINEAR MULTISTEP METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by WANG, RUUTH (page 6)
int numerical_solve_third_order(const gsl_vector * x, void *params, gsl_vector * f)
{
	//Save the coefficients that take into account variable step size.
	double alpha3n = ((struct rparams *) params)->alphastepforward;
	double gamma3n = ((struct rparams *) params)->gammastepforward;
	//The RHS of the evolution for third-order scheme.
	vector<double> rhs_evaluation = ((struct rparams *) params)->rhs_evaluation;
	//Also, save the current variables (y_in).
	vector<double> y_in = ((struct rparams *) params)->y_in;
	//Lengths of edges for the given cell.
	vector<double> length_edges = ((struct rparams *) params)->length_edges;
	//Step size.
	double delta_t = ((struct rparams *) params)->delta_t;
	//number of edges in the cell.
	int number_of_edges = length_edges.size();
	//Store submatrix of larger V matrix that is for this particular cell.
	vector<vector<double>> V_matrix = ((struct rparams *) params)->V_matrix;

	//Store each term of the numerical solution in the temporary variable f_term.
	double f_term = 0.0;

	for (int edge_index = 0; edge_index < number_of_edges; edge_index++){
		f_term = 0.0;
		f_term = alpha3n * alphaP / delta_t * exp(-y_in.at(edge_index)) * gsl_vector_get(x, edge_index);
		for (int secondary_edge_index = 0; secondary_edge_index < number_of_edges; secondary_edge_index++){
			f_term += gamma3n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * te * V_matrix[edge_index][secondary_edge_index] * gsl_vector_get(x, secondary_edge_index);
			if (max_concentration != 0.0){
				//add cmax term.
				f_term += gamma3n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * te * V_matrix[edge_index][secondary_edge_index] * log(max_concentration - exp(-gsl_vector_get(x, secondary_edge_index)) / length_edges.at(secondary_edge_index));
			}
		}
		f_term -= rhs_evaluation.at(edge_index);
		gsl_vector_set(f, edge_index, f_term);
	}

	return GSL_SUCCESS;
}

//This function is used to implement the third order imex scheme cell by cell.
//JACOBIAN
//Consult IMPLICIT-EXPLICIT METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by ASCHER, RUUTH, WETTON  (page 801-802)
int numerical_solve_third_order_df(const gsl_vector * x, void *params, gsl_matrix * J){
	//Save the coefficients that take into account variable step size.
	double alpha3n = ((struct rparams *) params)->alphastepforward;
	double gamma3n = ((struct rparams *) params)->gammastepforward;
	//The RHS of the evolution for third-order scheme.
	vector<double> rhs_evaluation = ((struct rparams *) params)->rhs_evaluation;
	//Also, save the current variables (y_in).
	vector<double> y_in = ((struct rparams *) params)->y_in;
	//Lengths of edges for the given cell.
	vector<double> length_edges = ((struct rparams *) params)->length_edges;
	//Step size.
	double delta_t = ((struct rparams *) params)->delta_t;
	int number_of_edges = length_edges.size();
	vector<vector<double>> V_matrix = ((struct rparams *) params)->V_matrix;
	double matrix_val = 0.0;

	//Build the offdiagonal elements (and put in initial piece of diagonal elements).
	for (int row_edge_index = 0; row_edge_index < number_of_edges; row_edge_index++){
		for (int column_edge_index = 0; column_edge_index < number_of_edges; column_edge_index++){
			matrix_val = gamma3n * exp(gsl_vector_get(x, row_edge_index) - y_in.at(row_edge_index)) * length_edges.at(row_edge_index) * te * V_matrix[row_edge_index][column_edge_index];
			if (max_concentration != 0.0){
				//Add cmax term.
				matrix_val += gamma3n * exp(gsl_vector_get(x, row_edge_index) - y_in.at(row_edge_index)) * length_edges.at(row_edge_index) * te * V_matrix[row_edge_index][column_edge_index] * exp(-gsl_vector_get(x, column_edge_index)) / (max_concentration * length_edges.at(column_edge_index) - exp(-gsl_vector_get(x, column_edge_index)));
			}
			gsl_matrix_set(J, row_edge_index, column_edge_index, matrix_val);
		}
	}

	//Now, add the additional terms along the diagonal term.
	for (int edge_index = 0; edge_index < number_of_edges; edge_index++){
		gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + alpha3n * alphaP / delta_t * exp(-y_in.at(edge_index)));
		for (int second_edge_index = 0; second_edge_index != number_of_edges; second_edge_index++){
			gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + gamma3n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * (te * V_matrix[edge_index][second_edge_index] * gsl_vector_get(x, second_edge_index)));
			if (max_concentration != 0.0){
				//add cmax term.
				gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + gamma3n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * (te * V_matrix[edge_index][second_edge_index] * log(max_concentration - exp(-gsl_vector_get(x, second_edge_index)) / length_edges.at(second_edge_index))));
			}
		}
	}

	return GSL_SUCCESS;
}

//Third-order envelope function
int numerical_solve_third_order_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J){
	numerical_solve_third_order(x, params, f);
	numerical_solve_third_order_df(x, params, J);

	return GSL_SUCCESS;
}

//This function is used to implement the second order imex scheme cell by cell.
//Consult IMPLICIT-EXPLICIT METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by ASCHER, RUUTH, WETTON  (page 807)
//Consult VARIABLE STEP-SIZE IMPLICIT-EXPLICIT LINEAR MULTISTEP METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by WANG, RUUTH (page 6)
int numerical_solve_second_order(const gsl_vector * x, void *params, gsl_vector * f)
{
	//Save the coefficients that take into account variable step size.
	double alpha2n = ((struct rparams *) params)->alphastepforward;
	double gamma2n = ((struct rparams *) params)->gammastepforward;
	//The RHS of the evolution for third-order scheme.
	vector<double> rhs_evaluation = ((struct rparams *) params)->rhs_evaluation;
	//Also, save the current variables (y_in).
	vector<double> y_in = ((struct rparams *) params)->y_in;
	//Lengths of edges for the given cell.
	vector<double> length_edges = ((struct rparams *) params)->length_edges;
	//Step size.
	double delta_t = ((struct rparams *) params)->delta_t;
	//number of edges in the cell.
	int number_of_edges = length_edges.size();
	//Store the matrix for the implicit terms.
	vector<vector<double>> V_matrix = ((struct rparams *) params)->V_matrix;

	//Store each term of the numerical solution in the temporary variable f_term.
	double f_term = 0.0;

	for (int edge_index = 0; edge_index < number_of_edges; edge_index++){
		f_term = 0.0;
		f_term = alpha2n * alphaP / delta_t * exp(-y_in.at(edge_index)) * gsl_vector_get(x, edge_index);
		for (int secondary_edge_index = 0; secondary_edge_index < number_of_edges; secondary_edge_index++){
			f_term += gamma2n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * te * V_matrix[edge_index][secondary_edge_index] * gsl_vector_get(x, secondary_edge_index);
			if (max_concentration != 0.0){
				//add cmax term.
				f_term += gamma2n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * te * V_matrix[edge_index][secondary_edge_index] * log(max_concentration - exp(-gsl_vector_get(x, secondary_edge_index)) / length_edges.at(secondary_edge_index));
			}
		}
		f_term -= rhs_evaluation.at(edge_index);
		gsl_vector_set(f, edge_index, f_term);
	}

	return GSL_SUCCESS;
}

//This function is used to implement the second order imex scheme cell by cell.
//JACOBIAN
//Consult IMPLICIT-EXPLICIT METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by ASCHER, RUUTH, WETTON  (page 801-802)
int numerical_solve_second_order_df(const gsl_vector * x, void *params, gsl_matrix * J){
	//Save the coefficients that take into account variable step size.
	double alpha2n = ((struct rparams *) params)->alphastepforward;
	double gamma2n = ((struct rparams *) params)->gammastepforward;
	//The RHS of the evolution for third-order scheme.
	vector<double> rhs_evaluation = ((struct rparams *) params)->rhs_evaluation;
	//Also, save the current variables (y_in).
	vector<double> y_in = ((struct rparams *) params)->y_in;
	//Lengths of edges for the given cell.
	vector<double> length_edges = ((struct rparams *) params)->length_edges;
	//Step size.
	double delta_t = ((struct rparams *) params)->delta_t;
	int number_of_edges = length_edges.size();
	//Store submatrix of larger V matrix that is for this particular cell.
	vector<vector<double>> V_matrix = ((struct rparams *) params)->V_matrix;
	double matrix_val = 0.0;

	//Build the offdiagonal elements (and put in initial piece of diagonal elements).
	for (int row_edge_index = 0; row_edge_index < number_of_edges; row_edge_index++){
		for (int column_edge_index = 0; column_edge_index < number_of_edges; column_edge_index++){
			matrix_val = gamma2n * exp(gsl_vector_get(x, row_edge_index) - y_in.at(row_edge_index)) * length_edges.at(row_edge_index) * te * V_matrix[row_edge_index][column_edge_index];
			if (max_concentration != 0.0){
				//Add cmax term.
				matrix_val += gamma2n * exp(gsl_vector_get(x, row_edge_index) - y_in.at(row_edge_index)) * length_edges.at(row_edge_index) * te * V_matrix[row_edge_index][column_edge_index] * exp(-gsl_vector_get(x, column_edge_index)) / (max_concentration * length_edges.at(column_edge_index) - exp(-gsl_vector_get(x, column_edge_index)));
			}
			gsl_matrix_set(J, row_edge_index, column_edge_index, matrix_val);
		}
	}

	//Now, add the additional terms along the diagonal term.
	for (int edge_index = 0; edge_index < number_of_edges; edge_index++){
		gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + alpha2n * alphaP / delta_t * exp(-y_in.at(edge_index)));
		for (int second_edge_index = 0; second_edge_index != number_of_edges; second_edge_index++){
			gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + gamma2n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * (te * V_matrix[edge_index][second_edge_index] * gsl_vector_get(x, second_edge_index)));
			if (max_concentration != 0.0){
				//add cmax term.
				gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + gamma2n * exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * (te * V_matrix[edge_index][second_edge_index] * log(max_concentration - exp(-gsl_vector_get(x, second_edge_index)) / length_edges.at(second_edge_index))));
			}
		}
	}

	return GSL_SUCCESS;
}

//Second-order envelope function
int numerical_solve_second_order_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J){
	numerical_solve_second_order(x, params, f);
	numerical_solve_second_order_df(x, params, J);

	return GSL_SUCCESS;
}

//This function is used to implement the first order imex scheme cell by cell.
//Consult IMPLICIT-EXPLICIT METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by ASCHER, RUUTH, WETTON  (page 807)
//Consult VARIABLE STEP-SIZE IMPLICIT-EXPLICIT LINEAR MULTISTEP METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by WANG, RUUTH (page 6)
int numerical_solve_first_order(const gsl_vector * x, void *params, gsl_vector * f){
	//The RHS of the evolution for third-order scheme.
	vector<double> rhs_evaluation = ((struct rparams *) params)->rhs_evaluation;
	//Also, save the current variables (y_in).
	vector<double> y_in = ((struct rparams *) params)->y_in;
	//Lengths of edges for the given cell.
	vector<double> length_edges = ((struct rparams *) params)->length_edges;
	//Step size.
	double delta_t = ((struct rparams *) params)->delta_t;
	//number of edges in the cell.
	int number_of_edges = length_edges.size();
	//Store the matrix for the implicit terms.
	vector<vector<double>> V_matrix = ((struct rparams *) params)->V_matrix;

	//Store each term of the numerical solution in the temporary variable f_term.
	double f_term = 0.0;
	for (int edge_index = 0; edge_index < number_of_edges; edge_index++){
		f_term = 0.0;
		f_term = alphaP / delta_t * exp(-y_in.at(edge_index)) * gsl_vector_get(x, edge_index);
		for (int secondary_edge_index = 0; secondary_edge_index < number_of_edges; secondary_edge_index++){
			f_term += exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * te * V_matrix[edge_index][secondary_edge_index] * gsl_vector_get(x, secondary_edge_index);
			if (max_concentration != 0.0){
				//Add cmax term.
				f_term += exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * te * V_matrix[edge_index][secondary_edge_index] * log(max_concentration - exp(-gsl_vector_get(x, secondary_edge_index)) / length_edges.at(secondary_edge_index));
			}
		}
		f_term -= rhs_evaluation.at(edge_index);
		gsl_vector_set(f, edge_index, f_term);
	}

	return GSL_SUCCESS;
}

//This function is used to implement the first order imex scheme cell by cell.
//JACOBIAN
//Consult IMPLICIT-EXPLICIT METHODS FOR TIME-DEPENDENT PARTIAL DIFFERENTIAL EQUATIONS by ASCHER, RUUTH, WETTON  (page 801-802)
int numerical_solve_first_order_df(const gsl_vector * x, void *params, gsl_matrix * J){
	//Also, save the current variables (y_in).
	vector<double> y_in = ((struct rparams *) params)->y_in;
	//Lengths of edges for the given cell.
	vector<double> length_edges = ((struct rparams *) params)->length_edges;
	//Step size.
	double delta_t = ((struct rparams *) params)->delta_t;
	int number_of_edges = length_edges.size();
	//Store submatrix of larger V matrix that is for this particular cell.
	vector<vector<double>> V_matrix = ((struct rparams *) params)->V_matrix;
	double matrix_val = 0.0;

	//Build the offdiagonal elements (and put in initial piece of diagonal elements).
	for (int row_edge_index = 0; row_edge_index < number_of_edges; row_edge_index++){
		for (int column_edge_index = 0; column_edge_index < number_of_edges; column_edge_index++){
			matrix_val = exp(gsl_vector_get(x, row_edge_index) - y_in.at(row_edge_index)) * length_edges.at(row_edge_index) * te * V_matrix[row_edge_index][column_edge_index];
			if (max_concentration != 0.0){
				//Add cmax term.
				matrix_val += exp(gsl_vector_get(x, row_edge_index) - y_in.at(row_edge_index)) * length_edges.at(row_edge_index) * te * V_matrix[row_edge_index][column_edge_index] * exp(-gsl_vector_get(x, column_edge_index)) / (max_concentration * length_edges.at(column_edge_index) - exp(-gsl_vector_get(x, column_edge_index)));
			}
			gsl_matrix_set(J, row_edge_index, column_edge_index, matrix_val);
		}
	}

	//Now, add the additional terms along the diagonal term.
	for (int edge_index = 0; edge_index < number_of_edges; edge_index++){
		gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + alphaP / delta_t * exp(-y_in.at(edge_index)));
		for (int second_edge_index = 0; second_edge_index != number_of_edges; second_edge_index++){
			gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * (te * V_matrix[edge_index][second_edge_index] * gsl_vector_get(x, second_edge_index)));
			if (max_concentration != 0.0){
				//add cmax term.
				gsl_matrix_set(J, edge_index, edge_index, gsl_matrix_get(J, edge_index, edge_index) + exp(gsl_vector_get(x, edge_index) - y_in.at(edge_index)) * length_edges.at(edge_index) * (te * V_matrix[edge_index][second_edge_index] * log(max_concentration - exp(-gsl_vector_get(x, second_edge_index)) / length_edges.at(second_edge_index))));
			}
		}
	}

	return GSL_SUCCESS;
}

//First-order envelope function
int numerical_solve_first_order_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J){
	numerical_solve_first_order(x, params, f);
	numerical_solve_first_order_df(x, params, J);

	return GSL_SUCCESS;
}

//This is the function for initializing the U-matrix for all of the ode solving functions (for calculation of explicit terms).
//Consult LUBENSKY PLOS 2012 paper for definition of U-matrix (supp info)
//This is no longer used in the imex solvers because this requires the initialization of a matrix with many zeros.
void Initialize_U_Matrix(gsl_matrix * U, void* params, const vector<int> &start_index_cells){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	//Store cell id 2 index.
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	//The dimension of the matrix that is of only one type is half the size of the matrix.
	int type_dimension = dim_proteins / 2.;
	//Initialize a vector to store the length list of edges.
	vector<double> edge_length_list;
	//Initialize a double to store the sum over edge lengths for a given cell.
	double sum_over_edge_lengths = 0.0;
	int start_index = 0;
	//number of edges in the given cell.
	int number_of_edges;
	//Helps to identify the index associated with the concentration on the other side of an edge.
	Cell* column_cell_pointer;
	vector<PCPEdge*> column_cell_ordered_edges;
	int column_cell_id;
	int column_cell_index;
	int cell_index = 0;
	int cell_id;
	//Helps to identify the index associated with the concentration on the other side of an edge.
	int beta_index;
	//Used for only adding j1, j2 contributions to neighboring cell elements.
	double current_element_value;
	//Used for finding a neighboring edge.
	int next_index;
	vector<PCPEdge*> edge_list_for_this_cell;
	PCPCell* pcpcit;

	//Loop through all cells (as below for V) in order to build the U matrix and the vector of concentrations.
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		pcpcit = dynamic_cast<PCPCell*>(*cit);
		//only proceed with calculations if the current cell is differentiated.
		if ((*pcpcit).is_Cell_with_PCP()){
			//Cast the cell iterator to a PCPCell pointer.
			PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
			//Store id, index for the given cell
			cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			//Find the start index in the U matrix for this given cell.
			start_index = start_index_cells.at(cell_index);
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
				if ((*(*edge_list_for_this_cell[i]).cell1) == *(*cit)){
					//Store the index for the cell on the other side of the edge.
					other_cell_index[i] = cell_id2index.at((*(*edge_list_for_this_cell[i]).cell2).id);
				}
				else{
					//Store the index for the cell on the other side of the edge.
					other_cell_index[i] = cell_id2index.at((*(*edge_list_for_this_cell[i]).cell1).id);
				}
			}
			//compute sum over edge lengths.
			sum_over_edge_lengths = 0.0;
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				sum_over_edge_lengths += edge_length_list[edge_index];
			}
			//Now access specific elements of U. Keep track of current absolute position in matrix by storing using start index.
			for (int edge_row_index = 0; edge_row_index != number_of_edges; edge_row_index++){
				for (int edge_column_index = 0; edge_column_index != number_of_edges; edge_column_index++){
					//For the edge referred to by the current column index, look to the other side of the edge.
					if ((*(*edge_list_for_this_cell[edge_column_index]).cell1) == *(*cit)){
						column_cell_pointer = (*edge_list_for_this_cell[edge_column_index]).cell2;
					}
					else {
						column_cell_pointer = (*edge_list_for_this_cell[edge_column_index]).cell1;
					}
					//Cast the cell on the other side of the current edge to a PCPCell pointer.
					PCPCell* pcp_column_cell_ptr = dynamic_cast<PCPCell*>(column_cell_pointer);
					column_cell_id = (*column_cell_pointer).id;
					column_cell_index = cell_id2index.at(column_cell_id);
					//For the cell on the other side of the edge, order the edges of that cell and find which edge corresponds to this edge.
					column_cell_ordered_edges.resize((*column_cell_pointer).edges.size());
					column_cell_ordered_edges = (*pcp_column_cell_ptr).FindOrderedEdges();
					//Search for which index corresponds to this edge in the neighboring cell.
					for (int i = 0; i != column_cell_ordered_edges.size(); i++){
						if (*column_cell_ordered_edges[i] == *edge_list_for_this_cell[edge_column_index]){
							//Store the index for the cell on the other side of the given edge.
							beta_index = i;
						}
					}
					if (edge_column_index != edge_row_index){
						//Consult equation S1, S2, S10, S11 of PLOS 2012 paper to see where this element comes from. Note that to compute the the term from V's for a given concentration, one sums over product of columns of V with the corresponding state vector.
						//Now use which cell to access the appropriate elements.
						//check whether this cell index is within the range of row index
						gsl_matrix_set(U, edge_row_index + start_index, edge_column_index + start_index, -epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U, edge_row_index + start_index, type_dimension + edge_column_index + start_index, 2 * jtwo / sum_over_edge_lengths + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths);
						gsl_matrix_set(U, type_dimension + edge_row_index + start_index, edge_column_index + start_index, 2 * jtwo / sum_over_edge_lengths + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths);
						gsl_matrix_set(U, type_dimension + edge_row_index + start_index, type_dimension + edge_column_index + start_index, -epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						if ((*pcp_column_cell_ptr).is_Cell_with_PCP()){
							gsl_matrix_set(U, edge_row_index + start_index, type_dimension + start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths);
							gsl_matrix_set(U, edge_row_index + start_index, start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths);
							gsl_matrix_set(U, type_dimension + edge_row_index + start_index, type_dimension + start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths);
							gsl_matrix_set(U, type_dimension + edge_row_index + start_index, start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths);
						}
						//Find the neighboring edge at the lower cell index, if you hit an index of -1, loop to highest index.
						next_index = (edge_row_index - 1) % number_of_edges;
						if (next_index == -1){
							next_index = number_of_edges - 1;
						}
						if ((edge_column_index == (edge_row_index + 1) % number_of_edges) || (edge_column_index == next_index)){
							current_element_value = gsl_matrix_get(U, edge_row_index + start_index, edge_column_index + start_index);
							gsl_matrix_set(U, edge_row_index + start_index, edge_column_index + start_index, current_element_value + jone / edge_length_list[edge_row_index]);
							current_element_value = gsl_matrix_get(U, edge_row_index + start_index, type_dimension + edge_column_index + start_index);
							gsl_matrix_set(U, edge_row_index + start_index, type_dimension + edge_column_index + start_index, current_element_value - jtwo / edge_length_list[edge_row_index]);
							current_element_value = gsl_matrix_get(U, type_dimension + edge_row_index + start_index, edge_column_index + start_index);
							gsl_matrix_set(U, type_dimension + edge_row_index + start_index, edge_column_index + start_index, current_element_value - jtwo / edge_length_list[edge_row_index]);
							current_element_value = gsl_matrix_get(U, type_dimension + edge_row_index + start_index, type_dimension + edge_column_index + start_index);
							gsl_matrix_set(U, type_dimension + edge_row_index + start_index, type_dimension + edge_column_index + start_index, current_element_value + jone / edge_length_list[edge_row_index]);
						}
					}
					else {
						gsl_matrix_set(U, edge_row_index + start_index, edge_column_index + start_index, epsnaught - jone / edge_length_list[edge_row_index] * 2 + jtwo / edge_length_list[edge_column_index] * 2 - epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U, edge_row_index + start_index, type_dimension + edge_column_index + start_index, -epsab + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths + 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U, type_dimension + edge_row_index + start_index, type_dimension + edge_column_index + start_index, epsnaught - jone / edge_length_list[edge_row_index] * 2 + jtwo / edge_length_list[edge_column_index] * 2 - epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U, type_dimension + edge_row_index + start_index, edge_column_index + start_index, -epsab + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths + 2 * jtwo / sum_over_edge_lengths);
						//if this edge is in contact with another differentiated cell
						if ((*pcp_column_cell_ptr).is_Cell_with_PCP()){
							gsl_matrix_set(U, edge_row_index + start_index, type_dimension + start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths + epsone);
							gsl_matrix_set(U, edge_row_index + start_index, start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths + epstwo);
							gsl_matrix_set(U, type_dimension + edge_row_index + start_index, type_dimension + start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths + epstwo);
							gsl_matrix_set(U, type_dimension + edge_row_index + start_index, start_index_cells.at(other_cell_index[edge_column_index]) + beta_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths + epsone);
						}
					}
				}
			}
		}
	}
}

//Instead of initializing a matrix for U terms (that is very sparse), define the evolution vector U separately (to avoid defining large matrix).
//This involves coupling within a cell and to nearest-neighboring cells.
void Define_Evolution_Vector_U(gsl_vector* evolution_vectorU, gsl_vector* concentration_vector, void* params, const vector<int> &start_index_cells){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	//Store cell id 2 index.
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	//The dimension of the matrix that is of only one type is half the size of the matrix.
	int type_dimension = dim_proteins / 2.;
	//Initialize a vector to store the length list of edges.
	vector<double> edge_length_list;
	//Initialize a double to store the sum over edge lengths for a given cell.
	double sum_over_edge_lengths = 0.0;
	int start_index = 0;
	//number of edges in the given cell.
	int number_of_edges;
	//Helps to identify the index associated with the concentration on the other side of an edge.
	Cell* column_cell_pointer;
	vector<PCPEdge*> column_cell_ordered_edges;
	int column_cell_id;
	int column_cell_index;
	int cell_index = 0;
	int cell_id;
	//Helps to identify the index associated with the concentration on the other side of an edge.
	int beta_index;
	//Used for only adding j1, j2 contributions to neighboring cell elements.
	double current_element_value;
	//Used for finding a neighboring edge.
	int next_index;
	vector<PCPEdge*> edge_list_for_this_cell;
	PCPCell* pcpcit;

	//Allocate the matrices and vectors for calculating UijCj (coupling between different concentrations with energetic couplings).
	int max_number_of_edges = 100;
	int max_number_of_rows = max_number_of_edges * 2;
	int max_number_of_columns = max_number_of_edges * 4;
	gsl_matrix* U_matrix_by_cell = gsl_matrix_calloc(max_number_of_rows, max_number_of_columns);
	gsl_vector* concentration_vector_by_cell = gsl_vector_calloc(max_number_of_columns);
	gsl_vector* evolution_SUBVECTOR = gsl_vector_calloc(max_number_of_rows);

	//Loop through all cells in order to build the U matrix.
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//Store id, index for the given cell
			cell_id = (*pcpcit).id;
			cell_index = cell_id2index.at(cell_id);
			//Find the start index in the U matrix for this given cell.
			start_index = start_index_cells.at(cell_index);
			edge_list_for_this_cell.resize((*(*cit)).edges.size());
			edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			number_of_edges = edge_list_for_this_cell.size();

			//Store U matrix sizes 
			//(every edge concentration A and B is coupled to every other concentration within the cells and to corresponding concentrations on the other side of the edges).
			int U_matrix_by_columns_size = number_of_edges * 4;
			int U_matrix_by_rows_size = number_of_edges * 2;
			//Store vectors for all concentrations (0 -> number_of_edges-1: concentrations A within the cell; number_of_edges -> 2*number_of_edges-1: concentrations B within the cell; 2*number_of_edges -> 3*number_of_edges-1: concentrations A neighboring the cell; 3*number_of_edges -> 4*number_of_edges-1: concentrations B neighboring the cell)
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(concentration_vector_by_cell, edge_index, gsl_vector_get(concentration_vector, start_index + edge_index));
				gsl_vector_set(concentration_vector_by_cell, number_of_edges + edge_index, gsl_vector_get(concentration_vector, start_index + type_dimension + edge_index));
			}

			//Build the list of edge lengths for this cell.
			vector<double> edge_length_list;
			edge_length_list.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
			}
			//This is used to store the index for the cell on the other side of all edges.
			vector<int> other_cell_index;
			other_cell_index.resize(number_of_edges);
			//build the vector of protein concentrations (following the same conventions as the U matrix)
			for (int i = 0; i != number_of_edges; i++){
				if ((*(*edge_list_for_this_cell[i]).cell1) == *(*cit)){
					//Store the index for the cell on the other side of the edge.
					other_cell_index[i] = cell_id2index.at((*(*edge_list_for_this_cell[i]).cell2).id);
				}
				else{
					//Store the index for the cell on the other side of the edge.
					other_cell_index[i] = cell_id2index.at((*(*edge_list_for_this_cell[i]).cell1).id);
				}
			}
			//compute sum over edge lengths.
			sum_over_edge_lengths = 0.0;
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				sum_over_edge_lengths += edge_length_list[edge_index];
			}
			//Now access specific elements of U. Keep track of current absolute position in matrix by storing using start index.
			for (int edge_row_index = 0; edge_row_index != number_of_edges; edge_row_index++){
				for (int edge_column_index = 0; edge_column_index != number_of_edges; edge_column_index++){
					//For the edge referred to by the current column index, look to the other side of the edge.
					if ((*(*edge_list_for_this_cell[edge_column_index]).cell1) == *(*cit)){
						column_cell_pointer = (*edge_list_for_this_cell[edge_column_index]).cell2;
					}
					else {
						column_cell_pointer = (*edge_list_for_this_cell[edge_column_index]).cell1;
					}
					//Cast the cell on the other side of the current edge to a PCPCell pointer.
					PCPCell* pcp_column_cell_ptr = dynamic_cast<PCPCell*>(column_cell_pointer);
					column_cell_id = (*column_cell_pointer).id;
					column_cell_index = cell_id2index.at(column_cell_id);
					//For the cell on the other side of the edge, order the edges of that cell and find which edge corresponds to this edge.
					column_cell_ordered_edges.resize((*column_cell_pointer).edges.size());
					column_cell_ordered_edges = (*pcp_column_cell_ptr).FindOrderedEdges();
					//Search for which index corresponds to this edge in the neighboring cell.
					for (int i = 0; i != column_cell_ordered_edges.size(); i++){
						if (*column_cell_ordered_edges[i] == *edge_list_for_this_cell[edge_column_index]){
							//Store the index for the cell on the other side of the given edge.
							beta_index = i;
						}
					}
					//If the cell in contact with this edge is also differentiated, store the concentration on the other side of the edge.
					if ((*pcp_column_cell_ptr).is_Cell_with_PCP()){
						gsl_vector_set(concentration_vector_by_cell, 2 * number_of_edges + edge_column_index, gsl_vector_get(concentration_vector, start_index_cells.at(other_cell_index[edge_column_index]) + beta_index));
						gsl_vector_set(concentration_vector_by_cell, 3 * number_of_edges + edge_column_index, gsl_vector_get(concentration_vector, type_dimension + start_index_cells.at(other_cell_index[edge_column_index]) + beta_index));
					}
					if (edge_column_index != edge_row_index){
						//Consult equation S1, S2, S10, S11 of PLOS 2012 paper to see where this element comes from. Note that to compute the the term from V's for a given concentration, one sums over product of columns of V with the corresponding state vector.
						//Now use which cell to access the appropriate elements.
						//check whether this cell index is within the range of row index
						gsl_matrix_set(U_matrix_by_cell, edge_row_index, edge_column_index, -epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U_matrix_by_cell, edge_row_index, number_of_edges + edge_column_index, 2 * jtwo / sum_over_edge_lengths + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths);
						gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, edge_column_index, 2 * jtwo / sum_over_edge_lengths + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths);
						gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, number_of_edges + edge_column_index, -epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						//if this edge is in contact with another differentiated cell, store the coupling to the neighboring concentration.
						if ((*pcp_column_cell_ptr).is_Cell_with_PCP()){
							gsl_matrix_set(U_matrix_by_cell, edge_row_index, 3 * number_of_edges + edge_column_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths);
							gsl_matrix_set(U_matrix_by_cell, edge_row_index, 2 * number_of_edges + edge_column_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths);
							gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, 3 * number_of_edges + edge_column_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths);
							gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, 2 * number_of_edges + edge_column_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths);
						}
						//Find the neighboring edge (within the same cell) at the lower cell index, if you hit an index of -1, loop to highest_index.
						//This penalizes sharp gradients in concentration.
						next_index = (edge_row_index - 1) % number_of_edges;
						if (next_index == -1){
							next_index = number_of_edges - 1;
						}
						if ((edge_column_index == (edge_row_index + 1) % number_of_edges) || (edge_column_index == next_index)){
							current_element_value = gsl_matrix_get(U_matrix_by_cell, edge_row_index, edge_column_index);
							gsl_matrix_set(U_matrix_by_cell, edge_row_index, edge_column_index, current_element_value + jone / edge_length_list[edge_row_index]);
							current_element_value = gsl_matrix_get(U_matrix_by_cell, edge_row_index, number_of_edges + edge_column_index);
							gsl_matrix_set(U_matrix_by_cell, edge_row_index, number_of_edges + edge_column_index, current_element_value - jtwo / edge_length_list[edge_row_index]);
							current_element_value = gsl_matrix_get(U_matrix_by_cell, number_of_edges + edge_row_index, edge_column_index);
							gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, edge_column_index, current_element_value - jtwo / edge_length_list[edge_row_index]);
							current_element_value = gsl_matrix_get(U_matrix_by_cell, number_of_edges + edge_row_index, number_of_edges + edge_column_index);
							gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, number_of_edges + edge_column_index, current_element_value + jone / edge_length_list[edge_row_index]);
						}
					}
					else {
						gsl_matrix_set(U_matrix_by_cell, edge_row_index, edge_column_index, epsnaught - jone / edge_length_list[edge_row_index] * 2 + jtwo / edge_length_list[edge_column_index] * 2 - epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U_matrix_by_cell, edge_row_index, number_of_edges + edge_column_index, -epsab + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths + 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, number_of_edges + edge_column_index, epsnaught - jone / edge_length_list[edge_row_index] * 2 + jtwo / edge_length_list[edge_column_index] * 2 - epsnaught * edge_length_list[edge_column_index] / sum_over_edge_lengths - 2 * jtwo / sum_over_edge_lengths);
						gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, edge_column_index, -epsab + epsab * edge_length_list[edge_column_index] / sum_over_edge_lengths + 2 * jtwo / sum_over_edge_lengths);
						//if this edge is in contact with another differentiated cell, store the coupling to the neighboring concentration.
						if ((*pcp_column_cell_ptr).is_Cell_with_PCP()){
							gsl_matrix_set(U_matrix_by_cell, edge_row_index, 3 * number_of_edges + edge_column_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths + epsone);
							gsl_matrix_set(U_matrix_by_cell, edge_row_index, 2 * number_of_edges + edge_column_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths + epstwo);
							gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, 3 * number_of_edges + edge_column_index, -epstwo * edge_length_list[edge_column_index] / sum_over_edge_lengths + epstwo);
							gsl_matrix_set(U_matrix_by_cell, number_of_edges + edge_row_index, 2 * number_of_edges + edge_column_index, -epsone * edge_length_list[edge_column_index] / sum_over_edge_lengths + epsone);
						}
					}
				}
			}
			//Only select the non-zero region of the evolution_SUBVECTOR. There is more memory allocated than necessary.
			gsl_vector_view evolution_SUBVECTOR_view = gsl_vector_subvector(evolution_SUBVECTOR, 0, U_matrix_by_rows_size);
			gsl_vector* store_result = &evolution_SUBVECTOR_view.vector;
			//Only select the non-zero region of the concentration_vector_by_cell. There is more memory allocated than necessary.
			gsl_vector_view concentration_view_by_cell = gsl_vector_subvector(concentration_vector_by_cell, 0, U_matrix_by_columns_size);
			//Only select the non-zero region of the U_matrix_sub. There is more memory allocated than necessary.
			gsl_matrix_view U_matrix_sub = gsl_matrix_submatrix(U_matrix_by_cell, 0, 0, U_matrix_by_rows_size, U_matrix_by_columns_size);
			//Multiply UijCj.
			gsl_blas_dgemv(CblasNoTrans, 1.0, &U_matrix_sub.matrix, &concentration_view_by_cell.vector, 0.0, store_result);

			//Loop through edges within the cell, store the evolution vector due to U terms for both protein A and B.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				//U terms for protein variables of type A within this cell.
				gsl_vector_set(evolution_vectorU, start_index + edge_index, gsl_vector_get(evolution_SUBVECTOR, edge_index));
				//U terms for protein variables of type B within this cell.
				gsl_vector_set(evolution_vectorU, start_index + edge_index + type_dimension, gsl_vector_get(evolution_SUBVECTOR, number_of_edges + edge_index));
			}
			//Zero out all of the temporary matrix and vectors.
			gsl_matrix_set_zero(U_matrix_by_cell);
			gsl_vector_set_zero(concentration_vector_by_cell);
			gsl_vector_set_zero(evolution_SUBVECTOR);
		}
	}
	//Free the gsl matrix and vectors for this function.
	gsl_matrix_free(U_matrix_by_cell);
	gsl_vector_free(evolution_SUBVECTOR);
	gsl_vector_free(concentration_vector_by_cell);
}

//This is the function for initializing the V-matrix for all of the ode solving functions (for calculation of explicit terms).
//The matrix dimension of V is half that of U since V only couples protein concentrations of the same type.
//Consult LUBENSKY PLOS 2012 paper for definition of V-matrix (supp info)
//This is no longer used in the imex solvers because this requires the initialization of a matrix with many zeros.
void Initialize_V_Matrix(gsl_matrix * V, void* params, const vector<int> &start_index_cells){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;
	int start_index;

	//Then iterate through the cells again in order to build the V matrix for proteins of type A
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
			cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			start_index = start_index_cells.at(cell_index);
			//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			vector<double> edge_length_list;
			edge_length_list.resize(edge_list_for_this_cell.size());
			//Store a list of edges lengths for this cell.
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
			}
			//compute sum over edge lengths.
			sum_over_edge_lengths = 0.0;
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				sum_over_edge_lengths += edge_length_list[edge_index];
			}
			for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
				for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
					if (edge_column_index != edge_row_index){
						//Consult equation S1, S2, S10, S11 of PLOS 2012 paper to see where this element comes from. Note that to compute the the term from V's for a given concentration, one sums over product of columns of V with the corresponding state vector.
						gsl_matrix_set(V, edge_row_index + start_index, edge_column_index + start_index, -edge_length_list[edge_column_index] / sum_over_edge_lengths);
					}
					else {
						gsl_matrix_set(V, edge_row_index + start_index, edge_column_index + start_index, 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths);
					}
				}
			}
		}
	}
}

//This is the function for initializing the V-matrix based on new edge lengths (for the implicit term->next step).
//The matrix dimension of V is half that of U since V only couples protein concentrations of the same type.
//Consult LUBENSKY PLOS 2012 paper for definition of V-matrix (supp info)
//This is no longer used in the imex solvers because this requires the initialization of a matrix with many zeros.
void Initialize_V_Matrix_New(gsl_matrix * V, void* params, const vector<int> &start_index_cells){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;
	int start_index;
	//Then iterate through the cells again in order to build the V matrix for proteins of type A
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
			cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			start_index = start_index_cells.at(cell_index);
			//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			vector<double> edge_length_list;
			edge_length_list.resize(edge_list_for_this_cell.size());
			//Store a list of edges lengths for this cell.
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				//Calculate new edge length.
				edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength_New();
			}
			//compute sum over edge lengths.
			sum_over_edge_lengths = 0.0;
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				sum_over_edge_lengths += edge_length_list[edge_index];
			}
			for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
				for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
					if (edge_column_index != edge_row_index){
						//Consult equation S1, S2, S10, S11 of PLOS 2012 paper to see where this element comes from. Note that to compute the the term from V's for a given concentration, one sums over product of columns of V with the corresponding state vector.
						gsl_matrix_set(V, edge_row_index + start_index, edge_column_index + start_index, -edge_length_list[edge_column_index] / sum_over_edge_lengths);
					}
					else {
						gsl_matrix_set(V, edge_row_index + start_index, edge_column_index + start_index, 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths);
					}
				}
			}
		}
	}
}

//MAKE SURE TO ALLOCATE edge_list_for_this_cell.size() X edge_list_for_this_cell.size() for all Initialize_V_Matrix_Cell functions.
//This requirement applies for function that take as input gsl_matrix* or vector<vector<double>>

//Take as input a gsl_matrix* that points to a matrix of zeros with the specified size.
//Store the V submatrix for the given cell in that gsl_matrix.
void Initialize_V_Matrix_Cell(gsl_matrix * V, void* params, const vector<int> &start_index_cells, Cell* cell_ptr){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;

	PCPCell* pcpcit = dynamic_cast<PCPCell*>(cell_ptr);
	if ((*pcpcit).is_Cell_with_PCP()){
		//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
		cell_id = (*pcpcit).id;
		cell_index = cell_id2index.at(cell_id);
		//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
		vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
		vector<double> edge_length_list;
		edge_length_list.resize(edge_list_for_this_cell.size());
		//Store a list of edges lengths for this cell.
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
		}
		//compute sum over edge lengths.
		sum_over_edge_lengths = 0.0;
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			sum_over_edge_lengths += edge_length_list[edge_index];
		}
		for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
			for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
				if (edge_column_index != edge_row_index){
					gsl_matrix_set(V, edge_row_index, edge_column_index, -edge_length_list[edge_column_index] / sum_over_edge_lengths);
				}
				else {
					gsl_matrix_set(V, edge_row_index, edge_column_index, 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths);
				}
			}
		}
	}
}

//Take as input a vector<vector<double>> that points to a matrix of zeros with the specified size.
//Store the V submatrix for the given cell in that vector<vector<double>>.
void Initialize_V_Matrix_Cell(vector<vector<double>> &V, void* params, const vector<int> &start_index_cells, Cell* cell_ptr){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;

	PCPCell* pcpcit = dynamic_cast<PCPCell*>(cell_ptr);
	if ((*pcpcit).is_Cell_with_PCP()){
		//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
		cell_id = (*pcpcit).id;
		cell_index = cell_id2index.at(cell_id);
		//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
		vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
		vector<double> edge_length_list;
		edge_length_list.resize(edge_list_for_this_cell.size());
		//Store a list of edges lengths for this cell.
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
		}
		//compute sum over edge lengths.
		sum_over_edge_lengths = 0.0;
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			sum_over_edge_lengths += edge_length_list[edge_index];
		}
		for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
			for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
				if (edge_column_index != edge_row_index){
					V.at(edge_row_index).at(edge_column_index) =  -edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
				else {
					V.at(edge_row_index).at(edge_column_index) = 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
			}
		}
	}
}

//Store the V submatrix for the given cell in V_matrix owned by the cell corresponding to cell_ptr.
void Initialize_V_Matrix_Cell(void* params, const vector<int> &start_index_cells, Cell* cell_ptr){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;
	int number_of_edges;

	PCPCell* pcpcit = dynamic_cast<PCPCell*>(cell_ptr);
	if ((*pcpcit).is_Cell_with_PCP()){
		//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
		cell_id = (*pcpcit).id;
		cell_index = cell_id2index.at(cell_id);
		//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
		vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
		vector<double> edge_length_list;
		number_of_edges = edge_list_for_this_cell.size();
		edge_length_list.resize(number_of_edges);
		//Store V in a C++ matrix  (not absolutely necessary and will probably change).
		double init_val = 0.0;
		vector<vector<double>> vec(number_of_edges, vector<double>(number_of_edges, init_val));
		//Use the V matrix based on new edge lengths.
		vector<vector<double>> V_matrix = vec;
		//Store a list of edges lengths for this cell.
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
		}
		//compute sum over edge lengths.
		sum_over_edge_lengths = 0.0;
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			sum_over_edge_lengths += edge_length_list[edge_index];
		}
		for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
			for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
				if (edge_column_index != edge_row_index){
					V_matrix.at(edge_row_index).at(edge_column_index) = -edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
				else {
					V_matrix.at(edge_row_index).at(edge_column_index) = 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
			}
		}
		(*pcpcit).V_matrix = V_matrix;
	}
}

//Take as input a gsl_matrix* that points to a matrix of zeros with the specified size.
//Store the V submatrix AT THE NEXT STEP for the given cell in that gsl_matrix.
void Initialize_V_Matrix_Cell_New(gsl_matrix * V, void* params, const vector<int> &start_index_cells, Cell* cell_ptr){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;

	PCPCell* pcpcit = dynamic_cast<PCPCell*>(cell_ptr);
	if ((*pcpcit).is_Cell_with_PCP()){
		//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
		cell_id = (*pcpcit).id;
		cell_index = cell_id2index.at(cell_id);
		//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
		vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
		vector<double> edge_length_list;
		edge_length_list.resize(edge_list_for_this_cell.size());
		//Store a list of edges lengths for this cell.
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
		}
		//compute sum over edge lengths.
		sum_over_edge_lengths = 0.0;
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			sum_over_edge_lengths += edge_length_list[edge_index];
		}
		for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
			for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
				if (edge_column_index != edge_row_index){
					gsl_matrix_set(V, edge_row_index, edge_column_index, -edge_length_list[edge_column_index] / sum_over_edge_lengths);
				}
				else {
					gsl_matrix_set(V, edge_row_index, edge_column_index, 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths);
				}
			}
		}
	}
}

//Take as input a vector<vector<double>> that points to a matrix of zeros with the specified size.
//Store the V submatrix AT THE NEXT STEP for the given cell in that vector<vector<double>>.
void Initialize_V_Matrix_Cell_New(vector<vector<double>> &V, void* params, const vector<int> &start_index_cells, Cell* cell_ptr){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;

	PCPCell* pcpcit = dynamic_cast<PCPCell*>(cell_ptr);
	if ((*pcpcit).is_Cell_with_PCP()){
		//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
		cell_id = (*pcpcit).id;
		cell_index = cell_id2index.at(cell_id);
		//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
		vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
		vector<double> edge_length_list;
		edge_length_list.resize(edge_list_for_this_cell.size());
		//Store a list of edges lengths for this cell.
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength();
		}
		//compute sum over edge lengths.
		sum_over_edge_lengths = 0.0;
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			sum_over_edge_lengths += edge_length_list[edge_index];
		}
		for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
			for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
				if (edge_column_index != edge_row_index){
					V.at(edge_row_index).at(edge_column_index) = -edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
				else {
					V.at(edge_row_index).at(edge_column_index) = 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
			}
		}
	}
}

//Store the V submatrix AT THE NEXT STEP for the given cell in V_matrix owned by the cell corresponding to cell_ptr.
void Initialize_V_Matrix_Cell_New(void* params, const vector<int> &start_index_cells, Cell* cell_ptr){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	int type_dimension = dim_proteins / 2;
	//Maybe keep track of how many edges each cell has(in some sort of a vector).Then, just operate block diagonally.
	double sum_over_edge_lengths = 0.0;
	int cell_id;
	int cell_index;
	int number_of_edges;

	PCPCell* pcpcit = dynamic_cast<PCPCell*>(cell_ptr);
	if ((*pcpcit).is_Cell_with_PCP()){
		//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell.
		cell_id = (*pcpcit).id;
		cell_index = cell_id2index.at(cell_id);
		//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
		vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
		vector<double> edge_length_list;
		number_of_edges = edge_list_for_this_cell.size();
		edge_length_list.resize(number_of_edges);
		//Store V in a C++ matrix  (not absolutely necessary and will probably change).
		double init_val = 0.0;
		vector<vector<double>> vec(number_of_edges, vector<double>(number_of_edges, init_val));
		//Use the V matrix based on new edge lengths.
		vector<vector<double>> V_matrix = vec;
		//Store a list of edges lengths for this cell.
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			edge_length_list[edge_index] = (*edge_list_for_this_cell[edge_index]).GetLength_New();
		}
		//compute sum over edge lengths.
		sum_over_edge_lengths = 0.0;
		for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
			sum_over_edge_lengths += edge_length_list[edge_index];
		}
		for (int edge_row_index = 0; edge_row_index != edge_list_for_this_cell.size(); edge_row_index++){
			for (int edge_column_index = 0; edge_column_index != edge_list_for_this_cell.size(); edge_column_index++){
				if (edge_column_index != edge_row_index){
					V_matrix.at(edge_row_index).at(edge_column_index) = -edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
				else {
					V_matrix.at(edge_row_index).at(edge_column_index) = 1 - edge_length_list[edge_column_index] / sum_over_edge_lengths;
				}
			}
		}
		//Store the V submatrix in (*pcpcit).V_matrix
		(*pcpcit).V_matrix = V_matrix;
	}
}

//Given a c++ matrix (square with side length matrix_length), store those values in V_GSL.
void convert_CMatrix_to_GSL(vector<vector<double>> &V_C, gsl_matrix* V_GSL, int matrix_length){
	for (int row_index = 0; row_index != matrix_length; row_index++){
		for (int column_index = 0; column_index != matrix_length; column_index++){
			gsl_matrix_set(V_GSL, row_index, column_index, V_C.at(row_index).at(column_index));
		}
	}
}

//MAKE SURE TO ALLOCATE dim_proteins
//Define the evolution vector from V terms (using y_vector, as well as V_matrix owned by all of the cells).
void Define_Evolution_Vector_V(gsl_vector* evolution_vectorV, gsl_vector* y_vector, void* params, const vector<int> &start_index_cells){
	//store cell_id2index.
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	//store number of variables for protein of type A (which is the same number of variables of type B).
	int type_dimension = (*(combinedparameters*)params).dimension_protein / 2;
	int start_index;
	//Initialize V_matrix_by_cell and evolution_SUBVECTOR to be larger than necessary.
	int max_number_of_edges = 100;
	gsl_matrix* V_matrix_by_cell = gsl_matrix_calloc(max_number_of_edges, max_number_of_edges);
	gsl_vector* evolution_SUBVECTOR = gsl_vector_calloc(max_number_of_edges);

	//Loop through all of the cells in the packing.
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		//Cast cell pointer.
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		//Check if cell has PCP.
		if ((*pcpcit).is_Cell_with_PCP()){
			//Store ordered edges.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			//Store number of edges.
			int number_of_edges = edge_list_for_this_cell.size();
			//Store index for the given cell given the cell id.
			int cell_index = cell_id2index.at((*pcpcit).id);
			//Now find where this cell would live in the y_vector.
			start_index = start_index_cells.at(cell_index);
			//Store vector views of type A protein and type B protein for this cell.
			gsl_vector_view yvectorA_view = gsl_vector_subvector(y_vector, start_index, number_of_edges);
			gsl_vector_view yvectorB_view = gsl_vector_subvector(y_vector, type_dimension + start_index, number_of_edges);
			//Store these vector views in a vector.
			gsl_vector* yvectorA = &yvectorA_view.vector;
			gsl_vector* yvectorB = &yvectorB_view.vector;
			//Take a piece of the too large evolution_SUBVECTOR above.
			gsl_vector_view evolution_SUBVECTOR_view = gsl_vector_subvector(evolution_SUBVECTOR, 0, number_of_edges);
			//Use this to store the result of Vlog(c)
			gsl_vector* store_result = &evolution_SUBVECTOR_view.vector;
			//Initialize the V_matrix for this cell.
			Initialize_V_Matrix_Cell(params, start_index_cells, *cit);
			//Convert the (*pcpcit).Vmatrix to a gsl matrix.
			convert_CMatrix_to_GSL((*pcpcit).V_matrix, V_matrix_by_cell, number_of_edges);
			//Now, now take the non-zero part.
			gsl_matrix_view V_matrix_sub = gsl_matrix_submatrix(V_matrix_by_cell, 0, 0, number_of_edges, number_of_edges);
			//Do Vijlog(cj) for this cell (type A)
			gsl_blas_dgemv(CblasNoTrans, 1.0, &V_matrix_sub.matrix, yvectorA, 0.0, store_result);
			//Store Vijlog(cj) for this cell (type A)
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(evolution_vectorV, start_index + edge_index, gsl_vector_get(evolution_SUBVECTOR, edge_index));
			}
			//Do Vijlog(cj) for this cell (type B)
			gsl_blas_dgemv(CblasNoTrans, 1.0, &V_matrix_sub.matrix, yvectorB, 0.0, store_result);
			//Store Vijlog(cj) for this cell (type B)
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(evolution_vectorV, type_dimension + start_index + edge_index, gsl_vector_get(evolution_SUBVECTOR, edge_index));	
			}
			//Zero out the temp matrix and the temp vector.
			gsl_matrix_set_zero(V_matrix_by_cell);
			gsl_vector_set_zero(evolution_SUBVECTOR);
		}
	}
	//Free temp matrix and temp vector memory.
	gsl_vector_free(evolution_SUBVECTOR);
	gsl_matrix_free(V_matrix_by_cell);
}

//MAKE SURE TO ALLOCATE dim_proteins
//Define the evolution vector from max terms. MAX CONCENTRATION
void Define_Evolution_Vector_Max_Conc(gsl_vector* evolution_vector_max_terms, gsl_vector* y_vector, gsl_vector* length_vector, void* params, const vector<int> &start_index_cells){
	//store cell_id2index.
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	//store number of variables for protein of type A (which is the same number of variables of type B).
	int type_dimension = (*(combinedparameters*)params).dimension_protein / 2;
	int start_index;
	//Initialize V_matrix_by_cell and evolution_SUBVECTOR to be larger than necessary.
	int max_number_of_edges = 100;
	gsl_matrix* V_matrix_by_cell = gsl_matrix_calloc(max_number_of_edges, max_number_of_edges);
	gsl_vector* evolution_SUBVECTOR = gsl_vector_calloc(max_number_of_edges);
	gsl_vector* log_cmax_minus_concentrationA = gsl_vector_calloc(max_number_of_edges);
	gsl_vector* log_cmax_minus_concentrationB = gsl_vector_calloc(max_number_of_edges);

	//Loop through all of the cells in the packing.
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		//Cast cell pointer.
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		//Check if cell has PCP.
		if ((*pcpcit).is_Cell_with_PCP()){
			//Store ordered edges.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			//Store number of edges.
			int number_of_edges = edge_list_for_this_cell.size();
			//Store index for the given cell given the cell id.
			int cell_index = cell_id2index.at((*pcpcit).id);
			//Now find where this cell would live in the y_vector.
			start_index = start_index_cells.at(cell_index);
			//Store vector views of type A protein and type B protein for this cell.
			gsl_vector_view yvectorA_view = gsl_vector_subvector(y_vector, start_index, number_of_edges);
			gsl_vector_view yvectorB_view = gsl_vector_subvector(y_vector, type_dimension + start_index, number_of_edges);
			//Store these vector views in a vector.
			gsl_vector* yvectorA = &yvectorA_view.vector;
			gsl_vector* yvectorB = &yvectorB_view.vector;

			//Store the vector view of length vector.
			gsl_vector_view length_vect_view = gsl_vector_subvector(length_vector, start_index, number_of_edges);
			//Store this view into a vector.
			gsl_vector* length_vect = &length_vect_view.vector;

			//Take a piece of the too large evolution_SUBVECTOR above.
			gsl_vector_view evolution_SUBVECTOR_view = gsl_vector_subvector(evolution_SUBVECTOR, 0, number_of_edges);
			gsl_vector* store_result = &evolution_SUBVECTOR_view.vector;

			gsl_vector_view log_cmax_minus_concentrationA_view = gsl_vector_subvector(log_cmax_minus_concentrationA, 0, number_of_edges);
			gsl_vector_view log_cmax_minus_concentrationB_view = gsl_vector_subvector(log_cmax_minus_concentrationB, 0, number_of_edges);

			gsl_vector* log_cmax_minus_concentrationA_vect = &log_cmax_minus_concentrationA_view.vector;
			gsl_vector* log_cmax_minus_concentrationB_vect = &log_cmax_minus_concentrationB_view.vector;
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(log_cmax_minus_concentrationA_vect, edge_index, log(max_concentration - exp(-gsl_vector_get(yvectorA, edge_index)) / gsl_vector_get(length_vect, edge_index)));
				gsl_vector_set(log_cmax_minus_concentrationB_vect, edge_index, log(max_concentration - exp(-gsl_vector_get(yvectorB, edge_index)) / gsl_vector_get(length_vect, edge_index)));
			}

			//Initialize the V_matrix for this cell.
			Initialize_V_Matrix_Cell(params, start_index_cells, *cit);

			//Convert the (*pcpcit).Vmatrix to a gsl matrix.
			convert_CMatrix_to_GSL((*pcpcit).V_matrix, V_matrix_by_cell, number_of_edges);

			//Now, now take the non-zero part.
			gsl_matrix_view V_matrix_sub = gsl_matrix_submatrix(V_matrix_by_cell, 0, 0, number_of_edges, number_of_edges);

			gsl_blas_dgemv(CblasNoTrans, 1.0, &V_matrix_sub.matrix, log_cmax_minus_concentrationA_vect, 0.0, store_result);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(evolution_vector_max_terms, start_index + edge_index, gsl_vector_get(evolution_SUBVECTOR, edge_index));
			}
			

			gsl_blas_dgemv(CblasNoTrans, 1.0, &V_matrix_sub.matrix, log_cmax_minus_concentrationB_vect, 0.0, store_result);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(evolution_vector_max_terms, type_dimension + start_index + edge_index, gsl_vector_get(evolution_SUBVECTOR, edge_index));
			}

			gsl_matrix_set_zero(V_matrix_by_cell);
			gsl_vector_set_zero(evolution_SUBVECTOR);
			gsl_vector_set_zero(log_cmax_minus_concentrationA);
			gsl_vector_set_zero(log_cmax_minus_concentrationB);
		}
	}
	//Free temp matrix and temp vector memory.
	gsl_vector_free(evolution_SUBVECTOR);
	gsl_matrix_free(V_matrix_by_cell);
	gsl_vector_free(log_cmax_minus_concentrationA);
	gsl_vector_free(log_cmax_minus_concentrationB);
}

//MAKE SURE TO ALLOCATE dim_proteins
// WARNING!: This requires that you use Initialize_V_Matrix_Cell(void* params, const vector<int> &start_index_cells, Cell* cell_ptr) to store V_matrix for all cells.
void Define_Vloglength(gsl_vector* Vlogj, gsl_vector* loglengthvector, void* params, const vector<int> &start_index_cells){
	//store cell_id2index.
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;
	//store number of variables for protein of type A (which is the same number of variables of type B).
	int type_dimension = (*(combinedparameters*)params).dimension_protein / 2;
	int start_index;
	//Initialize V_matrix_by_cell and evolution_SUBVECTOR to be larger than necessary.
	int max_number_of_edges = 100;
	gsl_matrix* V_matrix_by_cell = gsl_matrix_calloc(max_number_of_edges, max_number_of_edges);
	gsl_vector* Vloglength_vect = gsl_vector_calloc(max_number_of_edges);

	//Loop through all of the cells in the packing.
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		//Cast cell pointer.
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		//Check if cell has PCP.
		if ((*pcpcit).is_Cell_with_PCP()){
			//Store ordered edges.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			//Store number of edges.
			int number_of_edges = edge_list_for_this_cell.size();
			//Store index for the given cell given the cell id.
			int cell_index = cell_id2index.at((*pcpcit).id);
			//Now find where this cell would live in the y_vector.
			start_index = start_index_cells.at(cell_index);
			//Store vector views of lengths for this cell.
			gsl_vector_view loglengthvector_view = gsl_vector_subvector(loglengthvector, start_index, number_of_edges);
			gsl_vector* temp_loglengthvector = &loglengthvector_view.vector;
			gsl_vector_view Vlogj_SUBVECTOR_view = gsl_vector_subvector(Vloglength_vect, 0, number_of_edges);
			gsl_vector* store_result = &Vlogj_SUBVECTOR_view.vector;
			//Use the below code to initialize V matrix for this cell if not already computed.
			//Initialize_V_Matrix_Cell(V_matrix_by_cell, params, start_index_cells, *cit);

			//Convert the (*pcpcit).Vmatrix to a gsl matrix.
			convert_CMatrix_to_GSL((*pcpcit).V_matrix, V_matrix_by_cell, number_of_edges);
			//Now, now take the non-zero part.
			gsl_matrix_view V_matrix_sub = gsl_matrix_submatrix(V_matrix_by_cell, 0, 0, number_of_edges, number_of_edges);
			//Do Vijlog(length) for this cell
			gsl_blas_dgemv(CblasNoTrans, 1.0, &V_matrix_sub.matrix, temp_loglengthvector, 0.0, store_result);

			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(Vlogj, start_index + edge_index, gsl_vector_get(Vloglength_vect, edge_index));
			}

			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(Vlogj, type_dimension + start_index + edge_index, gsl_vector_get(Vloglength_vect, edge_index));
			}
			gsl_matrix_set_zero(V_matrix_by_cell);
			gsl_vector_set_zero(Vloglength_vect);
		}
	}
	gsl_vector_free(Vloglength_vect);
	gsl_matrix_free(V_matrix_by_cell);
}

//Define the vector of concentrations for multiplication with the U-matrix.
void Define_Concentration_Vector(gsl_vector * concentration_vector, void* params, const vector<double> &y, gsl_vector * longlengthvector){
	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;

	//For each index in the y vector and length vector, exponentiate the log variable and divide by the corresponding length.
	for (int ind = 0; ind != dim_proteins; ind++){
		gsl_vector_set(concentration_vector, ind, exp(-y.at(ind)) / gsl_vector_get(longlengthvector, ind));
	}
}

//Store evolve info (past mechanical variable states and mechanical forces) -> fully explicit implementation.
void store_evolve_info(void* params, double step_size){
	//Store the forces and mechanical state at the given time point (for the vertices).
	int i = 0;
	for (list<Vertex*>::iterator vit = tissue.vertices.begin(); vit != tissue.vertices.end(); vit++){
		tissue.f_mech_step_before.at(2).at(i) = (*(*vit)).SavedInfoForEvolve_TotalForce;
		tissue.past_mech_state.at(2).at(i) = (*(*vit)).r;
		i++;
	}
	//Store the forces and mechanical state at the given time point (for the box size).
	tissue.f_mech_step_before.at(2).at(i) = Vect(tissue.ForceLX, tissue.ForceLY);
	tissue.past_mech_state.at(2).at(i) = Vect(LX, LY);
	//This is also stored in the corresponding imex solver (first-, second-, and third- order), but if mech. dynamics are evolved without prot dynamics, it needs to be stored here.
	tissue.past_step_sizes.at(2) = step_size;
}

//Helper function for the first order mechanical evolution (euler step)
void helper_imex_first_order_mech(void* params, double step_size, ofstream* errorfile){
	//Store evolved vertex positions and box size.
	tissue.Integrate_Euiler(step_size, errorfile);

	//Save mechanical state and mechanical forces for the state before the euler step.
	store_evolve_info(params, step_size);
}

//Helper function for the second order mechanical evolution (in the spirit of variable step imex)
void helper_imex_second_order_mech(void* params, double step_size, ofstream* errorfile){
	//Store evolved vertex positions and box size.
	tissue.Integrate_IMEX_Second_Order(step_size, params, errorfile);

	//Save mechanical state and mechanical forces for the state before the second-order step.
	store_evolve_info(params, step_size);
}

//Helper function for the third order mechanical evolution (in the spirit of variable step imex)
void helper_imex_third_order_mech(void* params, double step_size, ofstream* errorfile){
	//Store evolved vertex positions and box size.
	tissue.Integrate_IMEX_Third_Order(step_size, params, errorfile);

	//Save mechanical state and mechanical forces for the state before the third-order step.
	store_evolve_info(params, step_size);
}

//This is the full function to numerically solve the PCP dynamics to first order.
vector<double> imex_solver_PCP_first_order(const vector<double> &y, double step_size, void* params, ofstream* errorfile){
	//Evolve the mechanical variables (based on euler explicit step). Also, store information from this step.
	if (evolve_mech_positions == 1){
		helper_imex_first_order_mech(params, step_size, errorfile);
	}

	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	//U is composed of four submatrices (A-A, A-B, B-A, B-B). Define this variable to store the number of protein concentrations variables for a given type (A or B).
	int type_dimension = dim_proteins / 2;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;

	//for storing information from this step.
	vector<double> g_this_step_vector;
	g_this_step_vector.resize(dim_proteins);
	vector<double> f_this_step_vector;
	f_this_step_vector.resize(dim_proteins);
	vector<double> y_this_step_vector;
	y_this_step_vector.resize(dim_proteins);

	//Initialize loop variables.
	int cell_id;
	int cell_index;
	int start_index;

	//This vector is important. It allows you to figure out what row or column index refers to a given cell's edges in the loop for defining U below.
	//Given some cell index, this vector will give you where the entries for this cell begin.
	vector<int> start_index_cells;
	vector<PCPEdge*> edge_list_for_this_cell;
	//Resize this vector such that it fits all cells in the tissue.
	start_index_cells.resize(tissue.cells.size());
	//Initialize start_index to zero.
	start_index = 0;
	//Use ind to refer to the current cell index.
	cell_index = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){

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

	//Store the lengths for all edges.
	gsl_vector * lengthvect = gsl_vector_alloc(type_dimension);
	gsl_vector * lengthvect_new = gsl_vector_alloc(type_dimension);
	int ind = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			int cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				gsl_vector_set(lengthvect, ind, (*edge_list_for_this_cell[edge_index]).GetLength());
				gsl_vector_set(lengthvect_new, ind, (*edge_list_for_this_cell[edge_index]).GetLength_New());
				ind++;
			}
		}
	}
	//This may seem redundant, but it is useful. Create a vector which is two copies of lengthvectforjacob stacked on each other.
	gsl_vector * longlengthvector = gsl_vector_alloc(dim_proteins);
	gsl_vector * longlengthvector_new = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(longlengthvector, ind, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector, ind + type_dimension, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector_new, ind, gsl_vector_get(lengthvect_new, ind));
		gsl_vector_set(longlengthvector_new, ind + type_dimension, gsl_vector_get(lengthvect_new, ind));
	}

	//Define the dimension of the matrix U.
	int U_dimension = dim_proteins;
	//Store the concentrations of edges that correspond to row/column indices of U.
	gsl_vector * concentration_vector = gsl_vector_alloc(U_dimension);
	//Store the vector that defines Uijcj for a given protein index i.
	gsl_vector * evolution_vectorU = gsl_vector_alloc(U_dimension);

	//Define concentration vector for multiplying the U matrix.
	Define_Concentration_Vector(concentration_vector, params, y, longlengthvector);
	//Define the vector for UijCj couplings.
	Define_Evolution_Vector_U(evolution_vectorU, concentration_vector, params, start_index_cells);

	//Define a vector for storing length*concentration=protein amount for type A and B (A, then B).
	gsl_vector * lengthtimesconcentration = gsl_vector_alloc(dim_proteins);
	//Copy longlengthvector into lengthtimesconcentration.
	gsl_vector_memcpy(lengthtimesconcentration, longlengthvector);
	//This multiplies concentration by length elementwise and stores in lengthtimesconcentration
	gsl_vector_mul(lengthtimesconcentration, concentration_vector);

	//Store the log of the length of each edge.
	gsl_vector * loglength = gsl_vector_alloc(type_dimension);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(loglength, ind, log(gsl_vector_get(lengthvect, ind)));
	}

	/*********This is where I begin building the right hand side of the third-order scheme*********/
	//Please note that GSL vector conventions forces me to directly multiply the vector with necessary factor.
	//When reading this code, pay attention to the gsl_vector_mul, gsl_vector_scale commands.

	//g contains the implicitly handled term. f handles the explicitly handled term.

	//Will store the current -log variable values.
	gsl_vector * y_this_step = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != dim_proteins; ind++){
		gsl_vector_set(y_this_step, ind, y.at(ind));
	}

	//Initialize g_this_step which we store -te*V*X_vector for both types of proteins.
	gsl_vector * g_this_step = gsl_vector_alloc(dim_proteins);
	gsl_vector * evolution_vectorV = gsl_vector_alloc(dim_proteins);
	Define_Evolution_Vector_V(evolution_vectorV, y_this_step, params, start_index_cells);
	//Scale by -Te.
	gsl_vector_scale(evolution_vectorV, -1. * te);

	//Add max concentration term.
	gsl_vector * max_conc_term = gsl_vector_alloc(dim_proteins);
	if (max_concentration != 0.0){
		Define_Evolution_Vector_Max_Conc(max_conc_term, y_this_step, lengthvect, params, start_index_cells);
		//Scale by -Te.
		gsl_vector_scale(max_conc_term, -1. * te);
	}

	for (int ind = 0; ind != dim_proteins; ind++){
		if (max_concentration != 0.0){
			gsl_vector_set(g_this_step, ind, gsl_vector_get(evolution_vectorV, ind) + gsl_vector_get(max_conc_term, ind));
		}
		else{
			gsl_vector_set(g_this_step, ind, gsl_vector_get(evolution_vectorV, ind));
		}
	}
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(g_this_step, longlengthvector);

	//Create a vector to store two identical copies of Vlogj_vector (one for type A and one for type B).
	gsl_vector * V_logj_full = gsl_vector_alloc(dim_proteins);
	Define_Vloglength(V_logj_full, loglength, params, start_index_cells);
	gsl_vector_scale(V_logj_full, te);

	//Store the explicit terms.
	gsl_vector * f_this_step = gsl_vector_alloc(dim_proteins);
	//Subtract the Vlogj for both type A and B
	gsl_vector_sub(f_this_step, V_logj_full);
	//Subtract the contribution for the U terms
	gsl_vector_sub(f_this_step, evolution_vectorU);
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(f_this_step, longlengthvector);

	//Now store g_this_step, f_this_step, y_this_step into c++ vector so that they can be used in future steps (third-order solver).
	for (int edge_index = 0; edge_index != dim_proteins; edge_index++){
		g_this_step_vector.at(edge_index) = gsl_vector_get(g_this_step, edge_index);
		f_this_step_vector.at(edge_index) = gsl_vector_get(f_this_step, edge_index);
		y_this_step_vector.at(edge_index) = gsl_vector_get(y_this_step, edge_index);
	}
	
	//Create a vector to store the right hand side of equation. (CONSULT PSEUDOCODE).
	gsl_vector * rhs_evolution = gsl_vector_alloc(dim_proteins);

	//Now scale y_this_step (cannot do before since it is used for exponential corrections)
	gsl_vector_scale(y_this_step, 1./ step_size * alphaP);
	//note lengthtimesconcentration is nothing but exp(-y)
	gsl_vector_mul(y_this_step, lengthtimesconcentration);

	//Add y_this_step to f_this_step for right hand side (this is for a forward-backward euler combo).
	gsl_vector_add(rhs_evolution, y_this_step);
	gsl_vector_add(rhs_evolution, f_this_step);


	//Create a vector for storing the updated protein variable values.
	vector<double> updated_protein_variables;
	updated_protein_variables.resize(dim_proteins);


	//Note: one cannot use a cell iterator in a parallelized for loop, must build pointer list.
	vector<Cell*> CellPointerList;
	CellPointerList.reserve(type_dimension);
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		CellPointerList.push_back(&(*(*cit)));
	}

	//Parallelize the for loop (solve PCP dynamics on a cell-by-cell basis).
#pragma omp parallel for schedule(guided)															//BEGIN parallel for
	//not sure if parallelization is worth it here (probably not).
	for (int i = 0; i < CellPointerList.size(); i++){
		Cell* cptr = CellPointerList[i];
		PCPCell* pcpcptr = dynamic_cast<PCPCell*>(cptr);
		//If the cell is differentiated,
		if ((*pcpcptr).is_Cell_with_PCP()){
			//This struct will store parameters for numerical solving.
			struct rparams p;
			p.delta_t = step_size;
			int iter = 0;
			int status;
			int cell_id = (*cptr).id;
			int cell_index = cell_id2index.at(cell_id);
			int start_index = start_index_cells[cell_index];
			int completion_code = 27;
			//Store the list of edges for a given cell.
			//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcptr).FindOrderedEdges();
			//Store the number of edges within this cell.
			int number_of_edges = edge_list_for_this_cell.size();
			vector<double> edge_lengths;
			edge_lengths.resize(number_of_edges);
			//define a length vector
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				//Use the new edge lengths.
				if (evolve_mech_positions == 1){
					edge_lengths.at(edge_index) = gsl_vector_get(lengthvect_new, edge_index + start_index);
				}
				else{
					edge_lengths.at(edge_index) = gsl_vector_get(lengthvect, edge_index + start_index);
				}
			}

			//Now assign this to rparams.
			p.length_edges = edge_lengths;

			gsl_vector_view rhs_evolution_viewA_this_cell = gsl_vector_subvector(rhs_evolution, start_index, number_of_edges);
			gsl_vector_view rhs_evolution_viewB_this_cell = gsl_vector_subvector(rhs_evolution, start_index + type_dimension, number_of_edges);
			////Then, store this view into an actual gsl vector.
			gsl_vector * rhs_evolution_type_A = &rhs_evolution_viewA_this_cell.vector;
			gsl_vector * rhs_evolution_type_B = &rhs_evolution_viewB_this_cell.vector;

			//Store the current variable log variables for protein type A.
			vector<double> y_inA_this_cell;
			y_inA_this_cell.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				y_inA_this_cell.at(edge_index) = y.at(edge_index + start_index);
			}

			//Store the current variable log variables for protein type B.
			vector<double> y_inB_this_cell;
			y_inB_this_cell.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				y_inB_this_cell.at(edge_index) = y.at(edge_index + start_index + type_dimension);
			}

			//Store the current rhs of the evolution equation for protein type A.
			vector<double> rhs_evolution_A;
			rhs_evolution_A.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				rhs_evolution_A.at(edge_index) = gsl_vector_get(rhs_evolution_type_A, edge_index);
			}

			//Store the current rhs of the evolution equation for protein type B.
			vector<double> rhs_evolution_B;
			rhs_evolution_B.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				rhs_evolution_B.at(edge_index) = gsl_vector_get(rhs_evolution_type_B, edge_index);
			}

			if (evolve_mech_positions == 1){
				Initialize_V_Matrix_Cell_New(params, start_index_cells, cptr);
			}
			else{
				Initialize_V_Matrix_Cell(params, start_index_cells, cptr);
			}

			//Now, assign the V_matrix to the parameters for the numerical solver.
			p.V_matrix = (*pcpcptr).V_matrix;
			//first solve for protein A 
			p.y_in = y_inA_this_cell;
			p.rhs_evaluation = rhs_evolution_A;

			//Initialize the solver and the function.
			gsl_multiroot_fdfsolver *s_A;
			gsl_multiroot_function_fdf f_A = { &numerical_solve_first_order, &numerical_solve_first_order_df, &numerical_solve_first_order_fdf, number_of_edges, &p };

			//Initialize the vector for storing initial values.
			gsl_vector * x_A = gsl_vector_alloc(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(x_A, edge_index, y_inA_this_cell.at(edge_index));
			}

			//Use hybrid solver.
			nonlinear_solver_PCP_first_order = gsl_multiroot_fdfsolver_hybridj;
			//Specify type of solver and dimension of the system to be solved.
			s_A = gsl_multiroot_fdfsolver_alloc(nonlinear_solver_PCP_first_order, number_of_edges);
			gsl_multiroot_fdfsolver_set(s_A, &f_A, x_A);

			double solver_tolerance = 1e-7;
			int max_iterations = 1000;
			try{
				//Iterate the solver until it reaches the appropriate tolerance or exceeds max iterations.
				do {

					status = gsl_multiroot_fdfsolver_iterate(s_A);

					if (status)
						break;

					//define a variable for this residual test.
					status = gsl_multiroot_test_residual(s_A->f, solver_tolerance);
					iter++;
				}
				//define variable for iter threshold
				while (status == GSL_CONTINUE || iter < max_iterations);
				throw status;
			}
			//Catch error from the root finding loop.
			catch (int e){
				if (e != completion_code){
					printf("error: %s\n", gsl_strerror(e));
				}
			}


			//Store the evolved protein variables.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				updated_protein_variables.at(start_index + edge_index) = gsl_vector_get(s_A->x, edge_index);
			}

			//IMPORTANT NOTE: YOU MUST DEALLOCATE THE VARIABLES FROM THE TYPE A SOLVER. (CANNOT REUSE FOR SECOND->WILL CAUSE MEMORY LEAK).
			gsl_vector_free(x_A);
			gsl_multiroot_fdfsolver_free(s_A);

			//Declare solver variables for type B.
			gsl_multiroot_fdfsolver *s_B;
			gsl_multiroot_function_fdf f_B = { &numerical_solve_first_order, &numerical_solve_first_order_df, &numerical_solve_first_order_fdf, number_of_edges, &p };

			//Declare vector for initial condition.
			gsl_vector * x_B = gsl_vector_alloc(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(x_B, edge_index, y_inB_this_cell.at(edge_index));
			}

			//Store parameters for solver.
			p.y_in = y_inB_this_cell;
			p.rhs_evaluation = rhs_evolution_B;

			//Declare hybrid solver type. (CONSULT GSL DOCUMENTATION).
			nonlinear_solver_PCP_first_order = gsl_multiroot_fdfsolver_hybridj;
			//Allocate the solver and set the relevant properties.
			s_B = gsl_multiroot_fdfsolver_alloc(nonlinear_solver_PCP_first_order, number_of_edges);
			gsl_multiroot_fdfsolver_set(s_B, &f_B, x_B);

			iter = 0;
			try{
				//Iterate the solver until it reaches the appropriate tolerance or exceeds max iterations.
				do {
					status = gsl_multiroot_fdfsolver_iterate(s_B);

					if (status)
						break;

					status = gsl_multiroot_test_residual(s_B->f, solver_tolerance);
					iter++;
				}
				while (status == GSL_CONTINUE || iter < max_iterations);
				throw status;
			}
			//Catch error from the root finding loop.
			catch (int e){
				if (e != completion_code){
					printf("error: %s\n", gsl_strerror(e));
				}
			}

			//Store the updated protein variable values.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				updated_protein_variables.at(start_index + edge_index + type_dimension) = gsl_vector_get(s_B->x, edge_index);
			}

			//IMPORTANT NOTE: YOU MUST DEALLOCATE THE VARIABLES FROM THE TYPE B SOLVER.
			gsl_vector_free(x_B);
			gsl_multiroot_fdfsolver_free(s_B);
		}
	}


	//Free the memory associated with all of the matrices and vectors used in gsl. Don't allow for memory leaks.
	gsl_vector_free(rhs_evolution);
	gsl_vector_free(concentration_vector);
	gsl_vector_free(evolution_vectorU);
	gsl_vector_free(lengthvect);
	gsl_vector_free(lengthvect_new);
	gsl_vector_free(lengthtimesconcentration);
	gsl_vector_free(longlengthvector);
	gsl_vector_free(longlengthvector_new);
	gsl_vector_free(loglength);
	gsl_vector_free(y_this_step);
	gsl_vector_free(g_this_step);
	gsl_vector_free(f_this_step);
	gsl_vector_free(V_logj_full);
	gsl_vector_free(evolution_vectorV);
	gsl_vector_free(max_conc_term);

	//Store information from this step in the last column of the vectors that store information from past steps.
	tissue.f_step_before[2] = f_this_step_vector;
	tissue.y_step_before[2] = y_this_step_vector;
	tissue.g_step_before[2] = g_this_step_vector;
	tissue.past_step_sizes[2] = step_size;

	return updated_protein_variables;
}

//This is the full function to numerically solve the PCP dynamics to third order.
vector<double> imex_solver_PCP_second_order(const vector<double> &y, double step_size, void* params, ofstream* errorfile){
	//Define coefficients for solving method.
	double gamma = 1. / 2.;
	double c = 1. / 8.;

	//Use second order step based on imex scheme (where the implicit variables are strictly zero).
	//Also, store information from this step.
	if (evolve_mech_positions == 1){
		helper_imex_second_order_mech(params, step_size, errorfile);
	}

	//info from previous time steps for higher order
	vector<vector<double>> g_step_before = tissue.g_step_before;
	vector<vector<double>> f_step_before = tissue.f_step_before;
	vector<vector<double>> y_step_before = tissue.y_step_before;

	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	//U is composed of four submatrices (A-A, A-B, B-A, B-B). Define this variable to store the number of protein concentrations variables for a given type (A or B).
	int type_dimension = dim_proteins / 2;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;

	//for storing information from this step.
	vector<double> g_this_step_vector;
	g_this_step_vector.resize(dim_proteins);
	vector<double> f_this_step_vector;
	f_this_step_vector.resize(dim_proteins);
	vector<double> y_this_step_vector;
	y_this_step_vector.resize(dim_proteins);

	//Declare loop variables.
	int cell_id;
	int cell_index;
	int start_index;

	//This vector is important. It allows you to figure out what row or column index refers to a given cell's edges in the loop for defining U below.
	//Given some cell index, this vector will give you where the entries for this cell begin.
	vector<int> start_index_cells;
	vector<PCPEdge*> edge_list_for_this_cell;
	//Resize this vector such that it fits all cells in the tissue.
	start_index_cells.resize(tissue.cells.size());
	//Initialize start_index to zero.
	start_index = 0;
	//Use ind to refer to the current cell index.
	cell_index = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){

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

	//Store the lengths for all edges.
	gsl_vector * lengthvect = gsl_vector_alloc(type_dimension);
	gsl_vector * lengthvect_new = gsl_vector_alloc(type_dimension);
	int ind = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			int cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				gsl_vector_set(lengthvect, ind, (*edge_list_for_this_cell[edge_index]).GetLength());
				gsl_vector_set(lengthvect_new, ind, (*edge_list_for_this_cell[edge_index]).GetLength_New());
				ind++;
			}
		}
	}
	//This may seem redundant, but it is useful. Create a vector which is two copies of lengthvectforjacob stacked on each other.
	gsl_vector * longlengthvector = gsl_vector_alloc(dim_proteins);
	gsl_vector * longlengthvector_new = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(longlengthvector, ind, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector, ind + type_dimension, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector_new, ind, gsl_vector_get(lengthvect_new, ind));
		gsl_vector_set(longlengthvector_new, ind + type_dimension, gsl_vector_get(lengthvect_new, ind));
	}

	//Define the dimension of the matrix U.
	int U_dimension = dim_proteins;
	//Store the concentrations of edges that correspond to row/column indices of U.
	gsl_vector * concentration_vector = gsl_vector_alloc(U_dimension);
	//Store the vector that defines Uijcj for a given protein index i.
	gsl_vector * evolution_vectorU = gsl_vector_alloc(U_dimension);

	//Define the vector of concentrations for multiplying the U matrix.
	Define_Concentration_Vector(concentration_vector, params, y, longlengthvector);
	//Define the evolution terms from U couplings.
	Define_Evolution_Vector_U(evolution_vectorU, concentration_vector, params, start_index_cells);

	//Define a vector for storing length*concentration=protein amount for type A and B (A, then B).
	gsl_vector * lengthtimesconcentration = gsl_vector_alloc(dim_proteins);
	//Copy longlengthvector into lengthtimesconcentration.
	gsl_vector_memcpy(lengthtimesconcentration, longlengthvector);
	//This multiplies concentration by length elementwise and stores in lengthtimesconcentration
	gsl_vector_mul(lengthtimesconcentration, concentration_vector);

	//Store the log of the length of each edge.
	gsl_vector * loglength = gsl_vector_alloc(type_dimension);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(loglength, ind, log(gsl_vector_get(lengthvect, ind)));
	}


	//Define ratios of step sizes.
	double omega_n_plus_1 = step_size / tissue.past_step_sizes[1];

	//Define coefficients for third order equations.
	double alpha0n = (2. * gamma - 1.) * pow(omega_n_plus_1, 2.) / (1. + omega_n_plus_1);
	double alpha1n = (1. - 2. * gamma) * omega_n_plus_1 - 1.;
	double alpha2n = (1. + 2. * gamma * omega_n_plus_1) / (1. + omega_n_plus_1);
	double beta0n = -gamma * omega_n_plus_1;
	double beta1n = 1. + gamma * omega_n_plus_1;
	double gamma0n = c / 2.;
	double gamma1n = 1. - gamma - (1. + 1. / omega_n_plus_1) * c / 2.;
	double gamma2n = gamma + c / 2. / omega_n_plus_1;


	/*********This is where I begin building the right hand side of the third-order scheme*********/
	//Please note that GSL vector conventions forces me to directly multiply the vector with necessary factor.
	//When reading this code, pay attention to the gsl_vector_mul, gsl_vector_scale commands.

	//g contains the implicitly handled term. f handles the explicitly handled term.

	//Will store the current -log variable values.
	gsl_vector * y_this_step = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != dim_proteins; ind++){
		gsl_vector_set(y_this_step, ind, y.at(ind));
	}

	//Initialize g_this_step which we store -te*V*X_vector for both types of proteins.
	gsl_vector * g_this_step = gsl_vector_alloc(dim_proteins);
	gsl_vector * evolution_vectorV = gsl_vector_alloc(dim_proteins);
	Define_Evolution_Vector_V(evolution_vectorV, y_this_step, params, start_index_cells);
	//Scale by -Te.
	gsl_vector_scale(evolution_vectorV, -1. * te);

	//Add max concentration term.
	gsl_vector * max_conc_term = gsl_vector_alloc(dim_proteins);
	if (max_concentration != 0.0){
		Define_Evolution_Vector_Max_Conc(max_conc_term, y_this_step, lengthvect, params, start_index_cells);
		//Scale by -Te.
		gsl_vector_scale(max_conc_term, -1. * te);
	}

	for (int ind = 0; ind != dim_proteins; ind++){
		if (max_concentration != 0.0){
			gsl_vector_set(g_this_step, ind, gsl_vector_get(evolution_vectorV, ind) + gsl_vector_get(max_conc_term, ind));
		}
		else{
			gsl_vector_set(g_this_step, ind, gsl_vector_get(evolution_vectorV, ind));
		}
	}
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(g_this_step, longlengthvector);

	//Create a vector to store two identical copies of Vloglength_vector (one for type A and one for type B).
	gsl_vector * V_logj_full = gsl_vector_alloc(dim_proteins);
	Define_Vloglength(V_logj_full, loglength, params, start_index_cells);
	gsl_vector_scale(V_logj_full, te);


	//Create a vector to store the right hand side of equation S12 from LUBENSKY PLOS 2012.
	gsl_vector * f_this_step = gsl_vector_alloc(dim_proteins);
	//Subtract the Vlogj for both type A and B
	gsl_vector_sub(f_this_step, V_logj_full);
	//Subtract the contribution for the U terms
	gsl_vector_sub(f_this_step, evolution_vectorU);
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(f_this_step, longlengthvector);

	//Now store g_this_step, f_this_step, y_this_step into c++ vector so that they can be used in future steps (third-order solver).
	for (int edge_index = 0; edge_index != dim_proteins; edge_index++){
		g_this_step_vector.at(edge_index) = gsl_vector_get(g_this_step, edge_index);
		f_this_step_vector.at(edge_index) = gsl_vector_get(f_this_step, edge_index);
		y_this_step_vector.at(edge_index) = gsl_vector_get(y_this_step, edge_index);
	}

	//Create a vector to store the right hand side of equation S12 from LUBENSKY PLOS 2012.
	gsl_vector * rhs_evolution = gsl_vector_alloc(dim_proteins);

	//Store past state of variables in a gsl vector.
	gsl_vector * y_nminus1 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(y_nminus1, index, y_step_before[1][index]);
	}

	//Store past state of explicit terms in a gsl vector.
	gsl_vector * f_nminus1 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(f_nminus1, index, f_step_before[1][index]);
	}

	//Store past state of implicit terms in a gsl vector.
	gsl_vector * g_nminus1 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(g_nminus1, index, g_step_before[1][index]);
	}

	//Consult third-order imex scheme for coefficient.
	gsl_vector_scale(y_nminus1, -1. * alpha0n / step_size * alphaP);
	//note lengthtimesconcentration is nothing but exp(-y)
	gsl_vector_mul(y_nminus1, lengthtimesconcentration);

	//Consult third-order imex scheme for coefficient.
	gsl_vector_scale(g_this_step, gamma1n);

	//Consult third-order imex scheme for coefficient.
	gsl_vector_scale(g_nminus1, gamma0n);
	gsl_vector * exponential_correction = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(exponential_correction, index, exp(y_step_before[1][index] - gsl_vector_get(y_this_step, index)));
	}
	gsl_vector_mul(g_nminus1, exponential_correction);

	//Consult third-order imex scheme for coefficient.
	gsl_vector_scale(f_this_step, beta1n);

	//Consult third-order imex scheme for coefficient.
	gsl_vector_scale(f_nminus1, beta0n);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(exponential_correction, index, exp(y_step_before[1][index] - gsl_vector_get(y_this_step, index)));
	}
	gsl_vector_mul(f_nminus1, exponential_correction);

	//Now scale y_this_step (cannot do before since it is used for exponential corrections)
	gsl_vector_scale(y_this_step, -1. * alpha1n / step_size * alphaP);
	//note lengthtimesconcentration is nothing but exp(-y)
	gsl_vector_mul(y_this_step, lengthtimesconcentration);

	//Store all of the terms for the right hand side of the evolution equation.
	gsl_vector_add(rhs_evolution, y_nminus1);
	gsl_vector_add(rhs_evolution, y_this_step);
	gsl_vector_add(rhs_evolution, g_nminus1);
	gsl_vector_add(rhs_evolution, g_this_step);
	gsl_vector_add(rhs_evolution, f_nminus1);
	gsl_vector_add(rhs_evolution, f_this_step);


	//Create a vector for storing the updated protein variable values.
	vector<double> updated_protein_variables;
	updated_protein_variables.resize(dim_proteins);

	//Note: one cannot use a cell iterator in a parallelized for loop, must build pointer list.
	vector<Cell*> CellPointerList;
	CellPointerList.reserve(type_dimension);
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		CellPointerList.push_back(&(*(*cit)));
	}

	//Parallelize the for loop (solve PCP dynamics on a cell-by-cell basis).
#pragma omp parallel for schedule(guided)															//BEGIN parallel for
	//not sure if parallelization is worth it here (probably not).
	for (int i = 0; i < CellPointerList.size(); i++){
		Cell* cptr = CellPointerList[i];
		PCPCell* pcpcptr = dynamic_cast<PCPCell*>(cptr);
		//If the cell is differentiated,
		if ((*pcpcptr).is_Cell_with_PCP()){
			//This struct will store parameters for numerical solving.
			struct rparams p;
			p.delta_t = step_size;
			p.alphastepforward = alpha2n;
			p.gammastepforward = gamma2n;
			int iter = 0;
			int status;
			int cell_id = (*cptr).id;
			int cell_index = cell_id2index.at(cell_id);
			int start_index = start_index_cells[cell_index];
			int completion_code = 27;
			//Store the list of edges for a given cell.
			//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcptr).FindOrderedEdges();
			//Store the number of edges within this cell.
			int number_of_edges = edge_list_for_this_cell.size();
			vector<double> edge_lengths;
			edge_lengths.resize(number_of_edges);
			//define a length vector
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				if (evolve_mech_positions == 1){
					//Use the new edge lengths for the next step term.
					edge_lengths.at(edge_index) = gsl_vector_get(lengthvect_new, edge_index + start_index);
				}
				else{
					edge_lengths.at(edge_index) = gsl_vector_get(lengthvect, edge_index + start_index);
				}
			}

			//Now assign this to rparams.
			p.length_edges = edge_lengths;

			gsl_vector_view rhs_evolution_viewA_this_cell = gsl_vector_subvector(rhs_evolution, start_index, number_of_edges);
			gsl_vector_view rhs_evolution_viewB_this_cell = gsl_vector_subvector(rhs_evolution, start_index + type_dimension, number_of_edges);
			//Then, store this view into an actual gsl vector.
			gsl_vector * rhs_evolution_type_A = &rhs_evolution_viewA_this_cell.vector;
			gsl_vector * rhs_evolution_type_B = &rhs_evolution_viewB_this_cell.vector;

			//Store the current log protein variables for protein of type A.
			vector<double> y_inA_this_cell;
			y_inA_this_cell.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				y_inA_this_cell.at(edge_index) = y.at(edge_index + start_index);
			}

			//Store the current log protein variables for protein of type B.
			vector<double> y_inB_this_cell;
			y_inB_this_cell.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				y_inB_this_cell.at(edge_index) = y.at(edge_index + start_index + type_dimension);
			}

			//Store the right hand side of the evolution of protein of type A.
			vector<double> rhs_evolution_A;
			rhs_evolution_A.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				rhs_evolution_A.at(edge_index) = gsl_vector_get(rhs_evolution_type_A, edge_index);
			}

			//Store the right hand side of the evolution of protein of type B.
			vector<double> rhs_evolution_B;
			rhs_evolution_B.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				rhs_evolution_B.at(edge_index) = gsl_vector_get(rhs_evolution_type_B, edge_index);
			}

			if (evolve_mech_positions == 1){
				Initialize_V_Matrix_Cell_New(params, start_index_cells, cptr);
			}
			else{
				Initialize_V_Matrix_Cell(params, start_index_cells, cptr);
			}

			//first solve for protein A, store necessary inputs for solver.
			p.V_matrix = (*pcpcptr).V_matrix;
			p.y_in = y_inA_this_cell;
			p.rhs_evaluation = rhs_evolution_A;

			//Initialize the solver and the function.
			gsl_multiroot_fdfsolver *s_A;
			gsl_multiroot_function_fdf f_A = { &numerical_solve_second_order, &numerical_solve_second_order_df, &numerical_solve_second_order_fdf, number_of_edges, &p };

			//Store the initial conditions for type A.
			gsl_vector * x_A = gsl_vector_alloc(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(x_A, edge_index, y_inA_this_cell.at(edge_index));
			}

			//Use a hybrid nonlinear root finder.
			nonlinear_solver_PCP_second_order = gsl_multiroot_fdfsolver_hybridj;
			//Allocate and set propertires of the solver.
			s_A = gsl_multiroot_fdfsolver_alloc(nonlinear_solver_PCP_second_order, number_of_edges);
			gsl_multiroot_fdfsolver_set(s_A, &f_A, x_A);

			double solver_tolerance = 1e-7;
			int max_iterations = 1000;
			try{
				//Iterate the solver until it reaches the appropriate tolerance or exceeds max iterations.
				do {

					status = gsl_multiroot_fdfsolver_iterate(s_A);

					if (status)
						break;

					status = gsl_multiroot_test_residual(s_A->f, solver_tolerance);
					iter++;
				}
				while (status == GSL_CONTINUE || iter < max_iterations);
				throw status;
			}
			//Catch error from the root finding loop.
			catch (int e){
				if (e != completion_code){
					printf("error: %s\n", gsl_strerror(e));
				}
			}

			//Store the updated protein variable values.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				updated_protein_variables.at(start_index + edge_index) = gsl_vector_get(s_A->x, edge_index);
			}

			//IMPORTANT NOTE: YOU MUST DEALLOCATE THE VARIABLES FROM THE TYPE A SOLVER. (CANNOT REUSE FOR SECOND->WILL CAUSE MEMORY LEAK).
			gsl_vector_free(x_A);
			gsl_multiroot_fdfsolver_free(s_A);

			//Store necessary parameters for B solver.
			p.y_in = y_inB_this_cell;
			p.rhs_evaluation = rhs_evolution_B;

			//Initialize the solver and the function.
			gsl_multiroot_fdfsolver *s_B;
			gsl_multiroot_function_fdf f_B = { &numerical_solve_second_order, &numerical_solve_second_order_df, &numerical_solve_second_order_fdf, number_of_edges, &p };

			//Store the initial conditions for type B.
			gsl_vector * x_B = gsl_vector_alloc(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(x_B, edge_index, y_inB_this_cell.at(edge_index));
			}
			//Use a hybrid nonlinear root finder.
			nonlinear_solver_PCP_second_order = gsl_multiroot_fdfsolver_hybridj;
			s_B = gsl_multiroot_fdfsolver_alloc(nonlinear_solver_PCP_second_order, number_of_edges);
			gsl_multiroot_fdfsolver_set(s_B, &f_B, x_B);

			iter = 0;
			try{
				//Iterate the solver until it reaches the appropriate tolerance or exceeds max iterations.
				do {
					status = gsl_multiroot_fdfsolver_iterate(s_B);

					if (status)
						break;

					status = gsl_multiroot_test_residual(s_B->f, solver_tolerance);
					iter++;
				}
				while (status == GSL_CONTINUE || iter < max_iterations);
				throw status;
			}
			//Catch error from the root finding loop.
			catch (int e){
				if (e != completion_code){
					printf("error: %s\n", gsl_strerror(e));
				}
			}

			//Store the updated protein variable values.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				updated_protein_variables.at(start_index + edge_index + type_dimension) = gsl_vector_get(s_B->x, edge_index);
			}

			//IMPORTANT NOTE: YOU MUST DEALLOCATE THE VARIABLES FROM THE TYPE A SOLVER.
			gsl_vector_free(x_B);
			gsl_multiroot_fdfsolver_free(s_B);
		}
	}


	//Free the memory associated with all of the matrices and vectors used in gsl. Don't allow for memory leaks.
	gsl_vector_free(rhs_evolution);
	gsl_vector_free(concentration_vector);
	gsl_vector_free(evolution_vectorU);
	gsl_vector_free(lengthvect);
	gsl_vector_free(lengthvect_new);
	gsl_vector_free(lengthtimesconcentration);
	gsl_vector_free(longlengthvector);
	gsl_vector_free(longlengthvector_new);
	gsl_vector_free(loglength);
	gsl_vector_free(V_logj_full);
	gsl_vector_free(evolution_vectorV);
	gsl_vector_free(max_conc_term);
	gsl_vector_free(y_this_step);
	gsl_vector_free(g_this_step);
	gsl_vector_free(f_this_step);
	gsl_vector_free(y_nminus1);
	gsl_vector_free(f_nminus1);
	gsl_vector_free(g_nminus1);
	gsl_vector_free(exponential_correction);

	//Store information from this step in the last column of the vectors that store information from past steps.
	tissue.f_step_before[2] = f_this_step_vector;
	tissue.y_step_before[2] = y_this_step_vector;
	tissue.g_step_before[2] = g_this_step_vector;
	tissue.past_step_sizes[2] = step_size;

	//Then, update the EvalCountProtein for purposes of documenting efficiency.
	(*(combinedparameters*)params).EvalCountProtein++;

	return updated_protein_variables;
}

//This is the full function to numerically solve the PCP dynamics to third order.
vector<double> imex_solver_PCP_third_order(const vector<double> &y, double step_size, void* params, ofstream* errorfile){
	//other choice for coefficients:
	//consult https://www.ideals.illinois.edu/bitstream/handle/2142/18708/vsimex_jcm_ideals.pdf?sequence=2
	//also consult Ascher 1995 (IMEX fixed step size)
	double gamma = 0.5;
	double theta = -1.25;
	double c = -0.52;

	//Use third order step based on imex scheme (where the implicit variables are strictly zero).
	//Also, store information from this step.
	if (evolve_mech_positions == 1){
		helper_imex_third_order_mech(params, step_size, errorfile);
	}

	//info from previous time steps for higher order
	vector<vector<double>> g_step_before = tissue.g_step_before;
	vector<vector<double>> f_step_before = tissue.f_step_before;
	vector<vector<double>> y_step_before = tissue.y_step_before;

	//Changes by Hayden for PCP: Please note that now that we are evolving protein concentrations we need to store the number of vertices and edges
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	//U is composed of four submatrices (A-A, A-B, B-A, B-B). Define this variable to store the number of protein concentrations variables for a given type (A or B).
	int type_dimension = dim_proteins / 2;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;

	//for storing information from this step.
	vector<double> g_this_step_vector;
	g_this_step_vector.resize(dim_proteins);
	vector<double> f_this_step_vector;
	f_this_step_vector.resize(dim_proteins);
	vector<double> y_this_step_vector;
	y_this_step_vector.resize(dim_proteins);

	//Initialize loop variables.
	int cell_id;
	int cell_index;
	int start_index;

	//This vector is important. It allows you to figure out what row or column index refers to a given cell's edges in the loop for defining U below.
	//Given some cell index, this vector will give you where the entries for this cell begin.
	vector<int> start_index_cells;
	vector<PCPEdge*> edge_list_for_this_cell;
	//Resize this vector such that it fits all cells in the tissue.
	start_index_cells.resize(tissue.cells.size());
	//Initialize start_index to zero.
	start_index = 0;
	//Use ind to refer to the current cell index.
	cell_index = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){

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

	//Store the lengths for all edges.
	gsl_vector * lengthvect = gsl_vector_alloc(type_dimension);
	gsl_vector * lengthvect_new = gsl_vector_alloc(type_dimension);
	int ind = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			int cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				gsl_vector_set(lengthvect, ind, (*edge_list_for_this_cell[edge_index]).GetLength());
				gsl_vector_set(lengthvect_new, ind, (*edge_list_for_this_cell[edge_index]).GetLength_New());
				ind++;
			}
		}
	}
	//This may seem redundant, but it is useful. Create a vector which is two copies of lengthvectforjacob stacked on each other.
	gsl_vector * longlengthvector = gsl_vector_alloc(dim_proteins);
	gsl_vector * longlengthvector_new = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(longlengthvector, ind, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector, ind + type_dimension, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector_new, ind, gsl_vector_get(lengthvect_new, ind));
		gsl_vector_set(longlengthvector_new, ind + type_dimension, gsl_vector_get(lengthvect_new, ind));
	}

	//Define the dimension of the matrix U.
	int U_dimension = dim_proteins;
	//Allocate the U-matrix as a gsl_matrix.
	gsl_matrix * U = gsl_matrix_alloc(U_dimension, U_dimension);
	//Store the concentrations of edges that correspond to row/column indices of U.
	gsl_vector * concentration_vector = gsl_vector_alloc(U_dimension);
	//Store the vector that defines Uijcj for a given protein index i.
	gsl_vector * evolution_vectorU = gsl_vector_alloc(U_dimension);

	//Initialize the U matrix.
	Initialize_U_Matrix(U, params, start_index_cells);

	//Define the vector of concentrations for multiplying the U matrix.
	Define_Concentration_Vector(concentration_vector, params, y, longlengthvector);

	//Multiply Uijcj->Build an evolution vector.
	gsl_blas_dgemv(CblasNoTrans, 1.0, U, concentration_vector, 0.0, evolution_vectorU);

	//Only need 2 * number of edges entries in this matrix (cell1 and cell2 for every edge) 
	//log terms only couple protein concentrations of the same type (specifying a 4 * number of edges X 4 * number of edges matrix would be wasteful).
	//Initialize a gsl_matrix for the V matrix (consult PLOS 2012 supplementary information).
	gsl_matrix * V = gsl_matrix_alloc(type_dimension, type_dimension);
	gsl_matrix * V_new = gsl_matrix_alloc(type_dimension, type_dimension);

	//Initialize the V matrix (coupling edges within the same cell->same type of protein)
	Initialize_V_Matrix(V, params, start_index_cells);
	//Initialize the V matrix based on new edge lengths.
	Initialize_V_Matrix_New(V_new, params, start_index_cells);

	//Define a vector for storing length*concentration=protein amount for type A and B (A, then B).
	gsl_vector * lengthtimesconcentration = gsl_vector_alloc(dim_proteins);
	//Copy longlengthvector into lengthtimesconcentration.
	gsl_vector_memcpy(lengthtimesconcentration, longlengthvector);
	//This multiplies concentration by length elementwise and stores in lengthtimesconcentration
	gsl_vector_mul(lengthtimesconcentration, concentration_vector);

	//Store the log of the length of each edge.
	gsl_vector * loglength = gsl_vector_alloc(type_dimension);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(loglength, ind, log(gsl_vector_get(lengthvect, ind)));
	}


	//Define ratios of step sizes.
	double omega_n_plus_2 = step_size / tissue.past_step_sizes[1];
	double omega_n_plus_1 = tissue.past_step_sizes[1] / tissue.past_step_sizes[0];

	//Define coefficients for third order equations.
	double alpha0n = -1. * pow(omega_n_plus_2, 2.) * pow(omega_n_plus_1, 3.) * (3. * pow(gamma, 2.) * omega_n_plus_2 + 2. * gamma * (1. - omega_n_plus_2) - 1.) / (1. + omega_n_plus_1) / (1. + omega_n_plus_1 * (1. + omega_n_plus_2));
	double alpha1n = -1. * pow(omega_n_plus_2,2.) * (gamma * omega_n_plus_1 * omega_n_plus_2 * (2. - 3. * gamma) + (1. - 2. * gamma) * (1. + omega_n_plus_1)) / (1. + omega_n_plus_2);
	double alpha2n = -1. * (omega_n_plus_1 * (1. + gamma*omega_n_plus_2) * (1. + omega_n_plus_2 * (3. * gamma - 2.)) + omega_n_plus_2 * (2. * gamma - 1.) + 1.) / (1. + omega_n_plus_1) - theta;
	double alpha3n = (1. + 2. * gamma*omega_n_plus_2 + omega_n_plus_1 * (1. + omega_n_plus_2*gamma)*(1. + 3. * gamma*omega_n_plus_2)) / (1. + omega_n_plus_2) / (1. + omega_n_plus_1*(1. + omega_n_plus_2)) + theta;
	double beta0n = pow(omega_n_plus_1,2.) * omega_n_plus_2 * (6. * gamma * (1. + gamma*omega_n_plus_2) + theta * (3. + 2. * omega_n_plus_2)) / 6. / (1. + omega_n_plus_1);
	double beta1n = -gamma * omega_n_plus_2 * (1. + omega_n_plus_1 * (1. + gamma*omega_n_plus_2)) - (theta * omega_n_plus_2 * (3. + omega_n_plus_1 * (3. + 2. * omega_n_plus_2))) / 6.;
	double beta2n = (1. + gamma * omega_n_plus_2) * (1. + omega_n_plus_1 * (1. + gamma * omega_n_plus_2)) / (1. + omega_n_plus_1) + theta * (1. + omega_n_plus_2 / 2. + omega_n_plus_2 * omega_n_plus_1 * (3. + 2. * omega_n_plus_2) / 6. / (1. + omega_n_plus_1));
	double gamma0n = theta * pow(omega_n_plus_1,2.) * omega_n_plus_2 * (3. + 2. * omega_n_plus_2) / 6. / (1. + omega_n_plus_1) - c;
	double gamma1n = (c * (1. + omega_n_plus_1) * (1. + omega_n_plus_1 * (1. + omega_n_plus_2)) - pow(omega_n_plus_1,2.) * pow(omega_n_plus_2,2.) * gamma * (1. - gamma)) / pow(omega_n_plus_1,2.) / (1. + omega_n_plus_2) - theta * omega_n_plus_2 * (3. + omega_n_plus_1 * (3. + 2. * omega_n_plus_2)) / 6.;
	double gamma2n = (pow(omega_n_plus_1,2.) * omega_n_plus_2 * (1. - gamma) * (1. + gamma*omega_n_plus_2) - c * (1. + omega_n_plus_1 * (1. + omega_n_plus_2))) / pow(omega_n_plus_1,2.) / omega_n_plus_2 + theta * (1. + omega_n_plus_2 / 2. + omega_n_plus_2 * omega_n_plus_1 * (3. + 2. * omega_n_plus_2) / 6. / (1. + omega_n_plus_1));
	double gamma3n = (pow(omega_n_plus_1,2.) * omega_n_plus_2 * gamma * (1. + gamma*omega_n_plus_2) + c * (1. + omega_n_plus_1)) / pow(omega_n_plus_1,2.) / omega_n_plus_2 / (1. + omega_n_plus_2);


	/*********This is where I begin building the right hand side of the third-order scheme*********/
	//Please note that GSL vector conventions forces me to directly multiply the vector with necessary factor.
	//When reading this code, pay attention to the gsl_vector_mul, gsl_vector_scale commands.

	//g contains the implicitly handled term. f handles the explicitly handled term.

	//Will store the current -log variable values.
	gsl_vector * y_this_step = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != dim_proteins; ind++){
		gsl_vector_set(y_this_step, ind, y.at(ind));
	}

	//Initialize g_this_step which we store -te*V*X_vector for both types of proteins.
	gsl_vector * g_this_step = gsl_vector_alloc(dim_proteins);
	gsl_vector * evolution_vectorV = gsl_vector_alloc(dim_proteins);
	Define_Evolution_Vector_V(evolution_vectorV, y_this_step, params, start_index_cells);
	//Scale by -Te.
	gsl_vector_scale(evolution_vectorV, -1. * te);

	//Add max concentration term.
	gsl_vector * max_conc_term = gsl_vector_alloc(dim_proteins);
	if (max_concentration != 0.0){
		Define_Evolution_Vector_Max_Conc(max_conc_term, y_this_step, lengthvect, params, start_index_cells);
		//Scale by -Te.
		gsl_vector_scale(max_conc_term, -1. * te);
	}

	for (int ind = 0; ind != dim_proteins; ind++){
		if (max_concentration != 0.0){
			gsl_vector_set(g_this_step, ind, gsl_vector_get(evolution_vectorV, ind) + gsl_vector_get(max_conc_term, ind));
		}
		else{
			gsl_vector_set(g_this_step, ind, gsl_vector_get(evolution_vectorV, ind));
		}
	}
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(g_this_step, longlengthvector);

	//The evolution vector will store sum(Vij*Xj) for protein A.
	gsl_vector * Vlogj_vector = gsl_vector_alloc(type_dimension);
	//look this up to make sure that it is doing what you want.
	gsl_blas_dgemv(CblasNoTrans, 1.0, V, loglength, 0.0, Vlogj_vector);
	//Scale this vector by TE
	gsl_vector_scale(Vlogj_vector, te);

	//Create a vector to store two identical copies of Vlogj_vector (one for type A and one for type B).
	gsl_vector * V_logj_full = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(V_logj_full, ind, gsl_vector_get(Vlogj_vector, ind));
		gsl_vector_set(V_logj_full, ind + type_dimension, gsl_vector_get(Vlogj_vector, ind));
	}


	//Create a vector to store the right hand side of equation S12 from LUBENSKY PLOS 2012.
	gsl_vector * f_this_step = gsl_vector_alloc(dim_proteins);
	//Subtract the Vlogj for both type A and B
	gsl_vector_sub(f_this_step, V_logj_full);
	//Subtract the contribution for the U terms
	gsl_vector_sub(f_this_step, evolution_vectorU);
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(f_this_step, longlengthvector);

	//Now store g_this_step, f_this_step, y_this_step into c++ vector so that they can be used in future steps (third-order solver).
	for (int edge_index = 0; edge_index != dim_proteins; edge_index++){
		g_this_step_vector.at(edge_index) = gsl_vector_get(g_this_step, edge_index);
		f_this_step_vector.at(edge_index) = gsl_vector_get(f_this_step, edge_index);
		y_this_step_vector.at(edge_index) = gsl_vector_get(y_this_step, edge_index);
	}

	//Create a vector to store the right hand side of equation S12 from LUBENSKY PLOS 2012.
	gsl_vector * rhs_evolution = gsl_vector_alloc(dim_proteins);

	//Convert the information from past steps into gsl_vectors
	gsl_vector * y_nminus2 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(y_nminus2, index, y_step_before[0][index]);
	}

	gsl_vector * y_nminus1 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(y_nminus1, index, y_step_before[1][index]);
	}

	gsl_vector * f_nminus2 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(f_nminus2, index, f_step_before[0][index]);
	}

	gsl_vector * f_nminus1 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(f_nminus1, index, f_step_before[1][index]);
	}

	gsl_vector * g_nminus2 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(g_nminus2, index, g_step_before[0][index]);
	}

	gsl_vector * g_nminus1 = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(g_nminus1, index, g_step_before[1][index]);
	}

	//Now, proceed to scale each of these terms with their corresponding factors and multiply by an exponential correction (CONSULT PSEUDOCODE EXPLANATION).
	gsl_vector_scale(y_nminus2, -1. * alpha0n / step_size * alphaP);
	//note lengthtimesconcentration is nothing but exp(-y)
	gsl_vector_mul(y_nminus2, lengthtimesconcentration);

	gsl_vector_scale(y_nminus1, -1. * alpha1n / step_size * alphaP);
	//note lengthtimesconcentration is nothing but exp(-y)
	gsl_vector_mul(y_nminus1, lengthtimesconcentration);
	
	gsl_vector_scale(g_this_step, gamma2n);

	gsl_vector_scale(g_nminus1, gamma1n);
	gsl_vector * exponential_correction = gsl_vector_alloc(dim_proteins);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(exponential_correction, index, exp(y_step_before[1][index] - gsl_vector_get(y_this_step, index)));
	}
	gsl_vector_mul(g_nminus1, exponential_correction);

	gsl_vector_scale(g_nminus2, gamma0n);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(exponential_correction, index, exp(y_step_before[0][index] - gsl_vector_get(y_this_step, index)));
	}
	gsl_vector_mul(g_nminus2, exponential_correction);

	gsl_vector_scale(f_this_step, beta2n); 

	gsl_vector_scale(f_nminus1, beta1n);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(exponential_correction, index, exp(y_step_before[1][index] - gsl_vector_get(y_this_step, index)));
	}
	gsl_vector_mul(f_nminus1, exponential_correction);

	gsl_vector_scale(f_nminus2, beta0n);
	for (int index = 0; index < dim_proteins; index++){
		gsl_vector_set(exponential_correction, index, exp(y_step_before[0][index] - gsl_vector_get(y_this_step, index)));
	}
	gsl_vector_mul(f_nminus2, exponential_correction);

	//Now scale y_this_step (cannot do before since it is used for exponential corrections)
	gsl_vector_scale(y_this_step, -1. * alpha2n / step_size * alphaP);
	//note lengthtimesconcentration is nothing but exp(-y)
	gsl_vector_mul(y_this_step, lengthtimesconcentration);

	//Store all of the terms for the right hand side of the evolution equation.
	gsl_vector_add(rhs_evolution, y_nminus2);
	gsl_vector_add(rhs_evolution, y_nminus1);
	gsl_vector_add(rhs_evolution, y_this_step);
	gsl_vector_add(rhs_evolution, g_nminus2);
	gsl_vector_add(rhs_evolution, g_nminus1);
	gsl_vector_add(rhs_evolution, g_this_step);
	gsl_vector_add(rhs_evolution, f_nminus2);
	gsl_vector_add(rhs_evolution, f_nminus1);
	gsl_vector_add(rhs_evolution, f_this_step);
	

	//Create a vector for storing the updated protein variable values.
	vector<double> updated_protein_variables;
	updated_protein_variables.resize(dim_proteins);

	//Note: one cannot use a cell iterator in a parallelized for loop, must build pointer list.
	vector<Cell*> CellPointerList;
	CellPointerList.reserve(type_dimension);
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		CellPointerList.push_back(&(*(*cit)));
	}


	//Parallelize the for loop (solve PCP dynamics on a cell-by-cell basis).
#pragma omp parallel for schedule(guided)															//BEGIN parallel for
	//not sure if parallelization is worth it here (probably not).
	for (int i = 0; i < CellPointerList.size(); i++){
		Cell* cptr = CellPointerList[i];
		PCPCell* pcpcptr = dynamic_cast<PCPCell*>(cptr);
		//If the cell is differentiated,
		if ((*pcpcptr).is_Cell_with_PCP()){
			//This struct will store parameters for numerical solving.
			struct rparams p;
			p.delta_t = step_size;
			p.alphastepforward = alpha3n;
			p.gammastepforward = gamma3n;
			int iter = 0;
			int status;
			int cell_id = (*cptr).id;
			int cell_index = cell_id2index.at(cell_id);
			int start_index = start_index_cells[cell_index];
			int completion_code = 27;
			//Store the list of edges for a given cell.
			//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcptr).FindOrderedEdges();
			//Store the number of edges within this cell.
			int number_of_edges = edge_list_for_this_cell.size();
			vector<double> edge_lengths;
			edge_lengths.resize(number_of_edges);
			//define a length vector
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				if (evolve_mech_positions == 1){
					//Use the new edge lengths for the next step term.
					edge_lengths.at(edge_index) = gsl_vector_get(lengthvect_new, edge_index + start_index);
				}
				else{
					edge_lengths.at(edge_index) = gsl_vector_get(lengthvect, edge_index + start_index);
				}
			}

			//Now assign this to rparams.
			p.length_edges = edge_lengths;

			gsl_vector_view rhs_evolution_viewA_this_cell = gsl_vector_subvector(rhs_evolution, start_index, number_of_edges);
			gsl_vector_view rhs_evolution_viewB_this_cell = gsl_vector_subvector(rhs_evolution, start_index + type_dimension, number_of_edges);

			//Store the current log protein variables for protein of type A.
			vector<double> y_inA_this_cell;
			y_inA_this_cell.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				y_inA_this_cell.at(edge_index) = y.at(edge_index + start_index);
			}
			
			//Store the current log protein variables for protein of type B.
			vector<double> y_inB_this_cell;
			y_inB_this_cell.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				y_inB_this_cell.at(edge_index) = y.at(edge_index + start_index + type_dimension);
			}

			//Store the right hand side of the evolution of protein of type A.
			vector<double> rhs_evolution_A;
			rhs_evolution_A.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				rhs_evolution_A.at(edge_index) = gsl_vector_get(&rhs_evolution_viewA_this_cell.vector, edge_index);
			}

			//Store the right hand side of the evolution of protein of type B.
			vector<double> rhs_evolution_B;
			rhs_evolution_B.resize(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				rhs_evolution_B.at(edge_index) = gsl_vector_get(&rhs_evolution_viewB_this_cell.vector, edge_index);
			}

			//extract a submatrix of V
			gsl_matrix_view V_submatrix_view;
			if (evolve_mech_positions == 1){
				V_submatrix_view = gsl_matrix_submatrix(V_new, start_index, start_index, number_of_edges, number_of_edges);
			}
			else{
				V_submatrix_view = gsl_matrix_submatrix(V, start_index, start_index, number_of_edges, number_of_edges);
			}

			//Store V in a C++ matrix  (not absolutely necessary and will probably change).
			double init_val = 0.0;
			vector<vector<double>> vec(number_of_edges, vector<double>(number_of_edges, init_val));
			//Use the V matrix based on new edge lengths.
			vector<vector<double>> V_matrix = vec;
			for (int first_edge_index = 0; first_edge_index != number_of_edges; first_edge_index++){
				for (int second_edge_index = 0; second_edge_index != number_of_edges; second_edge_index++){
					V_matrix[first_edge_index][second_edge_index] = gsl_matrix_get(&V_submatrix_view.matrix, first_edge_index, second_edge_index);
				}
			}

			//first solve for protein A, store necessary inputs for solver.
			p.V_matrix = V_matrix;
			p.y_in = y_inA_this_cell;
			p.rhs_evaluation = rhs_evolution_A;

			//Initialize the solver and the function.
			gsl_multiroot_fdfsolver *s_A;
			gsl_multiroot_function_fdf f_A = { &numerical_solve_third_order, &numerical_solve_third_order_df, &numerical_solve_third_order_fdf, number_of_edges, &p };
			
			//Store the initial conditions for type B.
			gsl_vector * x_A = gsl_vector_alloc(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(x_A, edge_index, y_inA_this_cell.at(edge_index));
			}

			//Use a hybrid nonlinear root finder.
			nonlinear_solver_PCP_third_order = gsl_multiroot_fdfsolver_hybridj;
			//Allocate and set propertires of the solver.
			s_A = gsl_multiroot_fdfsolver_alloc(nonlinear_solver_PCP_third_order, number_of_edges);
			gsl_multiroot_fdfsolver_set(s_A, &f_A, x_A);


			double solver_tolerance = 1e-10;
			int max_iterations = 1000;
			try{
				//Iterate the solver until it reaches the appropriate tolerance or exceeds max iterations.
				do {

				status = gsl_multiroot_fdfsolver_iterate(s_A);

				if (status)
					break;

				status = gsl_multiroot_test_residual(s_A->f, solver_tolerance);
				iter++;
				}
				while (status == GSL_CONTINUE || iter < max_iterations);
				throw status;
			}
			//Catch error from the root finding loop.
			catch (int e){
				if (e != completion_code){
					printf("error: %s\n", gsl_strerror(e));
				}
			}

			//Store the updated protein variable values.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				updated_protein_variables.at(start_index + edge_index) = gsl_vector_get(s_A->x, edge_index);
			}

			//IMPORTANT NOTE: YOU MUST DEALLOCATE THE VARIABLES FROM THE TYPE A SOLVER. (CANNOT REUSE FOR SECOND->WILL CAUSE MEMORY LEAK).
			gsl_vector_free(x_A);
			gsl_multiroot_fdfsolver_free(s_A);

			//Store necessary parameters for B solver.
			p.y_in = y_inB_this_cell;
			p.rhs_evaluation = rhs_evolution_B;

			//Initialize the solver and the function.
			gsl_multiroot_fdfsolver *s_B;
			gsl_multiroot_function_fdf f_B = { &numerical_solve_third_order, &numerical_solve_third_order_df, &numerical_solve_third_order_fdf, number_of_edges, &p };
			
			//Store the initial conditions for type B.
			gsl_vector * x_B = gsl_vector_alloc(number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_vector_set(x_B, edge_index, y_inB_this_cell.at(edge_index));
			}
			//Use a hybrid nonlinear root finder.
			nonlinear_solver_PCP_third_order = gsl_multiroot_fdfsolver_hybridj;
			s_B = gsl_multiroot_fdfsolver_alloc(nonlinear_solver_PCP_third_order, number_of_edges);
			gsl_multiroot_fdfsolver_set(s_B, &f_B, x_B);

			iter = 0;
			try{
				//Iterate the solver until it reaches the appropriate tolerance or exceeds max iterations.
				do {
					status = gsl_multiroot_fdfsolver_iterate(s_B);

					if (status)
						break;

					status = gsl_multiroot_test_residual(s_B->f, solver_tolerance);
					iter++;
				}
				while (status == GSL_CONTINUE || iter < max_iterations);
				throw status;
			}
			//Catch error from the root finding loop.
			catch (int e){
				if (e != completion_code){
					printf("error: %s\n", gsl_strerror(e));
				}
			}

			//Store the updated protein variable values.
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				updated_protein_variables.at(start_index + edge_index + type_dimension) = gsl_vector_get(s_B->x, edge_index);
			}
			
			//IMPORTANT NOTE: YOU MUST DEALLOCATE THE VARIABLES FROM THE TYPE A SOLVER.
			gsl_vector_free(x_B);
			gsl_multiroot_fdfsolver_free(s_B);
		}
	}

	//Free the memory associated with all of the matrices and vectors used in gsl. Don't allow for memory leaks.
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(V_new);
	gsl_vector_free(rhs_evolution);
	gsl_vector_free(concentration_vector);
	gsl_vector_free(evolution_vectorU);
	gsl_vector_free(lengthvect);
	gsl_vector_free(lengthvect_new);
	gsl_vector_free(lengthtimesconcentration);
	gsl_vector_free(longlengthvector);
	gsl_vector_free(longlengthvector_new);
	gsl_vector_free(loglength);
	gsl_vector_free(Vlogj_vector);
	gsl_vector_free(V_logj_full);
	gsl_vector_free(evolution_vectorV);
	gsl_vector_free(max_conc_term);
	gsl_vector_free(y_this_step);
	gsl_vector_free(g_this_step);
	gsl_vector_free(f_this_step);
	gsl_vector_free(y_nminus2);
	gsl_vector_free(y_nminus1);
	gsl_vector_free(f_nminus2);
	gsl_vector_free(f_nminus1);
	gsl_vector_free(g_nminus2);
	gsl_vector_free(g_nminus1);
	gsl_vector_free(exponential_correction);

	//Store information from this step in the last column of the vectors that store information from past steps.
	tissue.f_step_before[2] = f_this_step_vector;
	tissue.y_step_before[2] = y_this_step_vector;
	tissue.g_step_before[2] = g_this_step_vector;
	tissue.past_step_sizes[2] = step_size;

	return updated_protein_variables;
}

//! This solves the differential equation for the PCP protein dynamics using the IMEX scheme in the supplementary info of LUBENSKY PLOS 2012.
//! It is important to note that this scheme is still under investigation for whether it does a reasonable job of describing the dynamics (i.e., minimizing free energy).
vector<double> imex_solver_guillaume_PCP(const vector<double> &y, double step_size, void* params, ofstream* errorfile){
	//Evolve the mechanical variables (based on euler explicit step). Also, store information from this step.
	if (evolve_mech_positions == 1){
		helper_imex_first_order_mech(params, step_size, errorfile);
	}

	//Store the number of protein variables for evolution (protein A and protein B for each edge in a differentiated cell).
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	//U is composed of four submatrices (A-A, A-B, B-A, B-B). Define this variable to store the number of protein concentrations variables for a given type (A or B).
	int type_dimension = dim_proteins / 2;
	unordered_map<long, int> cell_id2index = (*(combinedparameters*)params).cell_id2index;

	//Initialize loop variables.
	int cell_id;
	int cell_index;
	int start_index;

	//This vector is important. It allows you to figure out what row or column index refers to a given cell's edges in the loop for defining U below.
	//Given some cell index, this vector will give you where the entries for this cell begin.
	vector<int> start_index_cells;
	vector<PCPEdge*> edge_list_for_this_cell;
	//Resize this vector such that it fits all cells in the tissue.
	start_index_cells.resize(tissue.cells.size());
	//Initialize start_index to zero.
	start_index = 0;
	//Use ind to refer to the current cell index.
	cell_index = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){

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

	//Store the lengths for all edges.
	gsl_vector * lengthvect = gsl_vector_alloc(type_dimension);
	gsl_vector * lengthvect_new = gsl_vector_alloc(type_dimension);
	int ind = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			int cell_id = (*(*cit)).id;
			cell_index = cell_id2index.at(cell_id);
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			for (int edge_index = 0; edge_index != edge_list_for_this_cell.size(); edge_index++){
				gsl_vector_set(lengthvect, ind, (*edge_list_for_this_cell[edge_index]).GetLength());
				gsl_vector_set(lengthvect_new, ind, (*edge_list_for_this_cell[edge_index]).GetLength_New());
				ind++;
			}
		}
	}
	//This may seem redundant, but it is useful. Create a vector which is two copies of lengthvectforjacob stacked on each other.
	gsl_vector * longlengthvector = gsl_vector_alloc(dim_proteins);
	gsl_vector * longlengthvector_new = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(longlengthvector, ind, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector, ind + type_dimension, gsl_vector_get(lengthvect, ind));
		gsl_vector_set(longlengthvector_new, ind, gsl_vector_get(lengthvect_new, ind));
		gsl_vector_set(longlengthvector_new, ind + type_dimension, gsl_vector_get(lengthvect_new, ind));
	}

	//Define the dimension of the matrix U.
	int U_dimension = dim_proteins;
	//Allocate the U-matrix as a gsl_matrix.
	gsl_matrix * U = gsl_matrix_calloc(U_dimension, U_dimension);
	//Store the concentrations of edges that correspond to row/column indices of U.
	gsl_vector * concentration_vector = gsl_vector_calloc(U_dimension);
	//Store the vector that defines Uijcj for a given protein index i.
	gsl_vector * evolution_vectorU = gsl_vector_calloc(U_dimension);

	//Initialize the U matrix.
	Initialize_U_Matrix(U, params, start_index_cells);

	//Define concentration vector for multiplying the U matrix.
	Define_Concentration_Vector(concentration_vector, params, y, longlengthvector);

	//Multiply Uijcj->Build an evolution vector.
	gsl_blas_dgemv(CblasNoTrans, 1.0, U, concentration_vector, 0.0, evolution_vectorU);

	//Only need 2 * number of edges entries in this matrix (cell1 and cell2 for every edge) 
	//log terms only couple protein concentrations of the same type (specifying a 4 * number of edges X 4 * number of edges matrix would be wasteful).
	//Initialize a gsl_matrix for the V matrix (consult PLOS 2012 supplementary information).
	gsl_matrix * V = gsl_matrix_alloc(type_dimension, type_dimension);
	gsl_matrix * V_new = gsl_matrix_alloc(type_dimension, type_dimension);

	//Initialize the V matrix (coupling edges within the same cell->same type of protein)
	Initialize_V_Matrix(V, params, start_index_cells);
	//Intialize V matrix using new edge length.
	Initialize_V_Matrix_New(V_new, params, start_index_cells);

	//Define a vector for storing length*concentration=protein amount for type A and B (A, then B).
	gsl_vector * lengthtimesconcentration = gsl_vector_calloc(dim_proteins);
	//Copy longlengthvector into lengthtimesconcentration.
	gsl_vector_memcpy(lengthtimesconcentration, longlengthvector);
	//This multiplies concentration by length elementwise and stores in lengthtimesconcentration
	gsl_vector_mul(lengthtimesconcentration, concentration_vector);

	//Store the log of the length of each edge.
	gsl_vector * loglength = gsl_vector_alloc(type_dimension);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(loglength, ind, log(gsl_vector_get(lengthvect, ind)));
	}


	//Will store the current -log variable values.
	gsl_vector * y_this_step = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != dim_proteins; ind++){
		gsl_vector_set(y_this_step, ind, y.at(ind));
	}

	//Initialize g_this_step which we store -te*V*X_vector for both types of proteins.
	gsl_vector * g_this_step = gsl_vector_alloc(dim_proteins);
	gsl_vector * evolution_vectorV = gsl_vector_alloc(dim_proteins);
	Define_Evolution_Vector_V(evolution_vectorV, y_this_step, params, start_index_cells);
	//Scale by -Te.
	gsl_vector_scale(evolution_vectorV, -1. * te);

	//Add max concentration term.
	gsl_vector * max_conc_term = gsl_vector_alloc(dim_proteins);
	Define_Evolution_Vector_Max_Conc(max_conc_term, y_this_step, lengthvect, params, start_index_cells);
	//Scale by -Te.
	gsl_vector_scale(max_conc_term, -1. * te);

	for (int ind = 0; ind != dim_proteins; ind++){
		gsl_vector_set(g_this_step, ind, gsl_vector_get(evolution_vectorV, ind) + gsl_vector_get(max_conc_term, ind));
	}
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(g_this_step, longlengthvector);

	//The evolution vector will store sum(Vij*Xj) for protein A.
	gsl_vector * Vlogj_vector = gsl_vector_alloc(type_dimension);
	//look this up to make sure that it is doing what you want.
	gsl_blas_dgemv(CblasNoTrans, 1.0, V, loglength, 0.0, Vlogj_vector);
	//Scale this vector by TE
	gsl_vector_scale(Vlogj_vector, te);

	//Create a vector to store two identical copies of Vlogj_vector (one for type A and one for type B).
	gsl_vector * V_logj_full = gsl_vector_alloc(dim_proteins);
	for (int ind = 0; ind != type_dimension; ind++){
		gsl_vector_set(V_logj_full, ind, gsl_vector_get(Vlogj_vector, ind));
		gsl_vector_set(V_logj_full, ind + type_dimension, gsl_vector_get(Vlogj_vector, ind));
	}

	//Store the explicit terms.
	gsl_vector * f_this_step = gsl_vector_alloc(dim_proteins);
	//Subtract the Vlogj for both type A and B
	gsl_vector_sub(f_this_step, V_logj_full);
	//Subtract the contribution for the U terms
	gsl_vector_sub(f_this_step, evolution_vectorU);
	//Multiply every element by the corresponding edge length.
	gsl_vector_mul(f_this_step, longlengthvector);

	//Create a vector to store the right hand side of equation. (CONSULT PSEUDOCODE).
	gsl_vector * rhs_evolution = gsl_vector_alloc(dim_proteins);

	//Now scale y_this_step (cannot do before since it is used for exponential corrections)
	gsl_vector_scale(y_this_step, 1. / step_size * alphaP);
	//note lengthtimesconcentration is nothing but exp(-y)
	gsl_vector_mul(y_this_step, lengthtimesconcentration);

	//Add y_this_step to f_this_step for right hand side (this is for a forward-backward euler combo).
	gsl_vector_add(rhs_evolution, y_this_step);
	gsl_vector_add(rhs_evolution, f_this_step);

	//Create a vector for storing the updated protein variable values.
	vector<double> updated_protein_variables;
	updated_protein_variables.resize(dim_proteins);

	//Loop through and build the inverse block matrix for each cell for each type of protein. Multiply corresponding elements of rhs_12 by the block matrix for each type and store the updated protein values.
	start_index = 0;
	double current_matrix_element;
	cell_index = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++){
		PCPCell* pcpcit = dynamic_cast<PCPCell*>(*cit);
		if ((*pcpcit).is_Cell_with_PCP()){
			//Store the list of edges for a given cell.
			//order edges: because the J1 and J2 terms are local couplings, need to store an order list of edges. Find this list from ordered vertices.
			vector<PCPEdge*> edge_list_for_this_cell = (*pcpcit).FindOrderedEdges();
			//Store the number of edges within this cell.
			int number_of_edges = edge_list_for_this_cell.size();
			//Initialize inverse block matrices for this cell for both types of proteins.
			gsl_matrix * inverse_block_matrix_A = gsl_matrix_calloc(number_of_edges, number_of_edges);
			gsl_matrix * inverse_block_matrix_B = gsl_matrix_calloc(number_of_edges, number_of_edges);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				gsl_matrix_set(inverse_block_matrix_A, edge_index, edge_index, 1 / (exp(-y.at(edge_index + start_index)) / step_size * alphaP + gsl_vector_get(lengthvect, edge_index + start_index) * te));
				gsl_matrix_set(inverse_block_matrix_B, edge_index, edge_index, 1 / (exp(-y.at(edge_index + start_index + type_dimension)) / step_size * alphaP + gsl_vector_get(lengthvect, edge_index + start_index) * te));
			}
			//Define a variable for storing the sum of lengths of the edges.
			double length_sum = 0.0;
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				length_sum = length_sum + gsl_vector_get(lengthvect, edge_index + start_index);
			}
			//Compute the sum of lk/(ck/deltat+te) for the protein of type A for this cell.
			double matrix_prefactor_sum_A = 0.0;
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				matrix_prefactor_sum_A = matrix_prefactor_sum_A + pow(gsl_vector_get(lengthvect, start_index + edge_index), 2) / (exp(-y.at(edge_index + start_index)) / step_size * alphaP + gsl_vector_get(lengthvect, edge_index + start_index) * te);
			}
			//Compute the sum of lk/(ck/deltat+te) for the protein of type B for this cell.
			double matrix_prefactor_sum_B = 0.0;
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				matrix_prefactor_sum_B = matrix_prefactor_sum_B + pow(gsl_vector_get(lengthvect, start_index + edge_index), 2) / (exp(-y.at(edge_index + start_index + type_dimension)) / step_size * alphaP + gsl_vector_get(lengthvect, edge_index + start_index) * te);
			}
			//This is the prefactor for the non deltaij elements (consult page 4 of supplementary info)
			double matrix_prefactor_A = pow(length_sum / te - matrix_prefactor_sum_A, -1);
			double matrix_prefactor_B = pow(length_sum / te - matrix_prefactor_sum_B, -1);
			for (int row_edge_index = 0; row_edge_index != number_of_edges; row_edge_index++){
				for (int column_edge_index = 0; column_edge_index != number_of_edges; column_edge_index++){
					current_matrix_element = gsl_matrix_get(inverse_block_matrix_A, row_edge_index, column_edge_index);
					gsl_matrix_set(inverse_block_matrix_A, row_edge_index, column_edge_index, current_matrix_element + matrix_prefactor_A * gsl_vector_get(lengthvect, start_index + row_edge_index) * gsl_vector_get(lengthvect, start_index + column_edge_index) / (exp(-y.at(row_edge_index + start_index)) / step_size * alphaP + gsl_vector_get(lengthvect, row_edge_index + start_index) * te) / (exp(-y.at(column_edge_index + start_index)) / step_size * alphaP + gsl_vector_get(lengthvect, column_edge_index + start_index) * te));
					current_matrix_element = gsl_matrix_get(inverse_block_matrix_B, row_edge_index, column_edge_index);
					gsl_matrix_set(inverse_block_matrix_B, row_edge_index, column_edge_index, current_matrix_element + matrix_prefactor_B * gsl_vector_get(lengthvect, start_index + row_edge_index) * gsl_vector_get(lengthvect, start_index + column_edge_index) / (exp(-y.at(row_edge_index + start_index + type_dimension)) / step_size * alphaP + gsl_vector_get(lengthvect, row_edge_index + start_index) * te) / (exp(-y.at(column_edge_index + start_index + type_dimension)) / step_size * alphaP + gsl_vector_get(lengthvect, column_edge_index + start_index) * te));
				}
			}
			//Define a vector view to extract the elements of rhs_s12 corresponding to edges of this cell.
			gsl_vector_view this_cell_type_A_view = gsl_vector_subvector(rhs_evolution, start_index, number_of_edges);
			gsl_vector_view this_cell_type_B_view = gsl_vector_subvector(rhs_evolution, start_index + type_dimension, number_of_edges);
			//Then, store this view into an actual gsl vector.
			gsl_vector * this_cell_type_A = &this_cell_type_A_view.vector;
			gsl_vector * this_cell_type_B = &this_cell_type_B_view.vector;
			//Initial a vector to store evolved protein variables of type A and B.
			gsl_vector * evolved_protein_A = gsl_vector_calloc(number_of_edges);
			gsl_vector * evolved_protein_B = gsl_vector_calloc(number_of_edges);
			//Store the inverse block times the rhs_s12 for protein A and B variables for this cell.
			gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_block_matrix_A, this_cell_type_A, 0.0, evolved_protein_A);
			gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_block_matrix_B, this_cell_type_B, 0.0, evolved_protein_B);
			for (int edge_index = 0; edge_index != number_of_edges; edge_index++){
				updated_protein_variables.at(start_index + edge_index) = gsl_vector_get(evolved_protein_A, edge_index);
				updated_protein_variables.at(start_index + edge_index + type_dimension) = gsl_vector_get(evolved_protein_B, edge_index);
			}
			//Free each inverse block matrix and free the evolved protein vectors that have been stored in updated_protein_variables.
			gsl_matrix_free(inverse_block_matrix_A);
			gsl_matrix_free(inverse_block_matrix_B);
			gsl_vector_free(evolved_protein_A);
			gsl_vector_free(evolved_protein_B);
			//Increment the start index by the number of edges in this cell so that you can move to the next cell.
			start_index += number_of_edges;
		}
		cell_index++;
	}

	//store by using gsl_vector_get.
	//Free the memory associated with all of the matrices and vectors used in gsl. Don't allow for memory leaks.
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(V_new);
	gsl_vector_free(rhs_evolution);
	gsl_vector_free(concentration_vector);
	gsl_vector_free(evolution_vectorU);
	gsl_vector_free(lengthvect);
	gsl_vector_free(lengthvect_new);
	gsl_vector_free(lengthtimesconcentration);
	gsl_vector_free(longlengthvector);
	gsl_vector_free(longlengthvector_new);
	gsl_vector_free(loglength);
	gsl_vector_free(Vlogj_vector);
	gsl_vector_free(V_logj_full);
	gsl_vector_free(evolution_vectorV);
	gsl_vector_free(max_conc_term);
	gsl_vector_free(y_this_step);
	gsl_vector_free(g_this_step);
	gsl_vector_free(f_this_step);

	return updated_protein_variables;
}

//ADAPTIVE STEP SIZE CONTROL checks estimated error with error tolerance (and scales step size accordingly).
bool adaptive_step_size_control(const vector<double> &y_error, const double error_tolerance_absolute, double* initial_step_size, void* params, double consistency_order){
	//THIS USES THE ALGORITHM DESCRIBED AT https://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html#Adaptive-Step_002dsize-Control
	//Initialize the flag for a failure in satisying norm tolerance.
	int step_norm_check = 0;
	//Initialize the flag that is output to show if satisying error tolerance.
	bool step_error_check = false;
	//step_error_check == false -> step size within error tolerance.
	//step_error_check == true -> step size needs to be scaled down.

	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	//Define bounds for the step size scale factor (consult zero-stability criterion in https://www.ideals.illinois.edu/bitstream/handle/2142/18708/vsimex_jcm_ideals.pdf?sequence=2)
	double lower_bound_scaling_factor = 0.5818;
	double upper_bound_scaling_factor = 1.1892;

	//Initialize the factor for scaling the step size.
	double scale_factor_double_step = 1.0;
	//Initialize a variable for storing the adjusted step size.
	double adjusted_step_size = *initial_step_size;
	double safety_factor = .9;
	double upper_fractional_tolerance = 1.0;
	double lower_fractional_tolerance = 0.1;

	//Check whether error exceeds tolerance.
	vector<double> excess_error;
	excess_error.resize(dim_proteins);
	//Check whether error is excessively below tolerance.
	vector<double> low_error;
	low_error.resize(dim_proteins);
	//Ratio of error to tolerance
	vector<double> error_ratio;
	error_ratio.resize(dim_proteins);
	//Store the max ratio of the error to the tolerance.
	double max_error_ratio = 0.0;
	//Check for excess error.
	double max_error_check = 0.0;

	//Subtract the tolerated from error from the observed error.
	for (int index = 0; index != dim_proteins; index++){
		excess_error.at(index) = y_error.at(index) - upper_fractional_tolerance * error_tolerance_absolute;
		//Store the ratio of error to tolerance.
		error_ratio.at(index) = y_error.at(index) / error_tolerance_absolute;
	}

	//Store the maximum error ratio.
	max_error_ratio = *max_element(error_ratio.begin(), error_ratio.end());
	//Store the maximum error check.
	max_error_check = *max_element(excess_error.begin(), excess_error.end());

	//If any of the components violated the error tolerance, set a factor for scaling the step size.
	if (max_error_check > 0.0){
		cout << "scale DOWN step due to double step error estimate" << endl;
		scale_factor_double_step = safety_factor * pow(max_error_ratio, -1. / consistency_order);
	}
	
	//Check if all of the errors are below 50% of the tolerance.
	for (int index = 0; index != dim_proteins; index++){
		low_error.at(index) = y_error.at(index) - lower_fractional_tolerance * error_tolerance_absolute;
	}

	//If all of the errors are less than half of the tolerance, define a factor for upscaling the step size.
	max_error_check = *max_element(low_error.begin(), low_error.end());
	if (max_error_check < 0.0){
		scale_factor_double_step = safety_factor * pow(max_error_ratio, -1. / (consistency_order + 1.));
	}

	//If the up scale factor is greater than the upper bound, snap the scale factor to the upper bound.
	if (scale_factor_double_step > upper_bound_scaling_factor){
		scale_factor_double_step = upper_bound_scaling_factor;
	}
	if (scale_factor_double_step < 1.00){
		step_error_check = true;
		cout << "step scaled down" << endl;
	}
	//If the scale factor is lower than the lower bound, snap the scale factor to the lower bound.
	if (scale_factor_double_step < lower_bound_scaling_factor){
		scale_factor_double_step = lower_bound_scaling_factor;
		cout << "Scaled step too abruptly. Step is scaled to:" << *initial_step_size * scale_factor_double_step << endl;
	}

	//Adjust the step size.
	adjusted_step_size = *initial_step_size * scale_factor_double_step;
	*initial_step_size = adjusted_step_size;

	//Return the flag for if the step size was too large.
	return step_error_check;
}

//Update the past information in the solver for both protein and mechanical dynamics. [2] -> [1], [1] -> [0], flush [2]
void update_past_information_solver(void* params){
	tissue.f_step_before[0] = tissue.f_step_before[1];
	tissue.g_step_before[0] = tissue.g_step_before[1];
	tissue.y_step_before[0] = tissue.y_step_before[1];
	tissue.past_step_sizes[0] = tissue.past_step_sizes[1];
	if (evolve_mech_positions == 1){
		tissue.f_mech_step_before[0] = tissue.f_mech_step_before[1];
		tissue.past_mech_state[0] = tissue.past_mech_state[1];
	}
	tissue.f_step_before[1] = tissue.f_step_before[2];
	tissue.g_step_before[1] = tissue.g_step_before[2];
	tissue.y_step_before[1] = tissue.y_step_before[2];
	tissue.past_step_sizes[1] = tissue.past_step_sizes[2];
	if (evolve_mech_positions == 1){
		tissue.f_mech_step_before[1] = tissue.f_mech_step_before[2];
		tissue.past_mech_state[1] = tissue.past_mech_state[2];
	}


	//Flush the temporary part of these matrices [2]
	double init_val = 0.0;
	int dim_protein = (*(combinedparameters*)params).dimension_protein;
	vector<double> vec(dim_protein, init_val);
	tissue.f_step_before[2] = vec;
	tissue.g_step_before[2] = vec;
	tissue.y_step_before[2] = vec;
	tissue.past_step_sizes[2] = init_val;
	if (evolve_mech_positions == 1){
		int dim_mechanical = (*(combinedparameters*)params).dim_mechanical;
		vector<Vect> vec_vect(dim_mechanical, Vect(init_val, init_val));
		tissue.f_mech_step_before[2] = vec_vect;
		tissue.past_mech_state[2] = vec_vect;
	}
}

//Back up information from solver before doing a calculation (just in case you need it).
void back_up_past_information_solver(){
	tissue.f_step_before_save = tissue.f_step_before;
	tissue.g_step_before_save = tissue.g_step_before;
	tissue.y_step_before_save = tissue.y_step_before;
	tissue.past_step_sizes_save = tissue.past_step_sizes;
	if (evolve_mech_positions == 1){
		tissue.f_mech_step_before_save = tissue.f_mech_step_before;
		tissue.past_mech_state_save = tissue.past_mech_state;
	}
}

//Restore past information after you have changed the information in a calculation.
void restore_past_information_solver(){
	tissue.f_step_before = tissue.f_step_before_save;
	tissue.g_step_before = tissue.g_step_before_save;
	tissue.y_step_before = tissue.y_step_before_save;
	tissue.past_step_sizes = tissue.past_step_sizes_save;
	if (evolve_mech_positions == 1){
		tissue.f_mech_step_before = tissue.f_mech_step_before_save;
		tissue.past_mech_state = tissue.past_mech_state_save;
	}
}

//IMEX SOLVER ADAPTIVE evolves the current protein log variables with the input function solver_fcn
bool imex_solver_adaptive(vector<double> &y, double* prot_stepSize, double* mech_stepSize, const double error_abs_tolerance_prot, const double error_abs_tolerance_mech, double* t, void* params, vector<double>(*solver_fcn) (const vector<double> &y, double step_size, void* params, ofstream* errorfile), ofstream* errorfile){
	//Store the size of the system to be solved.
	int dim_proteins = (*(combinedparameters*)params).dimension_protein;
	if (evolve_mech_positions == 1){
		//Check if the mechanical step is less than the protein step.
		if (*mech_stepSize < *prot_stepSize){
			cout << "mech step less than than prot step." << endl;
		}
	}

	//Store the current step size (make sure the total scale factor is not out of range).
	double initial_step_size;
	//Store the initial step size.
	if (evolve_mech_positions == 1){
		vector<double> step_sizes;
		step_sizes.resize(2);
		step_sizes.at(0) = *prot_stepSize; step_sizes.at(1) = *mech_stepSize;
		initial_step_size = *min_element(step_sizes.begin(), step_sizes.end());
	}
	else{
		initial_step_size = *prot_stepSize;
	}
	double current_step_size = initial_step_size;
	double current_step_size_mech = current_step_size;
	double used_step_size;
	//Define the size of a half step.
	double half_step = current_step_size / 2.;
	//Initialize the error check (combine error check for both the protein and mechanical steps).
	bool step_error_check = true;
	bool step_error_check_prot = true;
	bool step_error_check_mech = true;
	//These vectors are used for both half and full steps below as well as error tolerance and estimated error.
	vector<double> error_tolerance_prot;
	vector<double> error_vector_prot;
	vector<double> full_step_solution_prot;
	vector<double> second_half_step_solution_prot;
	vector<double> first_half_step_solution_prot;
	error_tolerance_prot.resize(dim_proteins);
	error_vector_prot.resize(dim_proteins);
	full_step_solution_prot.resize(dim_proteins);
	second_half_step_solution_prot.resize(dim_proteins);
	first_half_step_solution_prot.resize(dim_proteins);

	//dimension of mechanical systems
	int dimension_mechanical = tissue.vertices.size() + 1;
	//These vectors are used for both half and full steps below as well as error tolerance and estimated error.
	vector<Vect> initial_condition_mech;
	vector<Vect> error_tolerance_mech;
	vector<Vect> error_vector_mech;
	vector<Vect> full_step_solution_mech;
	vector<Vect> second_half_step_solution_mech;
	vector<Vect> first_half_step_solution_mech;
	initial_condition_mech.resize(dimension_mechanical);
	error_tolerance_mech.resize(dimension_mechanical);
	error_vector_mech.resize(dimension_mechanical);
	full_step_solution_mech.resize(dimension_mechanical);
	second_half_step_solution_mech.resize(dimension_mechanical);
	first_half_step_solution_mech.resize(dimension_mechanical);

	//Store a variable for storing the information from the first half step.
	combinedparameters half_step_parameters = (*(combinedparameters*)params);
	//The consistency order affects how the step size is scaled with observed error.
	double consistency_order;

	//The lower bound of the scaling factor for the step size (for zero-stability)
	double lower_bound_scaling_factor = 0.5818;

	//Check what solver is being used, and store the order of the method for step size control.
	if (solver_fcn == &imex_solver_PCP_third_order){
		consistency_order = 3.;
	}
	if (solver_fcn == &imex_solver_PCP_second_order){
		consistency_order = 2.;
	}
	if (solver_fcn == &imex_solver_PCP_first_order || solver_fcn == &imex_solver_guillaume_PCP){
		consistency_order = 1.;
	}

	//Store initial condition.
	if (evolve_mech_positions == 1){
		tissue.Store_Initial_Condition(initial_condition_mech);
	}

	//step_error_check == false -> step size within error tolerance.
	//step_error_check == true -> step size needs to be scaled down.
	//Step error check keeps track of whether step size was too large (and thus was scaled down)
	while (step_error_check == true){
		back_up_past_information_solver();

		//Reinitialize half_step_size (upon entering loop)
		half_step = current_step_size / 2.;

		//Take the first half step.
		first_half_step_solution_prot = solver_fcn(y, half_step, (void*)&half_step_parameters, errorfile);
		//Store the calculated results for the first half step for vertices and box size.
		if (evolve_mech_positions == 1){
			tissue.Store_New_Positions(first_half_step_solution_mech);
		}

		//Update the params for the first half step (now including the information from that step).
		update_past_information_solver((void*)&half_step_parameters);
		//Make the results from the first half step the initial condition for the next step.
		if (evolve_mech_positions == 1){
			tissue.Update_Old_Mechanical_Variables(first_half_step_solution_mech);
		}

		//Take the second half step.
		second_half_step_solution_prot = solver_fcn(first_half_step_solution_prot, half_step, (void*)&half_step_parameters, errorfile);
		update_past_information_solver((void*)&half_step_parameters);
		//Store the calculated results for the second half step for vertices and box size.
		if (evolve_mech_positions == 1){
			tissue.Store_New_Positions(second_half_step_solution_mech);
			//restore initial conditions:
			tissue.Update_Old_Mechanical_Variables(initial_condition_mech);
		}

		restore_past_information_solver();

		//Calculate the result using a full step.
		full_step_solution_prot = solver_fcn(y, current_step_size, params, errorfile);
		if (evolve_mech_positions == 1){
			//Store the calculated results for the second half step for vertices and box size.
			tissue.Store_New_Positions(full_step_solution_mech);
			//restore initial conditions:
			tissue.Update_Old_Mechanical_Variables(initial_condition_mech);
		}

		//Estimate error by subtracting two half steps from full step for protein variables.
		error_vector_prot.resize(dim_proteins);
		for (int index = 0; index != dim_proteins; index++){
			error_vector_prot.at(index) = fabs(full_step_solution_prot.at(index) - second_half_step_solution_prot.at(index));
		}

		if (evolve_mech_positions == 1){
			//Estimate error by subtracting two half steps from full step for mechanical variables (vertex positions and box size).
			error_vector_mech.resize(dimension_mechanical);
			for (int index = 0; index != dimension_mechanical; index++){
				error_vector_mech.at(index) = full_step_solution_mech.at(index) - second_half_step_solution_mech.at(index);
				error_vector_mech.at(index) = error_vector_mech.at(index).Abs();
			}
		}

		//Check whether step size is too large for protein variables. (if it is too large, this function will return step_error_check != false)
		step_error_check_prot = adaptive_step_size_control(error_vector_prot, error_abs_tolerance_prot, &current_step_size, params, consistency_order);

		if (evolve_mech_positions == 1){
			//Check whether step size is too large for mechanical variables. (if it is too large, this function will return step_error_check != false)
			step_error_check_mech = tissue.Double_Step_Size_Control(error_vector_mech, error_abs_tolerance_mech, &current_step_size_mech, consistency_order);
		}

		//Combine the error check for both the protein and mechanical evolution.
		if (evolve_mech_positions == 1){
			step_error_check = step_error_check_prot || step_error_check_mech;
		}
		else{
			step_error_check = step_error_check_prot;
		}

		//Note: it is necessary to calculate 2 * half_step because this is the step size used in the evolution above (current_step_size could have been scaled up after evolution).
		used_step_size = 2 * half_step;

		//if the step size was scaled down, reinitialize initial_step_size for reentering while loop.
		if (step_error_check == true){
			*prot_stepSize = current_step_size;
			if (evolve_mech_positions == 1){
				*mech_stepSize = current_step_size_mech;
			}
		}

		//Check whether the protein and mechanical step size need to be scaled up.
		//If either step size is above the max step size, force them to the maximum step size.
		if (current_step_size > *prot_stepSize) *prot_stepSize = current_step_size;
		if (*prot_stepSize > max_step_size) *prot_stepSize = max_step_size;
		if (evolve_mech_positions == 1){
			if (current_step_size_mech > *mech_stepSize) *mech_stepSize = current_step_size_mech;
			if (*mech_stepSize > max_step_size) *mech_stepSize = max_step_size;
		}
	}

	//If the step size was overall scaled more than the lower bound of the scaling factor.
	if (current_step_size / initial_step_size < lower_bound_scaling_factor){
		printf("overall step scale factor too low.");
	}

	//Increment the time by the utilized step size (once the step is successful).
	if (step_error_check == false){
		//It is important to increment the time by two half steps (because the variable current_step_size may have been scaled up).
		*t += used_step_size;
		//only assign these values if the step does not fail.
		for (int variable_index = 0; variable_index != dim_proteins; variable_index++){
			y.at(variable_index) = full_step_solution_prot.at(variable_index);
		}
		//NOTE: mechanical positions are updated outside of this function (within Evolve).
		//Now that you have completed the step successfully update the vectors that store information from past steps.
		update_past_information_solver(params);
		//CONVENTION:
		//before a step is taken, [0] store info from two steps before, [1] stores info from one step before.
		//after a step is taken, [0] store info from two steps before, [1] stores info from one step before, [2] stores information from the step that just happened.
		//then, within the imex_solver_adaptive, [1]->[0], [2]->[1] such that the information is updated for the next step. 
	}

	//temp for debugging: Output the time after this step (just so the user can have some clue).
	cout << "Time after this step:" << *t << endl;
	cout << "Step_size:" << current_step_size << endl;

	return step_error_check;
}


#endif