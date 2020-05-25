//! \file PCPConstants.h
//! This file contains declaration of global variables/parameters for the simulation.
//! The following variables are global and thus are accessible from any function in the simulation.
//! Note that this file only contains declaration of the global variables. Their values should be set in the Main.cpp file or through file/user input.
//! This only includes the variables needed for PCP dynamics (consult Constants.h for all other constants)

// include guard
#ifndef __PCPCONSTANTS_H_INCLUDED__
#define __PCPCONSTANTS_H_INCLUDED__

#include <string>
#include "gsl/gsl_multiroots.h"			//Use for multi-dimensional root-finding.
#include <random>
#include "VectClass.h"
//Changes by Hayden for PCP: Now I am going to use a dirichlet random number generation from gsl for
//initializing PCP concentrations on the cell edges.
#include "gsl\gsl_rng.h"

using namespace std;

//NOTE: It is necessary to change this to PCPMegaTissue or PCPTissue depending on the type of the simulation.
#define TISSUE_TYPE PCPPoisson

//Changes by Hayden for PCP: These are additional parameters needed to determine PCP dynamics.
//For further information on the meanings of these parameters, consult figure S5(A) in PLOS 2014 Lubensky paper.
extern double epsab;									//!< This modulates the favorability of having proteins of different types on the same side of a given edge
extern double epsnaught;								//!< This modulates the strength of a positive feedback: higher concentration of a protein generates a force that favors a higher concentration of that protein
extern double epsone;									//!< This modulates the favorability of having proteins of different types on opposite sides of a given edge
extern double epstwo;									//!< This modulates the favorability of having proteins of the same type on opposite sides of a given edge
extern double jone;										//!< The larger this number, the larger the tendency to favor equal concentration of the same protein along every edge.
extern double jtwo;										//!< This is the strength of the force that favors the segregration of proteins of different types within the same cell
extern double te;										//!< This is used to change the strength of the entropic redistribution of polarity proteins within a cell.
extern double GAMMA0;									//!< This is the tension in an edge (of a differentiated cell) that has no polarity protein.
extern double GAMMAC;									//!< This keeps track of the constant that is used to modulate edge tension linearly based on polarity protein concentration.
extern double time_for_renormalization;		 			//!< This keeps track of the next time at which the next renormalization occurs.
extern double increment_for_renormalization; 			//!< This is the regular simulation time interval at which the next renormalization occurs.
extern double kT1;										//!< This is a constant used in the assignment of polarity protein to a newly formed edge from resolving a fourfold vertex.
extern double alphaP;									//!< This is used to define the time scale of protein dynamics relative to the time scale of vertex evolution (which has an alpha by convention equal to 1). These alphas are inverse viscosities.
extern double cutline_pcp_initialization;				//!< In the case of PCP (not Mega), this specifies the cut line for initializing cells with PCP.
//Define type structure.
extern int standard_cell_type;							//!< Standard cell types (not any of the following three types)
extern int pcp_cell_type;								//!< Type that indicates that a cell has PCP
extern int checkerboard_type1;							//!< Define the two types that make up the checkerboard.
extern int checkerboard_type2;
extern int ticking_type;								//!< Type that indicates that a cell has a clock ticking for PCP (Poisson)
extern int intersect_type;								//!< Type that indicates the original front of PCP propagation

//Separately controlled protein step size and mechanical step size.
extern double prot_step_size;								//!< Step size for protein integration
extern double mech_step_size;								//!< Step size for mechanical integration

//Parameters for gsl_odeiv2 ODE integrator for PCP evolution (look at www.gnu.org/software/gsl for documentation and\or possible values)
extern double error_tolerance_absolute_protein;						//!< Absolute tolerance for adaptive step-size. See GSL documentation for details.

//Random number generator for PCP initial conditions
extern gsl_rng* polarity_random_generator;							//!< Random number generator used for randomized polarity protein assignment.

//Solver for nonlinear root finder (PCP Evolution)
extern const gsl_multiroot_fdfsolver_type* nonlinear_solver_PCP_first_order;		//!< nonlinear solver for first order solving protein dynamics
extern const gsl_multiroot_fdfsolver_type* nonlinear_solver_PCP_second_order;		//!< nonlinear solver for second order solving protein dynamics
extern const gsl_multiroot_fdfsolver_type* nonlinear_solver_PCP_third_order;		//!< nonlinear solver for third order solving protein dynamics

//Check whether T1 or resolve has led to a topological change in a differentiated cell.
extern bool reinitialize_IMEX_solver;

//Force a maximum step size (this is not absolutely necessary for low protein tolerances, but is helpful)
extern double max_step_size;

//Define a counter for storing how many PCP differentiation events have occurred.
extern int differentiation_count;
extern int partial_column_insertion;

//Define the standard deviation and mean for the pcp waiting time (from start of pcp clock).
extern double mean_pcpwaitingtime;
extern double standard_dev_pcpwaitingtime;

//Store last differentiation time.
extern double last_diff_time;

//New constants for PCP - Chiral Simulations
extern int evolve_mech_positions;
extern int simulation_type;
extern double delta_protein_amount;
extern double dirichlet_noise_param;
extern double max_concentration;						//!< This is used to set the maximum concentration on an edge.
														//!< If max concentration is zero, the code drops the max concentration term from the free energy.

#endif