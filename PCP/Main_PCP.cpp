// The main function of the vertex model is found here
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector> 
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <time.h>
#include <unordered_map>

//Base class files
#include "TissueClass.h"
#include "Constants.h"
//PCP header files
#include "PCPTissueClass.h"
#include "PCPEdgeClass.h"
#include "PCPConstants.h"

using namespace std;

//Default settings of the constants if the file does not specify	constants given a null or zero value must be specified in the file
//If they are not declared with a default value their will be a unresolved external error
string InputFileName = "NULL";
string OutputFileName = "NULL";
string FileType = "Null";
double LX = 0;
double LY = 0;
double LXnew = 0;
double LYnew = 0;
double MINAREA = 0;
double MINLENGTH = 0;
double NEWEDGELENGTH = 0;
int SEED = 2;										//Note: setting seed to 0 or 1 gives unwanted behavior
double K = 0;
double GAMMA = 0.0;
double INITIALSTEPSIZE = 0;
double prot_step_size = INITIALSTEPSIZE;
double mech_step_size = INITIALSTEPSIZE;
double TIMEINTERVAL = 0;
double TOTALTIME = 0;
double alphaM = 100.0;
double TISSUE_ALPHA;
double TISSUE_EXT_STRESS_X;
double TISSUE_EXT_STRESS_Y;
mt19937_64 generator;
double CURRENTTIME;

//Changes by Hayden for PCP: declare variables needed for PCP dynamics.
//These are read in from the input file.
double epsab;
double epsnaught;
double epsone;
double epstwo;
double jone;
double jtwo;
double te;
double GAMMA0;
double GAMMAC;
double time_for_renormalization;
double increment_for_renormalization;
double kT1;
double alphaP;
double cutline_pcp_initialization = 0.0;
//Defining cell types.
int standard_cell_type = 0;
int pcp_cell_type = 3;								//!< Type that indicates that a cell has PCP
int ticking_type = 4;								//!< Type that indicates that a cell has a clock ticking for PCP (Poisson)
int intersect_type = 6;								//!< Type that indicates the original front of PCP propagation
int checkerboard_type1 = 1;
int checkerboard_type2 = 2;

//Store bool for reinitialize IMEX solver, needs to be external in the case of a T1/Resolve at the end of a call of Evolve.
bool reinitialize_IMEX_solver = true;
//Force a maximum step size (this is not absolutely necessary for low protein tolerances, but is helpful)
double max_step_size;

//Definitions for the integrator
//Using separate evolvers for evolution of vertices and evolution of protein concentrations
//protein label indicates that they are used for protein solver, mechanical solver uses no special label.
double error_tolerance_absolute_mechanical = 1e-3;
//For explicit treatment of protein dynamics, change this to an explicit solver like gsl_odeiv2_step_rkf45
double error_tolerance_absolute_protein = 1e-3;

//Declare solver types for first, second, and third order root finders.
const gsl_multiroot_fdfsolver_type* nonlinear_solver_PCP_first_order;
const gsl_multiroot_fdfsolver_type* nonlinear_solver_PCP_second_order;
const gsl_multiroot_fdfsolver_type* nonlinear_solver_PCP_third_order;

//New constants for PCP - Chiral Simulations
int evolve_mech_positions = 0;
int simulation_type = 0;
double delta_protein_amount = 0.1;
double max_concentration = 0.95;
double dirichlet_noise_param = 1.0;
//Tell the code whether you want to read in protein concentrations from the input file.
bool includeConcentrations = false;

//UNUSED:
double standard_dev_pcpwaitingtime = 0.0;
double mean_pcpwaitingtime = 0.0;
double last_diff_time = 0.0;

//Changes by Hayden for PCP
//for random number generation of polarity proteins
gsl_rng* polarity_random_generator = gsl_rng_alloc(gsl_rng_mt19937);

//To convert this Main.cpp to work for Tissue instead of PCPTissue, use RK4_Adaptive.
string INTEGRATORNAME = "IMEX";				//Not currently used within PCPTissue (only for output)

//To change this code to work for Tissue rather than PCPTissue, use Tissue tissue;
PCPTissue tissue;

void main(){
	//ask for the imput file name
	cout << "Enter the file name containing the initial configuration:" << endl;
	getline(cin, InputFileName);

	//Seed both random number generators
	//Uses seed defined in input file.
	srand(SEED);
	gsl_rng_set(polarity_random_generator, SEED);

	//Convert read input data not to take includeConcentrations if you want a Tissue instead of PCPTissue.
	//Read the input data and set any other constants like A0 or cell type that need to be set manually
	tissue.ReadInputData(InputFileName, includeConcentrations);

	//Define the protein and mechanical step size.
	prot_step_size = INITIALSTEPSIZE;
	mech_step_size = INITIALSTEPSIZE;

	//Store the first time point for both renormalization and PCP generation and differentiation.
	time_for_renormalization = CURRENTTIME + increment_for_renormalization;

	//Declare strings for file names.
	string filename = InputFileName.substr(0, InputFileName.size() - 4);
	string OutputFileName1 = filename + "_Out" + ".txt";
	string OutputFileName2 = filename + "_OutAppend" + ".txt";
	string OutputForDebug = filename + "_debug" + ".txt";

	//open file for output append
	ofstream datafile;
	datafile.open(OutputFileName2);
	ofstream datafile_debug;
	datafile_debug.open(OutputForDebug);

	//Evolve the tissue and output the results
	while (CURRENTTIME < TOTALTIME){
		tissue.WriteOutputDataAppend(&datafile);
		tissue.Evolve(TIMEINTERVAL, &mech_step_size, &prot_step_size, &datafile_debug);
		tissue.WriteOutputData4Matlab(OutputFileName1);
	}
	tissue.WriteOutputDataAppend(&datafile);
	datafile.close();
	datafile_debug.close();

	//Free random number generator for protein concentrations.
	gsl_rng_free(polarity_random_generator);

	//Write output data (including virtually all information necessary for describing the state of the system).
	tissue.WriteOutputData4Matlab(OutputFileName1);
}
