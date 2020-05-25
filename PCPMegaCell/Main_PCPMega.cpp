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
#include "PCPMegaTissueClass.h"
#include "PCPMegaCellClass.h"
//Additional header file for MegaPCP
#include "MegaCellConstants.h"
//Header file for gsl matrices
#include "gsl/gsl_matrix.h"

using namespace std;

//Defult settings of the constants if the file does not specify	constants given a null or zero value must be specified in the file
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
double alphaM = 1.0;
double TISSUE_ALPHA = 1.0;
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
double additive_tension_for_partners;
double time_for_propagation;
double increment_for_propagation;
double time_for_renormalization;
double increment_for_renormalization;
double kT1;
double alphaP;
double cutline_pcp_initialization;
//Defining cell types.
int standard_cell_type = 0;
int pcp_cell_type = 3;								//!< Type that indicates that a cell has PCP
int ticking_type = 4;								//!< Type that indicates that a cell has a clock ticking for PCP (Poisson)
int intersect_type = 6;								//!< Type that indicates the original front of PCP propagation
int center_cell_type = 5;							//!< Type of the central mega cell.
int outer_cell_type = 7;							//!< Type of outer mega cell.
int checkerboard_type1 = 1;
int checkerboard_type2 = 2;
//Define upper and lower cutoff for defining intermediate protein amounts.
double upper_protein_cutoff = 0.90;
double lower_protein_cutoff = 0.10;
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

//Fix defects: will use this to test hypothesis on half column formation
vector<PCPCell*> fix_defects;
double time_for_defect_fix;
//Count how many lines of cells have differentiated.
int differentiation_count = 1;
int partial_column_insertion = 2;

//Changes by Hayden for PCP
//for random number generation of polarity proteins
gsl_rng* polarity_random_generator = gsl_rng_alloc(gsl_rng_mt19937);

//MEGA CELL CONSTANTS
double factor_increased_tension_at_inner_boundary;	//Note: This allows for increased tension at the inner boundary.
//Attempt to compress packing.
double inner_pressure = 0.25;
double outer_pressure = 0.0;
int is_voronoi_packing;
//Initialize a vector for storing the origin of the inner mega cell.
Vect origin_mega_cell = Vect(0.0, 0.0);

//To convert this Main.cpp to work for Tissue instead of PCPTissue, use RK4_Adaptive.
string INTEGRATORNAME = "IMEX";				//Not currently used within PCPTissue (only for output)

//To change this code to work for Tissue rather than PCPTissue, use Tissue tissue;
PCPMegaTissue tissue;

void main(){
	//ask for the imput file name
	cout << "Enter the file name containing the initial configuration:" << endl;
	getline(cin, InputFileName);

	//Tell the code whether you want to read in protein concentrations from the input file.
	bool includeConcentrations = false;
	//Convert read input data not to take includeConcentrations if you want a Tissue instead of PCPTissue.
	//Read the input data and set any other constants like A0 or cell type that need to be set manually
	tissue.ReadInputData(InputFileName, includeConcentrations);

	//Define the protein and mechanical step size.
	prot_step_size = INITIALSTEPSIZE;
	mech_step_size = INITIALSTEPSIZE;
	//Store the first time point for both renormalization and PCP generation and differentiation.
	time_for_propagation = CURRENTTIME + increment_for_propagation;
	time_for_renormalization = CURRENTTIME + increment_for_renormalization;
	time_for_defect_fix = CURRENTTIME + increment_for_propagation * 9. / 10.;

	//Seed both random number generators
	//Uses seed defined in input file.
	srand(SEED);
	gsl_rng_set(polarity_random_generator, SEED);

	//Declare strings for file names.
	string filename = InputFileName.substr(0, InputFileName.size() - 4);
	string OutputFileName1 = filename + "_Out" + ".txt";
	string OutputFileName2 = filename + "_OutAppend" + ".txt";
	string OutputAfterInitialization = filename + "_initializedpolarityline" + ".txt";
	string OutputForDebug = filename + "_debug" + ".txt";

	//open file for output append
	ofstream datafile;
	datafile.open(OutputFileName2);
	ofstream datafile_debug;
	datafile_debug.open(OutputForDebug);

	//Output the protein concentration after initialization (whether or not they are randomly generated within the code).
	tissue.WriteOutputData4Matlab(OutputAfterInitialization);

	//Evolve the tissue and output the results
	while (CURRENTTIME < TOTALTIME){
		tissue.WriteOutputDataAppend(&datafile);
		//If you want a Tissue instead of a PCPTissue, don't take &prot_step_size as input.
		//Evolve both protein variables and mechanical variables.
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
