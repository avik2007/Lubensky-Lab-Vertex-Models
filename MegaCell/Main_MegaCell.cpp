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

#include "CellClass.h"
#include "EdgeClass.h"
#include "VertexClass.h"
#include "TissueClass.h"
#include "Constants.h"
#include "MegaCellClass.h"
#include "MegaTissueClass.h"
#include "MegaCellConstants.h"

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
double MINAREA = 0.01;							//Note: minimum area for apoptosis in MEGACELL SIMULATION.
double MINLENGTH = 0;
double NEWEDGELENGTH = 0;
int SEED = 2;									//Note: setting seed to 0 or 1 gives unwanted behavior
double K = 0;
double GAMMA = 0.0;
double INITIALSTEPSIZE = 0;
double mech_step_size = INITIALSTEPSIZE;
double TIMEINTERVAL = 0;
double TOTALTIME = 0;
double alphaM = 1.0;
double TISSUE_ALPHA = 1.0;
double TISSUE_EXT_STRESS_X;
double TISSUE_EXT_STRESS_Y;
mt19937_64 generator;
double CURRENTTIME;
double error_tolerance_absolute_mechanical = 1e-3;
double factor_increased_tension_at_inner_boundary = 2.0;	//Note: This allows for increased tension at the inner boundary.

//Attempt to compress packing: exert increased pressure at the inner cell, no pressure at the outer cell.
double inner_pressure = 0.25;
double outer_pressure = 0.0;
int is_voronoi_packing = 0;
int center_cell_type = 5;								//!< Type of the central mega cell.
int outer_cell_type = 7;								//!< Type of outer mega cell.

//Declare global variable for vector (origin of the central mega cell).
Vect origin_mega_cell;

//Solve with RK4 adaptive.
string INTEGRATORNAME = "RK4_Adaptive";

//Create a megatissue object.
MegaTissue tissue;

void main(){
	//ask for the imput file name
	cout << "Enter the file name containing the initial configuration:" << endl;
	getline(cin, InputFileName);

	//Read the input data and set any other constants like A0 or cell type that need to be set manually
	tissue.ReadInputData(InputFileName);

	//Define the protein and mechanical step size.
	mech_step_size = INITIALSTEPSIZE;

	//Seed both random number generators
	//Uses seed defined in input file.
	srand(SEED);

	//Declare strings for output file names.
	string filename = InputFileName.substr(0, InputFileName.size() - 4);
	string OutputFileName1 = filename + "_Out" + ".txt";
	string OutputFileName2 = filename + "_OutAppend" + ".txt";
	string OutputAfterInitialization = filename + "_initializedpolarityline" + ".txt";
	string OutputForDebug = filename + "_debug" + ".txt";

	//Output the protein concentration after initialization (whether or not they are randomly generated within the code).
	tissue.WriteOutputData4Matlab(OutputAfterInitialization);

	//open file for output append
	ofstream datafile;
	datafile.open(OutputFileName2);
	ofstream datafile_debug;
	datafile_debug.open(OutputForDebug);

	//Evolve the tissue and output the results
	while (CURRENTTIME < TOTALTIME){
		tissue.WriteOutputDataAppend(&datafile);
		//If you want a Tissue instead of a PCPTissue, don't take &prot_step_size as input.
		//Evolve both protein variables and mechanical variables.
		tissue.Evolve(TIMEINTERVAL, &mech_step_size, &datafile_debug);
	}
	tissue.WriteOutputDataAppend(&datafile);
	datafile.close();
	datafile_debug.close();

	//Write output data (including virtually all information necessary for describing the state of the system).
	tissue.WriteOutputData4Matlab(OutputFileName1);
}
