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
#include <numeric>
#include <map>
#include <time.h>
#include <unordered_map>
#include <chrono>
#include <random>
#include <ctime>
#include <omp.h>
#include <typeinfo>
#include <cassert> 
#include <stddef.h>

#include "CellClass.h"
#include "EdgeClass.h"
#include "VertexClass.h"
#include "TissueClass.h"
#include "AGCell.h"
#include "AGTissue.h"
#include "Constants.h"
#include "AGConstants.h"
#include "AGPerimConstants.h"
#include "AGPerimTissue.h"



using namespace std;


//Defult settings of the constants if the file does not specify. Constants given a null or zero value must be specified in the file If they are not declared with a default value their will be a unresolved external error
string InputFileName = "NULL";
string INTEGRATORNAME = "RK4_Adaptive"; // // "Euiler";
double LX = 0;
double LY = 0;
double MINAREA = 0;
double A0 = 1;
double MINLENGTH = 0.000001; // for agperim vertex, this is appropriate, leave it as such
double NEWEDGELENGTH = 0.0001;
int SEED = 5;								//Note: setting seed to 0 or 1 gives unwanted behavior
double K = 5e-6;
double GAMMA = 10;
double INITIALSTEPSIZE = 0.001;
double TOTALTIME = 10.0;
double CURRENTTIME = 0.0;
double TIMEINTERVAL = 1.0;
double TISSUE_ALPHA = 500.0;
double alphaM = 1;
double TISSUE_EXT_STRESS = 0;
double TISSUE_EXT_STRESS_X = 0;
double TISSUE_EXT_STRESS_Y = 0;
double TISSUE_EXT_STRESS_X_MAX = 0;
double TISSUE_EXT_STRESS_Y_MAX = 0;
double error_tolerance_absolute_mechanical = 0.2;
double LXnew = 0;
double LYnew = 0;
double TAU = 1;
double TAU_STRESS = 5;
double A0SKEW = 1;
double LX0 = 0;
double LY0 = 0;
double T1PERCENT = 0;
int MITOSISANGLE = 0;
int TOTALT1 = 0;
double THETAAN = 0;
double ANSTRENGTH = 0;
double TCJSTRENGTH = 0;
bool MITOSISON = false;
bool T1ON = true;
bool APOPTOSISON = true;
bool WRITEMITOSIS = false;
bool VERBOSEON = true;
bool EQUILIBRIUMMITOSIS = false;
double P0 = 0;
double LAMBDA = 0.1;
double DEFORM_NUM = 0.001;
double C4 = 1;
int energyListSize = 50;
double equilibriumtolerance = 0.0000000000001;

//This is the generator for random numbers to use call "distribution(generator)" 
mt19937_64 generator;
uniform_real_distribution<double> distribution(0.0, 1.0);
chrono::time_point<std::chrono::system_clock> starttime, endtime;

//Constants
AGPerimTissue tissue;


//*****************************************************************************************************
//********************         MAIN                    ************************************************
//*****************************************************************************************************

// MainAGPerim_claude.cpp -- a copy of MainAGPerim.cpp with minimal changes, written (2026-10-15) to
// re-equilibrate the quartic (C4 != 1) data sets as a SLURM array job after the C4 shadowing fix in
// AGTissue.h. Changes from the original, all taken from the user's own Great Lakes version
// (CompressionInducedFluidizationData/CIF_CodeScriptsToRunInGreatLakes/MainAGPerim.cpp):
//   1. Input from the command line instead of cin + a Windows path:
//          ./agperim_equilibration_claude <raw input file> <output folder>
//      (the output folder keeps the run from writing next to the input data).
//   2. No per-time-interval _out.txt dump (the Great Lakes version also has it commented out).
//   3. The _EQ.txt file is written only if the tissue equilibrated before TOTALTIME, with the
//      "initialEnergy / finalEnergy" footer the Great Lakes version writes (the original repo
//      version always wrote _EQ.txt and had no footer).
// The evolution itself (the while loop calling tissue.Evolve) is unchanged.

int main(int argc, char* argv[])
{
	if (argc < 3) {
		cout << "Usage: " << argv[0] << " <raw input file> <output folder>" << endl;
		return 1;
	}
	InputFileName = argv[1];
	string OutputFolder = argv[2];
	cout << "Reading argv[1] (name) it is:  " << InputFileName << endl;
	starttime = std::chrono::system_clock::now();

	//Read the input file
	bool readfile;
	readfile = tissue.ReadInputData(InputFileName);
	if (!readfile){
		cout << "Ending because file was not readable" << endl;
		return 0;

	}
	//Seed the random number generator
	mt19937_64 generator(SEED);
	uniform_real_distribution<double> distribution(0.0, 1.0);
	

	//Set up the output files: <output folder>/<input file name without folder and without ".txt">
	string basename = InputFileName.substr(InputFileName.find_last_of('/') + 1);
	basename = basename.substr(0, basename.size() - 4);
	string filename = OutputFolder + "/" + basename;
	string OutputFileName2 = filename + "_error.txt";
	
	ofstream eqfile, errorfile;
	errorfile.open(OutputFileName2);
	
	//******************************INITIAlIZE***********************************************************************
	//**********************************************************************************************************
	Vect o(0, 0);
	for (list<Vertex*>::iterator it = tissue.vertices.begin(); it != tissue.vertices.end(); it++){
		Vertex* v = &(*(*it));
		(*v).SavedInfoForEvolve_MechForce = o; //Set the initial forces on vertices to 0
	}

	double stepSize = INITIALSTEPSIZE; 
	//***********************************************************************************************
	//***************************EVOLVE**********************************************************
	//*********************************************************************************************
	cout << endl << "*******************************************************************************" << endl;
	cout << "Evolving with step size" << INITIALSTEPSIZE << endl;
	cout << "*******************************************************************************" << endl;

	CURRENTTIME = 0;
	double initialEnergy = tissue.CalculateFreeEnergy();
	errorfile << initialEnergy << endl;
	
	// above time variables allow us to pick how often we write to the appendfile
	while ((CURRENTTIME < TOTALTIME) && !(tissue.equilibrated))
	{
		cout << "CURRENTTIME " << CURRENTTIME << endl;
		tissue.Evolve(TIMEINTERVAL, &stepSize, &errorfile);
		// equilibrated gets updated in Evolve
	}
	double finalEnergy = tissue.CalculateFreeEnergy();
	if (tissue.IsEquilibrated())
	{
		string OutputFileName3 = filename + "_EQ.txt";
		eqfile.open(OutputFileName3);
		tissue.WriteOutputDataAppend(&eqfile);
		eqfile << "initialEnergy: " << initialEnergy << "\t finalEnergy: " << finalEnergy << endl;
		eqfile << "number of cells: " << tissue.cells.size() << endl;
		eqfile << "This tissue equilibrated prior to exiting evolve." << endl;
	}
	else
	{
		cout << "NOT EQUILIBRATED by TOTALTIME = " << TOTALTIME << "; no _EQ.txt written" << endl;
	}

	///////******************************End***********************************************************************
	///////**********************************************************************************************************
	endtime = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endtime - starttime;
	cout << "ELAPSED TIME: " << elapsed_seconds.count() / 60.0 << "  min\n";
	
	return 0;
}
