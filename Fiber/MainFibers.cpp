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
#include "Constants.h"
#include "FiberConstants.h"
#include "FiberTissueClass.h"
#include "FiberCellClass.h"
#include "FiberEdgeClass.h"
#include "NucleationSite.h"
#include "FiberClass.h"

using namespace std;

//Defult settings of the constants if the file does not specify. Constants given a null or zero value must be specified in the file If they are not declared with a default value their will be a unresolved external error
string InputFileName = "NULL";
string INTEGRATORNAME = "Euiler";
string MYOSINDYNAMICS = "Constant";
double LX = 0;
double LY = 0;
double MINAREA = 0;
double A0 = 0;
double MINLENGTH = 1;
double NEWEDGELENGTH = 5;
int SEED = 5;								//Note: setting seed to 0 or 1 gives unwanted behavior
double K = 5e-6;
double GAMMA = 20;
double INITIALSTEPSIZE = 0.001;
double TOTALTIME = 10.0;
double CURRENTTIME = 0.0;
double TIMEINTERVAL = 1.0;
int PLOTNUMBER = 0;
double MAXAREA = 0;
double TISSUE_ALPHA = 50.0;
double alphaM = 1;
double TISSUE_EXT_STRESS = 0;
double TISSUE_EXT_STRESS_X = 0;
double TISSUE_EXT_STRESS_Y = 0;
double TISSUE_EXT_STRESS_X_MAX = 0;
double TISSUE_EXT_STRESS_Y_MAX = 0;
double LAMBDA = 0;
double error_tolerance_absolute_mechanical = 0.2;
double LXnew = 0;
double LYnew = 0;
double TAU = 1;
double TAU_STRESS = 10;
double A0SKEW = 1;
bool MAKEFIBERS = 0;
double SPRINGCONST = 1;
double DTENSIONCONST = 1;
double ANGLECONST = 0.5;
double FIB_BIRTHRATE = 0;
double FIB_INITPROB = 0;
double STRAINDIST;
double LX0 = 0;
double LY0 = 0;
bool StrainY = true;
double RELAXINTERVAL = 100;
int MITOSISANGLE = 0;
double BETAZ=0;
double BETAGAMMA=1;
double B=0;
int TOTALT1 = 0; 
bool CFO = 1; 
mt19937_64 generator;
uniform_real_distribution<double> distribution;
double MAXGAMMA = 11; 
double MINPERPFORCE = 5; 
double distance_to_next_vert_cutoff = 0.05;			//added by Hayden (in microns); at fixed distance from next vertex, move fiber to next edge
double starting_distance_from_first_vert = 0.0000;	//added by Hayden (in microns); start the fiber (on the new edge), directly on top of vertex

//Relavent parameters (all length units are um and all time units are min) 
double INITIALVEL = 0.1;			// initial fiber velocity in um/min
double SLOWVEL = 0.1;				// fiber velocity after passing a vertex in um/min //CHANGE BY HAYDEN: MADE TWO VELOCITIES THE SAME
double NUCLEATIONRATE = 0.01;		// Nucleation rate per neucleating site in min-1

//CHANGE BY HAYDEN: TEST THE MOST EXTREME CASE
double INITIALDEATHRATE = 0.04;		// deathrate for fibers that have not passed a vertex in min-1
double PASSEDDEATHRATE = 0.2;		// death rate for fibers that have passed a vertex in min-1
//CHANGE BY HAYDEN: TEST THE MOST EXTREME CASE

double fiber_nuc_cutoff_for_death = 0.9861; //added by Hayden
double fiber_nuc_cutoff_for_nucleation = 0.7449; //added by Hayden

//Constants
FiberTissue tissue;

//*****************************************************************************************************
//********************         MAIN                    ************************************************
//*****************************************************************************************************

// BEGINNING OF MAIN FUNCTION FOR RUNNING A BATCH OF FILES (IN RELEASE)
int main(int argc, char* argv[]){  
	cout << "You made it into main congratulations!" << endl;
	const char* temp = argv[1];
	InputFileName = temp;
	cout << "Reading argv[1] (name) it is:  " << InputFileName << endl;

//// BEGINNING OF MAIN FUNCTION FOR RUNNING ONE SPECIFIC FILE AND DEBUGGING
//int main(){
//	cout << "You made it into main congratulations!" << endl;
//	InputFileName = "C:\\Users\\nunley.UMROOT\\Desktop\\Fiber Simulation\\24_11_VMformat.txt";
	
	//***************************************************************************************************************
	//******************************SETUP FILE HANDELING*************************************************************
	//***************************************************************************************************************
	//Set up the output files 
	string OutputFileName1, OutputFileName2, OutputFileName3, OutputFileName4;
	string filename = InputFileName.substr(0, InputFileName.size() - 4);

	//FILE FOR OUTPUTTING APPENDED FILE OF STATE OF SIMULATION AT SEVERAL TIMEPOINTS
	OutputFileName1 = filename + "_out.txt";

	////FILE FOR OUTPUTTING STATISTICS OF FIBER LIFETIMES
	//OutputFileName2 = filename + "_stats.txt";
	////FILE FOR OUTPUTTING FINAL STATE OF SIMULATION
	//OutputFileName3 = filename + "_final.txt";

	ofstream datafile, endfile, otherfile;
	datafile.open(OutputFileName1);

	//endfile.open(OutputFileName2);
	//otherfile.open(OutputFileName3);

	bool readfile;
	readfile = tissue.ReadInputData(InputFileName);
	if (!readfile){
		cout << "Ending because file was not readable" << endl;
		return 0;
	}

	////////***************************INITIALIZE**********************************************************
	SEED =847288;
	mt19937_64 generator(SEED);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	tissue.InitializeFromData();			//removes the boundary 
	tissue.CreateNucleationSites();			//creates fiber verticies from their parent verticies 
	tissue.RecordCellSize();				//record whether a cell is above or below median size for this packing
	tissue.WriteOutputData(&datafile);		//write output data (in appended file).

	////***************************EVOLVE**********************************************************
	//CURRENTTIME = 0; 
	//TOTALTIME = 180;
	//while (CURRENTTIME < TOTALTIME){
	//	cout << endl << "CURRENTTIME  " << CURRENTTIME << "  N fibers  " << tissue.fibers.size() << endl; 
	//  // step size for moving fiber positions
	//	double dt_fib = 0.005;
	//  // time interval (at which to output state of simulation
	//	double tstep_tot_fib = 2.0;
	//  // evolve the fibers
	//	tissue.EvolveFibers(dt_fib, tstep_tot_fib);
	//  // write to appended file (state of simulation)
	//	tissue.WriteOutputData(&datafile);
	//}
	// // write statistics of fiber lifetimes
	//tissue.WriteOutputStatistics(&endfile);
	// // write the final state of simulation to file.
	//tissue.WriteOutputData(&otherfile);

	datafile.close();

	//endfile.close();
	//otherfile.close();


	//***************************END**********************************************************
	cout << "DONE";
	return 0;
}
