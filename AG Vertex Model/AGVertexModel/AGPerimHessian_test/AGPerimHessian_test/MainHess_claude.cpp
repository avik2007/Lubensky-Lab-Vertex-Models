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
#include <Eigen/Dense>
#include <complex> 
#include <iomanip>

#include "CellClass.h"
#include "EdgeClass.h"
#include "VertexClass.h"
#include "TissueClass.h"
#include "Constants.h"
#include "AGConstants.h"
#include "AGPerimConstants.h"
#include "TissueClass.h"
#include "VectClass.h"
#include "AGTissue.h"
#include "AGPerimTissue.h"
#include "AGPerimEdge.h"
#include "AGCell.h"


using namespace Eigen;
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
int energyListSize = 2;
double equilibriumtolerance = 0.00000000001;

//This is the generator for random numbers to use call "distribution(generator)" 
mt19937_64 generator;
uniform_real_distribution<double> distribution(0.0, 1.0);
chrono::time_point<std::chrono::system_clock> starttime, endtime;

//Constants
AGPerimTissue tissue;


//*****************************************************************************************************
//********************         MAIN                    ************************************************
//*****************************************************************************************************

// MainHess_claude.cpp -- a copy of MainHessianTest.cpp with minimal changes, written to run
// the Hessian shear modulus as a SLURM array job on Great Lakes. Changes from the original:
//   1. Input comes from the command line instead of cin + a Windows path:
//          ./agperim_hess_claude <input EQ file> <output folder>
//      (this is the original's own commented-out argv block, re-enabled, plus an output folder
//      so nothing is written next to the input data).
//   2. NO evolution: the shear modulus is computed on the EQ tissue exactly as it is read in.
//      The original re-equilibrated to a 1e-17 tolerance (with T1s on) first; that is removed so
//      the Hessian is taken at the same state the squeeze benchmark started from.
//   3. Debugging output removed: the _out.txt / _EQ.txt / _H.txt / _P.txt dumps, the Hessian
//      symmetry printout, and the CalculateShearModulusForPythonCalculation call (Python path only).
//   4. Bug fix: Smd previously called the single-precision function a second time; it now calls
//      the double-precision calculateShearModulusWithEigenDecompositiond.
//   5. Each cell's clockwise vertex list (SavedInfoForEvolve_CWVerts) is set before the modulus
//      is computed, as Evolve's first step would have done (needed because of change 2).
//   6. Output: <output folder>/<input basename>_SM.txt with
//          line 1 = double-precision shear modulus (Smd)   <- what calculateShearModulusfromSMfiles.m reads
//          line 2 = single-precision shear modulus (Smf)
//      Both are the library's normalized value Sm = G * Nc / (K (LX LY)^2); a value of exactly 0
//      means the library's floppy-state early return. The library's own diagnostic text still goes
//      to <output folder>/<input basename>_error.txt, as in the original.

int main(int argc, char* argv[])
{
	if (argc < 3) {
		cout << "Usage: " << argv[0] << " <input EQ file> <output folder>" << endl;
		return 1;
	}
	InputFileName = argv[1];
	string OutputFolder = argv[2];
	cout << "Reading argv[1] (name) it is:  " << InputFileName << endl;
	starttime = std::chrono::system_clock::now();

	//Read the input file
	bool readfile;
	readfile = tissue.ReadInputData(InputFileName);
	if (!readfile) {
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
	string OutputFileNameSM = filename + "_SM.txt";

	ofstream errorfile, SMfile;
	errorfile.open(OutputFileName2);

	//******************************INITIAlIZE***********************************************************************
	//**********************************************************************************************************
	Vect o(0, 0);
	for (list<Vertex*>::iterator it = tissue.vertices.begin(); it != tissue.vertices.end(); it++) {
		Vertex* v = &(*(*it));
		(*v).SavedInfoForEvolve_MechForce = o; //Set the initial forces on vertices to 0
	}

	// Per-cell set-up that the first step of Evolve normally performs (AGPerimTissue.cpp, Evolve):
	// store each cell's clockwise vertex list. The Hessian and shear-modulus functions read
	// SavedInfoForEvolve_CWVerts, so without evolution it must be filled here, otherwise they
	// work on an empty/stale list and return a wrong modulus (found in the 2026-10-15 smoke test).
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++) {
		(*cit)->SavedInfoForEvolve_CWVerts = (*cit)->GetClockWiseVertices();
	}

	//******************************SHEAR MODULUS (no evolution: the EQ tissue as read in)**************************
	double Smf, Smd;
	Smf = tissue.calculateShearModulusWithEigenDecomposition(&errorfile);
	Smd = tissue.calculateShearModulusWithEigenDecompositiond(&errorfile);
	SMfile.open(OutputFileNameSM);
	SMfile << std::setprecision(20) << Smd << endl;
	SMfile << std::setprecision(20) << Smf << endl;
	SMfile.close();

	endtime = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endtime - starttime;
	cout << "Smd = " << std::setprecision(12) << Smd << "  Smf = " << Smf << endl;
	cout << "ELAPSED TIME: " << elapsed_seconds.count() / 60.0 << "  min\n";

	return 0;
}
