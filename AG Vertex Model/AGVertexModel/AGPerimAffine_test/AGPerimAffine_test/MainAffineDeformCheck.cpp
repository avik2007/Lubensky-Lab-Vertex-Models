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
#include "AGStiffConstants.h""


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

//Uncomment the following lines to run the the code on the comand line as ./vertexmodel nameofimputfile
//int main(int argc, char* argv[]){
//	cout << "You made it into main congradulations!" << endl;
//	const char* temp = argv[1];
//	InputFileName = temp;
//	cout << "Reading argv[1] (name) it is:  " << InputFileName << endl;

int main()
{
	//Path to the folder that your initial files are in 
	//string PathName = "C:\\Users\\abmondal\\Desktop\\";
	string PathName = "C:\\Users\\abmondal\\Documents\\Lubensky-Lab-Vertex-Models\\MATLAB Data Folder\\";


	cout << "Enter name of imput file: ";
	cin >> InputFileName;
	InputFileName = PathName + InputFileName;
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


	//Set up the output files 
	string filename = InputFileName.substr(0, InputFileName.size() - 4);
	//string OutputFileName1 = filename + "_out.txt";
	//string OutputFileName2 = filename + "_error.txt";
	string OutputFileName1 = filename + "energytracking.txt";
	//ofstream datafile, eqfile, errorfile, digraphfile, hfile;
	ofstream datafile;
	datafile.open(OutputFileName1);
	//errorfile.open(OutputFileName2);

	cout << "I made output files!" << endl;
	//******************************INITIAlIZE***********************************************************************
	//**********************************************************************************************************
	Vect o(0, 0);
	for (list<Vertex*>::iterator it = tissue.vertices.begin(); it != tissue.vertices.end(); it++) {
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

	//cout << LX << " " << LY << endl;
	CURRENTTIME = 0;
	double gammas[10] = { 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01 };
	double energies[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	for (int i = 0; i < 10; i++)
	{
		tissue.DeformTissue(1 + gammas[i], 1 / (1 + gammas[i]));
		energies[i] = tissue.CalculateFreeEnergy() / (LX * LY);
		tissue.DeformTissue(1 / (1 + gammas[i]), 1 + gammas[i]);
		datafile << gammas[i] << " " << energies[i] << endl;
	}
	
	endtime = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endtime - starttime;
	cout << "ELAPSED TIME: " << elapsed_seconds.count() / 60.0 << "  min\n";



	return 0;
}

