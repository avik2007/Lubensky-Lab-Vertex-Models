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


#include "CellClass_TL.h"
#include "EdgeClass_TL.h"
#include "VertexClass_TL.h"
#include "TissueClass_TL.h"
#include "Constants_TL.h"

using namespace std;

//Defult settings of the constants if the file does not specify	constants given a null or zero value must be specified in the file
//If they are not declared with a default value their will be a unresolved external error
string InputFileName = "NULL";
string OutputFileName = "NULL";
string FileType = "Null";
string FileSpecification = "NULL";
double LX = 0;
double LY = 0;
double LXnew = 0;
double LYnew = 0;
double MINAREA = 0;
double MINLENGTH = 0;
double NEWEDGELENGTH = 0;
int SEED = 2;								//Note: setting seed to 0 or 1 gives unwanted behavior
double GAMMA1 = 0;
double GAMMA2 = 0;
double INITIALSTEPSIZE = 0;
double mech_step_size = INITIALSTEPSIZE;
double TIMEINTERVAL = 0;
double TOTALTIME = 0;
double alphaM = 1.0;
double TISSUE_ALPHA;
double TISSUE_EXT_STRESS_X;
double TISSUE_EXT_STRESS_Y;
mt19937_64 generator;
double CURRENTTIME;
double FORCETIME;
double T1TOTAL;
double APOPTOTAL;
double T1count;
double TOTALMITOSES;

double lnew = sqrt(1250*4/sqrt(3)); //length of equilateral triangle of area A0=1250
double rl_limit_type1 = lnew/3; //set to 0 for no variability of rest lengths
double rl_limit_type2 = lnew/3;
double sc_limit_type1 = 0; //set to 0 for no variability of spring constants
double sc_limit_type2 = 0;
//double v_limit = lnew/5;


double error_tolerance_absolute_mechanical = 1e-3;

//To convert this Main.cpp to work for Tissue instead of PCPTissue, use RK4_Adaptive.
string INTEGRATORNAME = "RK4_Adaptive";

//To change this code to work for Tissue rather than PCPTissue, use Tissue tissue;
Tissue tissue;
chrono::time_point<std::chrono::system_clock> starttime, endtime;

//Uncomment for running on flux
//int main(int argc, char* argv[]){  
	//const char* temp = argv[1];
	//InputFileName = temp;
	//cout << "Reading argv[1] (name) it is:  " << InputFileName << endl;
	//temp = argv[2];
	//FileSpecification = temp;
	//cout << "Specification " << FileSpecification << " to file name noted, boss. You can count on me, boss." << endl;


//Uncomment for running on desktop 	
int main(){
	cout << "Enter the file name containing the initial configuration:" << endl;
	InputFileName = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\TriLattNbC_0600.txt";
	FileSpecification = "ApoptosisTest_noApop";

	
	//Will time the program
	//starttime = std::chrono::system_clock::now();

	tissue.ReadInputData(InputFileName);
	//Seed both random number generators
	//Uses seed defined in input file.
	//srand(SEED);
	
	tissue.AssignRestLengths(rl_limit_type1, rl_limit_type2);
	tissue.AssignSpringConstants(sc_limit_type1, sc_limit_type2);
	//tissue.CalculateA0();
	tissue.CalculateA0Random();

	//Define the protein and mechanical step size.
	mech_step_size = INITIALSTEPSIZE;

	
	string filename = InputFileName.substr(0, InputFileName.size() - 4);
	//string OutputFileName1 = FileSpecification + "_" + filename + "_Out" + ".txt";
	//string OutputFileName1 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\randomseedtest_TriLattNbC_0600_Out.txt";
	//string OutputFileName2 = FileSpecification + "_" + filename + "_OutAppend" + ".txt";
	//string OutputFileName2 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\test_TriLattNbC_60000_OutAppend.txt";
	string OutputForDebug =  FileSpecification + "_" + filename + "_debug" + ".txt"; //FileSpecification + "_" + 
	//string OutputFileNameDensity1 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\test_TriLattNbC_60000_DensityOut.txt";
	string OutputFileNameDensity1 = FileSpecification + "_" + filename + "_DensityOut.txt";
	//string OutputFileNameDensity1errorcheck = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\test_TriLattNbC_60000_DensityOutErrorCheck.txt";
	string OutputFileNameDensity1errorcheck = FileSpecification + "_" + filename + "_DensityOutErrorCheck.txt";
	//string OutputFileNameDensity2 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\test_TriLattNbC_60000_DensityOutAppend.txt";
	//string OutputFileNameDensity2 = FileSpecification + "_" + filename + "_DensityOutAppend.txt";

	

	//open file for output append
	//ofstream datafile;
	//datafile.open(OutputFileName2);
	ofstream datafile_debug;
	datafile_debug.open(OutputForDebug);
	//ofstream densityfile;
	//densityfile.open(OutputFileNameDensity2);

	
	
	
	//Evolve the tissue and output the results
	int t = 0;

	/////// density fluctuation inputs
	int circ_numOut = 1000;
	int rad_numOut = 100;
	bool logdist = true;
	int circ_numAppend = 500;
	int rad_numAppend = 20;
	int n = 200;
	//////////////////////////////

	while (CURRENTTIME < TOTALTIME) {


		//tissue.AreaCalc_Append(circ_numAppend, rad_numAppend, logdist, &densityfile, t, n);
		//tissue.IntersectionTest(t, n);

		t = t + 1;
		//tissue.WriteOutputDataAppend(&datafile);
		//If you want a Tissue instead of a PCPTissue, don't take &prot_step_size as input.
		//Evolve both protein variables and mechanical variables.
		tissue.Evolve(TIMEINTERVAL, &mech_step_size, &datafile_debug);
		//tissue.Mitosis();
		//if (CURRENTTIME == 10) {
		//	tissue.vertexNoise(v_limit);
		//}
		if (tissue.RelaxedPacking(0.000001) == true) {
			cout << "TISSUE BELOW RELAXATION THRESHOLD" << endl;
			break;
		}
		
		
	}
	

	//tissue.WriteOutputDataAppend(&datafile);
	//datafile.close();
	datafile_debug.close();
	//densityfile.close();



	//Write output data (including virtually all information necessary for describing the state of the system).
	//tissue.WriteOutputData4Matlab(OutputFileName1);
	//tissue.CircleDensityTest_Out(circ_numOut, rad_numOut, logdist, vertexCount, OutputFileNameDensity1);
	//Will time the program
	starttime = std::chrono::system_clock::now();
	tissue.AreaCalc(circ_numOut, rad_numOut, logdist, OutputFileNameDensity1, OutputFileNameDensity1errorcheck);
	//tissue.AreaCalcTEST();

	//Print how long it took 
	endtime = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endtime - starttime;
	cout << "ELAPSED TIME: " << elapsed_seconds.count() / 60.0 << "  min\n";
	
	return 0;
}
