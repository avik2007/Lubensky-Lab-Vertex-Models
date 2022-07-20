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

#include "CellClassClone.h"
#include "EdgeClassClone.h"
#include "VertexClassClone.h"
#include "TissueClassClone.h"
#include "ConstantsClone.h"

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
double K = 0;
double GAMMA = 0.0;
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
double lnew;
double lold = pow((2 * 1250 ) / (3 * sqrt(3)), .5);
double A0old;
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
	//cout << "Specification " << FileSpecification << " to file name noted, boss. You can count on me." << endl;


//Uncomment for running on desktop 	
int main(){
	cout << "Enter the file name containing the initial configuration:" << endl;
	InputFileName = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\NbC_0500_rng_101.txt";
	FileSpecification = "TEMP";
	
	//Will time the program
	starttime = std::chrono::system_clock::now();
	
	tissue.ReadInputData(InputFileName);
	//tissue.OriginalA0Input();

	//Define the protein and mechanical step size.
	mech_step_size = INITIALSTEPSIZE;

	//Seed both random number generators
	//Uses seed defined in input file.
	srand(SEED);
	
	string filename = InputFileName.substr(0, InputFileName.size() - 4);
	//string OutputFileName1 = FileSpecification + "_" + filename + "_Out" + ".txt";
	//string OutputFileName1 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\NewLengthScale_DTriNbC_25000_rng_1_Out.txt";
	//string OutputFileName2 = FileSpecification + "_" + filename + "_OutAppend" + ".txt";
	//string OutputFileName2 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\NbC_0500_rng_101_OutAppend.txt";
	string OutputForDebug =  FileSpecification + "_" + filename + "_debug" + ".txt"; //FileSpecification + "_" + 
	//string OutputfileNameCLONE = FileSpecification + "_" + filename + "_CloneInfoAppend.txt";
	//string OutputfileNameCLONE = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\NbC_0500_rng_101_CloneInfoAppend.txt";
	string OutputFileNameDensity1 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\test_NbC_0500_rng_1_DensityOut.txt";
	//string OutputFileNameDensity1 = FileSpecification + "_" + filename + "_DensityOut.txt";
	string OutputFileNameDensity2 = "C:\\Users\\jakobths\\source\\repos\\VertexModel\\VertexModel\\test_NbC_0500_rng_1_DensityOutAppend.txt";
	//string OutputFileNameDensity2 = FileSpecification + "_" + filename + "_DensityOutAppend.txt";

	

	//open file for output append
	//ofstream datafile;
	//datafile.open(OutputFileName2);
	ofstream datafile_debug;
	datafile_debug.open(OutputForDebug);
	//ofstream clonefile;
	//clonefile.open(OutputfileNameCLONE);
	ofstream densityfile;
	densityfile.open(OutputFileNameDensity2);
	
	
	
	//Evolve the tissue and output the results
	int t = 0;

	/////// density fluctuation inputs
	int circ_numAppend = 100;
	int rad_numAppend = 10;
	int circ_numOut = 100;
	int rad_numOut = 10;
	int n = 1;
	bool logdist = true;
	bool vertexCount = false;
	//////////////////////////////

	while (CURRENTTIME < TOTALTIME) {

		if (t == 0) {
			tissue.NotClone(); //Tells all current cells that they are not clones
			tissue.MarkClones(); //Randomly assigns cells as clones
		}

		//////////// Density Fluctuation Contribution
		tissue.CircleDensityTest_OutAppend(circ_numAppend, rad_numAppend, t, n, logdist, vertexCount, &densityfile);
		////////////////////////////////////////

		
		t = t + 1;
		//tissue.WriteOutputDataAppend(&datafile);
		//tissue.WriteOutputCloneInfoAppend(&clonefile);
		//If you want a Tissue instead of a PCPTissue, don't take &prot_step_size as input.
		//Evolve both protein variables and mechanical variables.
		tissue.SetNewA0();
		tissue.Evolve(TIMEINTERVAL, &mech_step_size, &datafile_debug);
		//tissue.Mitosis();
		
		//stops the packing if it exceeds a given size
		//if (tissue.CellLimit() == true) {
		//	break;
		//}
		

	}
	
	cout << "Total Amount of T1 Exhanges: " << T1TOTAL << endl;
	cout << "Total Amount of Mitoses: " << TOTALMITOSES << endl;
	cout << "Total Amount of Apoptoses: " << APOPTOTAL << endl;

	//tissue.WriteOutputDataAppend(&datafile);
	//tissue.WriteOutputCloneInfoAppend(&clonefile);
	//datafile.close();
	datafile_debug.close();
	//clonefile.close();
	densityfile.close();


	//Write output data (including virtually all information necessary for describing the state of the system).
	//tissue.WriteOutputData4Matlab(OutputFileName1); //this does not yet include the Clone stuff. unneeededed for now
	tissue.CircleDensityTest_Out(circ_numOut, rad_numOut, logdist, vertexCount, OutputFileNameDensity1);


	//Print how long it took 
	endtime = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endtime - starttime;
	cout << "ELAPSED TIME: " << elapsed_seconds.count() / 60.0 << "  min\n";
	
	return 0;
}
