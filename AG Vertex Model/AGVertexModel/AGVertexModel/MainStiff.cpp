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
#include "AGConstants.h"
#include "AGPerimConstants.h"
#include "TissueClass.h"
#include "VectClass.h"
#include "AGTissue.h"
#include "AGPerimTissue.h"
#include "AGPerimEdge.h"
#include "AGCell.h"
#include "AGStiffTissue.h"
#include "AGStiffCell.h"
#include "AGStiffConstants.h"
#include "AGStiffEdge.h"


using namespace std;

//Defult settings of the constants if the file does not specify. Constants given a null or zero value must be specified in the file If they are not declared with a default value their will be a unresolved external error
string InputFileName = "NULL";
string INTEGRATORNAME = "RK4_Adaptive";
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
double TISSUE_ALPHA = 50.0;
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
double TAU_STRESS = 10;
double A0SKEW = 1;
double LX0 = 0;
double LY0 = 0;
double T1PERCENT = 0;
int MITOSISANGLE = 0;
int TOTALT1= 0; 
bool APOPTOSISON = true;
bool EQUILIBRIUMMITOSIS = true;
bool MITOSISON = true;
bool T1ON = true;
bool VERBOSEON = true;
bool WRITEMITOSIS = true;
double ANSTRENGTH = 0;
double C4 = 0;
double LAMBDA = 10;
double TCJSTRENGTH = 10;
double THETAAN = 0;
double Kstdratio = 0;
double A0stdratio = 0;
double LAMBDAstdratio = 0;


//This is the generator for random numbers to use call "distribution(generator)" 
mt19937_64 generator;
uniform_real_distribution<double> distribution(0.0, 1.0);
chrono::time_point<std::chrono::system_clock> starttime, endtime;

//Constants
AGStiffTissue tissue;

//*****************************************************************************************************
//********************         MAIN                    ************************************************
//*****************************************************************************************************

//Uncomment the following lines to run the the code on the comand line as ./vertexmodel nameofimputfile
//int main(int argc, char* argv[]){
//	cout << "You made it into main congradulations!" << endl;
//	const char* temp = argv[1];
//	InputFileName = temp;
//	cout << "Reading argv[1] (name) it is:  " << InputFileName << endl;

int main(){

	//Path to the folder that your initial files are in 
	string PathName = "C:\\Users\\abmondal\\Documents\\Lubensky-Lab-Vertex-Models\\MATLAB Data Folder\\StiffnessModelGammaSweep\\"; // CHANGE THIS TO MATCH YOUR FOLDER

	cout << "Enter name of imput file: ";
	cin >> InputFileName;
	InputFileName = PathName + InputFileName;
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


	//Set up the output files 
	string OutputFileName1, OutputFileName2;
	string filename = InputFileName.substr(0, InputFileName.size() - 4);
	OutputFileName1 = filename + "_out.txt";
	string OutputFileName3 = filename + "_error.txt";
	ofstream datafile, eqfile;
	datafile.open(OutputFileName1);
	cout << "I made output files!" << endl;

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


	ofstream errorfile;
	errorfile.open(OutputFileName3);
	errorfile << "Entering equilibration stage.\n";
	cout << "Entering equilibration stage.\n";
	double initialEnergy = tissue.CalculateFreeEnergy();
	tissue.WriteOutputDataAppend(&datafile);
	errorfile << "Updating parameters to vary A0 and/or K. \n";
	double currA0 = 0;
	tissue.A0 = currA0;
	tissue.A0std = tissue.A0stdratio * (LX * LY / tissue.cells.size());
	tissue.constantA0 = false;
	tissue.AssignPreferredAreas(tissue.constantA0, tissue.A0, tissue.A0std);
	tissue.constantK = false;
	
	tissue.constantLAMBDA = false;
	
	tissue.meanK = K;
	tissue.Kstd = tissue.Kstdratio * K;
	tissue.AssignAreaModuli(tissue.constantK, tissue.meanK, tissue.Kstd);
	tissue.meanLAMBDA = LAMBDA; 
	tissue.LAMBDAstd = tissue.LAMBDAstdratio * LAMBDA;
	tissue.AssignPerimeterModuli(tissue.constantLAMBDA, tissue.meanLAMBDA, tissue.LAMBDAstd);
	INTEGRATORNAME = "RK4_Adaptive";
	CURRENTTIME = 0;
	cout << "Initial Polydispersity = " << tissue.calculatePolydispersity() << endl;
	
	cout << "CURRENTTIME" << CURRENTTIME << endl;
	tissue.equilibriumtolerance = 1e-7;
	while (!(tissue.IsEquilibrated()))
	{
		cout << "CURRENTTIME" << CURRENTTIME << endl;
		tissue.Evolve(TIMEINTERVAL, &stepSize, &errorfile);
		cout << "exited evolve" << endl;
		tissue.WriteOutputDataAppend(&datafile);
		cout << "CURRENTTIME" << CURRENTTIME << endl;
	}
	errorfile << "Left evolve loop.\n";
	errorfile << "Reporting max and min A0 and K's now\n" << endl;
	vector<Vertex*> cwverts;
	double area, stifffunc, pm;
	double maxa = 0;
	double mina = 1000000000000;
	double maxs = 0;
	double mins = 1000000000000;
	double maxp = 0;
	double minp = 1000000000000;
	double areaavg = 0;
	for (list<Cell*>::iterator cit = tissue.cells.begin(); cit != tissue.cells.end(); cit++)
	{
		AGStiffCell* agit = dynamic_cast<AGStiffCell*>(*cit);
		cwverts = (*(*cit)).SavedInfoForEvolve_CWVerts;							//Get the cw vertices in the cell from the preassigned list 
		area = (*(*cit)).GetArea(cwverts);
		areaavg += area / tissue.cells.size();
		if (maxa < area)
		{
			maxa = area;
		}
		if (mina > area)
		{
			mina = area;
		}
		stifffunc = (*agit).KAM;
		if (stifffunc < mins)
		{
			mins = stifffunc;
		}
		if (stifffunc > maxs)
		{
			maxs = stifffunc;
		}
		pm = (*agit).PM;
		if (pm < minp)
		{
			minp = pm;
		}
		if (pm > maxp)
		{
			maxp = pm;
		}

	}


	double finalEnergy = tissue.CalculateFreeEnergy();
	tissue.WriteOutputDataAppend(&datafile);
	datafile << "------------------------------------------------------------\n \n";
	string OutputFileName4 = filename + "_EQ.txt";
	eqfile.open(OutputFileName4);
	tissue.WriteOutputDataAppend(&eqfile);
	cout << "Final Polydispersity = " << tissue.calculatePolydispersity() << endl;



	errorfile << "Stiffness ratio = " << maxs / mins << endl;
	errorfile << "Max K = " << maxs << endl;
	errorfile << "Min K = " << mins << endl;
	errorfile << "Area ratio = " << maxa / mina << endl;
	errorfile << "Max area = " << maxa << endl;
	errorfile << "Min area = " << mina << endl;
	errorfile << "Max lambda = " << maxp << endl;
	errorfile << "Min lambda = " << minp << endl;
	errorfile << "lambda/kabar ratio = " << tissue.meanLAMBDA / (tissue.meanK * areaavg);


	
	///////******************************End***********************************************************************
	///////**********************************************************************************************************
	endtime = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endtime - starttime;
	cout << "ELAPSED TIME: " << elapsed_seconds.count() / 60.0 << "  min\n";
	return 0;
}

