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
#include <csignal>	// signal(), sig_atomic_t: save the state when SLURM signals a time limit
#include <iomanip>	// setprecision for the chunk ledger

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
//   4. (2026-10-16) If the tissue did NOT equilibrate by TOTALTIME, its final state is written to
//      <name>_NE.txt, following the user's Great Lakes MainAGDebug.cpp (same footer and message).
//      The _NE.txt file is a complete input file (WriteOutputDataAppend writes all globals, cells
//      and vertices), so it can be fed straight back in as <input file> to continue the
//      equilibration for another TOTALTIME from where it stopped. Without this, an unequilibrated
//      run's state was lost when the program exited.
//   5. (2026-10-17) Chunked, restartable equilibration (protocol in
//      vm_claude_plans/2026-10-17_chunked-equilibration-protocol.md):
//      a. Optional 3rd argument <chunk>: run this many time units, overriding the file's TOTALTIME
//         (the raw inputs carry TOTALTIME up to 1e5, longer than a 72 h SLURM job can reach).
//      b. SIGUSR1 / SIGTERM (sent by SLURM before the time limit, or by scancel) set a flag. The
//         loop stops after the current Evolve(TIMEINTERVAL) call and the state is saved to _NE.txt,
//         so a timeout never loses the run.
//      c. Cumulative time: the CURRENTTIME written to _EQ/_NE.txt is the input file's CURRENTTIME
//         (0 for a raw input) plus the time run in this chunk. The evolution itself still starts
//         each chunk at CURRENTTIME = 0, exactly as before.
//      d. Stable names: any trailing "_NE" is stripped from the input name, so continuing from
//         X_ELS_3_NE.txt writes X_ELS_3_NE.txt / X_ELS_3_EQ.txt again (not X_ELS_3_NE_NE.txt).
//         _error.txt is appended to, with a header line per chunk.
//      e. One line per chunk is appended to <name>_chunks.tsv (the tissue's history).
// The evolution itself (the while loop calling tissue.Evolve) is unchanged.

// Set by the signal handler; checked by the evolve loop after every Evolve(TIMEINTERVAL) call.
// sig_atomic_t + volatile is the only kind of variable a signal handler may safely write.
volatile sig_atomic_t stopRequested = 0;
void requestStop(int) { stopRequested = 1; }

int main(int argc, char* argv[])
{
	if (argc < 3) {
		cout << "Usage: " << argv[0] << " <input file (raw or _NE.txt)> <output folder> [chunk length]" << endl;
		return 1;
	}
	InputFileName = argv[1];
	string OutputFolder = argv[2];
	cout << "Reading argv[1] (name) it is:  " << InputFileName << endl;
	starttime = std::chrono::system_clock::now();
	// Install the stop handlers before the (long) evolution starts.
	signal(SIGUSR1, requestStop);
	signal(SIGTERM, requestStop);

	//Read the input file
	bool readfile;
	readfile = tissue.ReadInputData(InputFileName);
	if (!readfile){
		cout << "Ending because file was not readable" << endl;
		return 0;

	}
	// Chunk length (5a) and the tissue's cumulative time at the start of this chunk (5c).
	if (argc >= 4) { TOTALTIME = atof(argv[3]); }
	double startCumulativeTime = CURRENTTIME;
	cout << "Chunk length " << TOTALTIME << ", cumulative time at start " << startCumulativeTime << endl;
	//Seed the random number generator
	mt19937_64 generator(SEED);
	uniform_real_distribution<double> distribution(0.0, 1.0);
	

	//Set up the output files: <output folder>/<input file name without folder and without ".txt">
	string basename = InputFileName.substr(InputFileName.find_last_of('/') + 1);
	basename = basename.substr(0, basename.size() - 4);
	// 5d: strip any trailing "_NE" so all chunks of one tissue share one set of file names.
	while (basename.size() > 3 && basename.compare(basename.size() - 3, 3, "_NE") == 0)
		basename = basename.substr(0, basename.size() - 3);
	string filename = OutputFolder + "/" + basename;
	string OutputFileName2 = filename + "_error.txt";

	ofstream eqfile, errorfile;
	errorfile.open(OutputFileName2, ios::app);	// 5d: append, one header line per chunk
	errorfile << "##### CHUNK from " << InputFileName << " at cumulative time " << startCumulativeTime
		<< ", chunk length " << TOTALTIME << endl;
	
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
	while ((CURRENTTIME < TOTALTIME) && !(tissue.equilibrated) && !stopRequested)	// 5b: also stop on signal
	{
		cout << "CURRENTTIME " << CURRENTTIME << endl;
		tissue.Evolve(TIMEINTERVAL, &stepSize, &errorfile);
		// equilibrated gets updated in Evolve
	}
	double finalEnergy = tissue.CalculateFreeEnergy();
	double chunkTime = CURRENTTIME;
	bool isEQ = tissue.IsEquilibrated();
	// Final slope of the equilibrium checker (same quantity IsEquilibrated tests), for the ledger.
	double finalSlope = (tissue.energyList.size() > 1) ? getLinearRegressionSlope(tissue.timeList, tissue.energyList) : 0;
	// 5c: write the cumulative time into the saved state.
	CURRENTTIME = startCumulativeTime + chunkTime;
	if (isEQ)
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
		// Not equilibrated: save the final state as <name>_NE.txt (as in the user's MainAGDebug.cpp)
		// so a follow-up run can restart from it instead of from the raw input.
		cout << "NOT EQUILIBRATED by TOTALTIME = " << TOTALTIME << "; final state written to _NE.txt" << endl;
		string OutputFileName3 = filename + "_NE.txt";
		eqfile.open(OutputFileName3);
		tissue.WriteOutputDataAppend(&eqfile);
		eqfile << "initialEnergy: " << initialEnergy << "\t finalEnergy: " << finalEnergy << endl;
		eqfile << "number of cells: " << tissue.cells.size() << endl;
		eqfile << "This tissue DID not equilibrate prior to exiting solve. This the final recorded state." << endl;
		if (stopRequested) { eqfile << "Stopped by signal (SLURM time limit or scancel) before the end of the chunk." << endl; }
	}
	// 6 (2026-10-31): flag a run whose integrator fell back to fixed-step RK4 (AGPerimTissue::Evolve backtrack).
	// The fallback is permanent within a process, so the rest of this chunk ran at a step that could only shrink.
	bool rk4Fallback = (INTEGRATORNAME != "RK4_Adaptive");
	if (rk4Fallback)
	{
		eqfile << "WARNING: RK4_Adaptive switched to fixed-step RK4 during this run (final step size " << stepSize
			<< "). Equilibration may have been throttled and an EQ verdict may be false; consider re-running from this file." << endl;
		cout << "WARNING: RK4_Adaptive switched to fixed-step RK4 during this run (final step size " << stepSize << ")" << endl;
	}
	eqfile.close();

	// 5e: one ledger line per chunk. Header is written only when the ledger is new.
	string ledgerName = filename + "_chunks.tsv";
	bool ledgerIsNew = !ifstream(ledgerName).good();
	ofstream ledger(ledgerName, ios::app);
	if (ledgerIsNew) { ledger << "input\tstart_t\tend_t\tinitialE\tfinalE\tfinal_slope\toutcome\twall_min\trk4_fallback\tfinal_stepsize" << endl; }
	std::chrono::duration<double> chunk_seconds = std::chrono::system_clock::now() - starttime;
	ledger << std::setprecision(17) << InputFileName << "\t" << startCumulativeTime << "\t" << CURRENTTIME << "\t"
		<< initialEnergy << "\t" << finalEnergy << "\t" << finalSlope << "\t"
		<< (isEQ ? "EQ" : (stopRequested ? "SIGNAL" : "NE")) << "\t" << chunk_seconds.count() / 60.0 << "\t"
		<< (rk4Fallback ? "WARNING_RK4_FALLBACK" : "none") << "\t" << stepSize << endl;

	///////******************************End***********************************************************************
	///////**********************************************************************************************************
	endtime = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = endtime - starttime;
	cout << "ELAPSED TIME: " << elapsed_seconds.count() / 60.0 << "  min\n";
	
	return 0;
}
