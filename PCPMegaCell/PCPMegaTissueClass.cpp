/**	\file	PCPMegaTissueClass.cpp
*	Implementation of PCPMegaTissue class.
*	Contains function definitions associated with PCPMegaTissue objects.
*	Depends on base class Tissue and on subclass PCPTissue.
*/

#include <cstdlib>
#include <list>
#include <vector> 
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <ctime>
#include <random>
#include <iostream>
#include <cassert>					//for error checking/assertion
#include <iomanip>					//for changing output precision

//Base class files
#include "Constants.h"
#include "MegaCellConstants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
//PCP subclass files
#include "PCPEdgeClass.h"
#include "PCPTissueClass.h"
#include "PCPConstants.h"
//PCP Mega subclass files
#include "PCPMegaCellClass.h"
#include "PCPMegaTissueClass.h"

using namespace std;

//Initialize PCP concentrations within a line of cells in the tissue around the central MEGACELL.
void PCPMegaTissue::AddDifferentiatedCellColumn(){
	list<Edge*>::iterator edgeIter;
	//Vector for first cells to get PCP.
	vector<Cell*> cell_along_boundary;
	//Find cells that are of type 0, that border the central cell (of type 5).
	for (edgeIter = edges.begin(); edgeIter != edges.end(); edgeIter++){
		if ((*(*(*edgeIter)).cell1).type == center_cell_type && (*(*(*edgeIter)).cell2).type == standard_cell_type) {
			cell_along_boundary.push_back((*(*edgeIter)).cell2); 	//push back a pointer to the edge if it has a yflag and is at the boundary of propagation.
		}
		if ((*(*(*edgeIter)).cell1).type == standard_cell_type && (*(*(*edgeIter)).cell2).type == center_cell_type) {
			cell_along_boundary.push_back((*(*edgeIter)).cell1); 	//push back a pointer to the edge if it has a yflag and is at the boundary of propagation.
		}
	}
	//sort the pointers.
	std::sort(cell_along_boundary.begin(), cell_along_boundary.end());
	//Only keep unique cells
	auto last = std::unique(cell_along_boundary.begin(), cell_along_boundary.end());
	cell_along_boundary.erase(last, cell_along_boundary.end());

	//Convert vector<Cell*> to vector<PCPCell*>
	vector<PCPCell*> pcp_cell_along_boundary;
	pcp_cell_along_boundary.resize(cell_along_boundary.size());
	for (int i = 0; i != cell_along_boundary.size(); i++){
		//Cast the pointer to PCPCell*
		PCPCell* pcpcellptr = dynamic_cast<PCPCell*>(cell_along_boundary.at(i));
		//Sets cell type equal to one ->differentiated and initializes polarity proteins.
		pcp_cell_along_boundary.at(i) = pcpcellptr;
	}
	//Order these cells by angle about the central megacell.
	vector<PCPCell*> ordered_pcp_cell_along_boundary = (*this).order_cells_angle(pcp_cell_along_boundary);
	//Store ordered angles about the central megacell.
	vector<double> ordered_cell_angles = (*this).generate_cell_angle_vector(pcp_cell_along_boundary);
	for (int i = 0; i != cell_along_boundary.size(); i++){
		//Sets cell type equal to one->differentiated and initializes polarity proteins.
		(*ordered_pcp_cell_along_boundary.at(i)).HelperEvolve_InitializePolarityProteins();
	}
}

//WARNING!!! This function assumes that cell.SavedInfoForEvolve_CWVerts has been updated correctly
//Calculates the force on every vertex and stores the result in vertex.Tissue subclasses have their own implementation
//This version of the function does not evolve the boxsize.
void PCPMegaTissue::ForceFunction(ofstream* errorfile){
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double area;
	vector<Vertex*> cwverts;
	Vect FT, FP;

	//CURRENTLY THIS ASSUMES THAT NO EDGE CROSSES THE OUTER DOMAIN BOUNDARY.
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){			//loop through all of the cells in the tissue
		if ((*(*cit)).type == center_cell_type)
			origin_mega_cell = (*(*cit)).GetCenter();
	}

	//First loop through all the cells and find the pressure WARNING!! cell.SavedInfoForEvolve_CWVerts must already be correctly assigned
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
		cwverts = (*(*cit)).SavedInfoForEvolve_CWVerts;							//Get the cw vertices in the cell from the preassigned list 
		area = (*(*cit)).GetArea(cwverts);										//Get the area quickly using the pre assigned cw vertices
		(*(*cit)).SavedInfoForEvolve_Pressure = (*(*cit)).GetPressure(area);	//Get the pressure
	}

	//Calculate the force on each vertex from the pressure and tension
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
		//Get the force from the tension on the edge
		FT = (*(*vit)).GetNetForce_Tension();
		
		//Get the force from the pressure on the edge
		FP = (*(*vit)).GetNetForce_Pressure();
		
		//Add the forces on this edge to the vertex
		(*(*vit)).SavedInfoForEvolve_MechForce = FT + FP;
	}

	//For MegaCell simulations, there is no need to scale the box size: the concentric packing is surrounded by a bounding mega cell.
	double GammaX, GammaY;
	GammaX = 0.0; GammaY = 0.0;

	//Add the force from the tissue stress to every vertex
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
		(*(*vit)).SavedInfoForEvolve_TotalForce.x = (*(*vit)).SavedInfoForEvolve_MechForce.x + (alphaM / TISSUE_ALPHA) * (*(*vit)).r.x * GammaX;
		(*(*vit)).SavedInfoForEvolve_TotalForce.y = (*(*vit)).SavedInfoForEvolve_MechForce.y + (alphaM / TISSUE_ALPHA) * (*(*vit)).r.y * GammaY;
	}

	//Update LY and LX
	ForceLX = GammaX * LX;
	ForceLY = GammaY * LY;
}

//Evolve that includes special functions for mega cell simulation (for y junction formation)
void PCPMegaTissue::Evolve(double interval, double* prot_stepSize, double* mech_stepSize, ofstream* errorfile){
	//find the end point of this evolve->Closest f
	double endtime = floor((CURRENTTIME + interval) / interval) * interval;
	bool flag0, flag1, flag2, flag3;
	bool reinitialize_mechanical = true;
	bool reinitialize_protein = true;
	//Changes by Hayden for PCP: protein and vertex positions are evolved separately (keep track of two times). After evolving protein concentrations, vertex positions are evolved to the corresponding time point.
	int dim_protein;												//The dimension of the mechanical system and the protein system.
	int status;
	vector<double> y_protein;										//Current log variable values.
	//Store parameters for evolution for mechanical and protein in params.
	combinedparameters params;										//parameters for protein evolution.
	//Variable used to calculate protein free energy
	double protein_free_energy;
	//Keep track of total steps taken.
	int total_steps_taken = 0;
	//Specify which order you would like to use: 1, 2, or 3.
	int method_order = 2;


	while (CURRENTTIME < endtime){
		if (time_for_propagation < CURRENTTIME){
			(*this).HelperWriteOutputData_Write2File(errorfile, true);
			//Increment to set next time for propagation of differentiation.
			time_for_propagation += increment_for_propagation;
			//flag to specify whether one should delay to break defects.
			bool delay = false;
			//flag to specify whether one will insert a partial column.
			bool partial_column_insertion = false;
			(*this).HelperEvolve_GeneratePCPProtein(delay, partial_column_insertion);
			//Reinitialize IMEX solver
			reinitialize_IMEX_solver = true;
			reinitialize_mechanical = true;
			reinitialize_protein = true;
		}

		if (time_for_defect_fix < CURRENTTIME){
			(*this).HelperWriteOutputData_Write2File(errorfile, true);
			//Increment to set next time for propagation of differentiation.
			time_for_defect_fix += increment_for_propagation;
			//differentiate cells at the boundary
			(*this).HelperEvolve_GeneratePCPProtein_Defect();
			//Reinitialize IMEX solver
			reinitialize_IMEX_solver = true;
			reinitialize_mechanical = true;
			reinitialize_protein = true;
		}

		//If topology has changed (T1, apoptosis, mitosis, or anything that changed the dimension of the system N and/or changed the ordering of vertices in list), then we need to reinitialize ODE objects
		if (reinitialize_mechanical){
			for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
				(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
			}
			params.vert_id2index = CreateVertexId2IndexMap();		//Get vertex & cell id to index map
			params.cell_id2index = CreateCellId2IndexMap();
			int dim_mechanical = vertices.size() + 1;
			params.dim_mechanical = dim_mechanical;
			//higher order methods:
			double init_val = 0.0;
			vector<vector<Vect>> vec_vect(3, vector<Vect>(params.dim_mechanical, Vect(init_val, init_val)));
			if (reinitialize_IMEX_solver){
				//Initialize the vectors for storing information from past steps.
				f_mech_step_before = vec_vect;
				past_mech_state = vec_vect;
				past_step_sizes.resize(3);
				*mech_stepSize = INITIALSTEPSIZE;
			}
			reinitialize_mechanical = false;
		} //End of initialization

		//If topology has changed (T1, apoptosis, mitosis, or anything that changed the dimension of the system N and/or changed the ordering of vertices in list), then we need to reinitialize ODE objects
		if (reinitialize_protein){
			//Changes by Hayden for PCP: Define a variable to store the number of edges in differentiated cells in the tissue.
			int M = (*this).Find_Number_Differentiated_Edges();
			params.dimension_protein = 2 * M; //Number of dynamical variables associated with protein concentrations.
			dim_protein = params.dimension_protein;
			y_protein.resize(dim_protein);

			//Initialize protein vector with current log variable values.
			(*this).HelperEvolve_Initialize_Protein_Vector(y_protein);

			//higher order methods:
			double init_val = 0.0;
			//The convention for params.y_step_before, params.f_step_before, params.g_step_before, params.past_step_sizes
			//before a step is taken, [0] store info from two steps before, [1] stores info from one step before.
			//after a step is taken, [0] store info from two steps before, [1] stores info from one step before, [2] stores information from the step that just happened.
			//then, within the imex_solver_adaptive, [1]->[0], [2]->[1] such that the information is updated for the next step.
			vector<vector<double>> vec(3, vector<double>(dim_protein, init_val));
			vector<vector<Vect>> vec_vect(3, vector<Vect>(params.dim_mechanical, Vect(init_val, init_val)));

			//TEMP FOR DEBUGGING: only do first order steps for the first two steps. (after that just use information from past steps even if those steps occurred in the last Evolve.
			if (reinitialize_IMEX_solver){
				//Initialize the vectors for storing information from past steps.
				y_step_before = vec;
				f_step_before = vec;
				g_step_before = vec;
				past_step_sizes.resize(3);
				*prot_stepSize = INITIALSTEPSIZE;
			}
			reinitialize_protein = false;		//Finally set reinitialize protein to false
			reinitialize_IMEX_solver = false;		//Set flag for topological changes between instances of evolve to false.
		}

		//STEP 2: Integrate 
		//Loop through the vertices and record the current position in rold 
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++) { (*(*vit)).rold = (*(*vit)).r; }
		LXold = LX;				//Record the old box size
		LYold = LY;

		//Specify whether vertex positions will also evolve.
		bool evolve_mech_positions = true;
		//Evolve protein and vertex positions with IMEX Solver.
		(*this).HelperEvolve_IMEX(y_protein, prot_stepSize, mech_stepSize, error_tolerance_absolute_protein, error_tolerance_absolute_mechanical, &CURRENTTIME, (void*)&params, evolve_mech_positions, method_order, errorfile);

		//Use this for only evolving mechanical positions.
		//(*this).Integrate_IMEX_Without_Protein(mech_stepSize, error_tolerance_absolute_mechanical, &CURRENTTIME, (void*)&params, errorfile);

		//Update protein variables with new protein variables.
		(*this).HelperEvolve_Update_Protein_Variables(y_protein);

		if (evolve_mech_positions){
			//Step 3: Update
			LX = LXnew;					//Update the current box size
			LY = LYnew;
			(*this).HelperEvolve_UpdatePosition((void*)&params);  //Update the vertex positions (this automaticaly updates edge flags if a vertex has crossed the boundary)
		}

		total_steps_taken++;					//Keep track of total steps taken.

		//STEP 4: look for topological changes and set flags if any have occured 
		//No apoptosis in my simulations
		flag1 = 0;
		//Somehow let the code know if either of these occurred on differentiated edges.
		if (evolve_mech_positions){
			flag2 = (*this).T1(errorfile);
			flag3 = (*this).ResolveVertices(errorfile);
		}
		else{
			flag2 = 0;
			flag3 = 0;
		}
		if (flag1 || flag2 || flag3){
			cout << "reinitialize protein" << endl;
			reinitialize_protein = true;
			reinitialize_mechanical = true;		//initialize mechanical and protein systems separately.
			reinitialize_IMEX_solver = true;	//use external variable to check whether we need to reinitialize_IMEX_solver (because may need to reinitialize upon reentering evolve).
		}
	}

	cout << total_steps_taken << endl;

	cout << "TISSUE EVOLVED TO TIME " << CURRENTTIME << "  LX = " << LX << "  LY = " << LY << endl;

}

//HELPER FUNCTIONS

//! Helper function for checking the cell input conforms to expectations.
void PCPMegaTissue::HelperReadInputData_CheckTypes(){
	list<Cell*>::iterator cit;
	int number_of_central_cells = 0;
	int number_of_outer_cells = 0;
	int expected_number_of_central_cells = 1;
	int expected_number_of_outer_cells = 1;
	for (cit = cells.begin(); cit != cells.end(); cit++){	//look through all of the cells
		if ((*(*cit)).type == center_cell_type){
			number_of_central_cells++;
		}
		if ((*(*cit)).type == outer_cell_type){
			number_of_outer_cells++;
		}
		//Unidentified cell type (exit code).
		if ((*(*cit)).type != standard_cell_type && (*(*cit)).type != center_cell_type && (*(*cit)).type != outer_cell_type){
			cout << "Unexpected cell type. Code will exit." << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (number_of_central_cells != expected_number_of_central_cells){
		cout << "There is more than one central cell. Code will exit." << endl;
		exit(EXIT_FAILURE);
	}
	if (number_of_outer_cells != expected_number_of_outer_cells){
		cout << "There is more than one outer cell. Code will exit." << endl;
		exit(EXIT_FAILURE);
	}
}

// Reads the data and sets all of the global variables
// Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
void PCPMegaTissue::HelperReadInputData_SetGlobalValues(ifstream* datafile){

	string line, word, key, value;		//temporary variables to store lines and words 
	string type;						//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline;				//ssline is a string stream so we can comvert strings to stringstreams
	map <string, string> globals;



	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//

	//Read until you get to the word Global signifying the begining of the global variables
	while (line.compare("Global:") != 0){
		getline(*datafile, line);
	}

	//Now read the global variables and make a map with the keywords
	while (line.compare("endGlobal") != 0){
		getline(*datafile, line);							//get the next line
		ssline.clear();
		ssline << line;										//put it into a stream
		ssline >> key;										//the first word is the key
		ssline >> value;									//the second is the value of the key
		if (value != ""){									//if the item actually has a value
			globals[key] = value;								//stick the item and it's value in the globals map
		}
	}

	//Now assign the constants their values
	//makes sure the key is actually part of the map before giving the constant a value if the key is not defined it will have the value given at the begining of main
	// the atof is there to convert the string to a double
	if (globals.count("LX") == 1){ LX = atof(globals["LX"].c_str()); }
	else if (globals.count("LX") == 0){ cout << "WARNING: LX is not given in input file, used default value " << LX << endl; }
	else { cout << "ERROR: Global variable LX has multiple definition in input file, used default value LX = " << LX << endl; }

	if (globals.count("LY") == 1){ LY = atof(globals["LY"].c_str()); }
	else if (globals.count("LY") == 0){ cout << "WARNING: LY is not given in input file, used default value " << LY << endl; }
	else { cout << "ERROR: Global variable LX has multiple definition in input file, used default value LY = " << LY << endl; }

	if (globals.count("MINAREA") == 1){ MINAREA = atof(globals["MINAREA"].c_str()); }
	else if (globals.count("MINAREA") == 0){ cout << "WARNING: MINAREA is not given in input file, used default value " << MINAREA << endl; }
	else { cout << "ERROR: Global variable MINAREA has multiple definition in input file, used default value MINAREA = " << MINAREA << endl; }

	if (globals.count("MINLENGTH") == 1){ MINLENGTH = atof(globals["MINLENGTH"].c_str()); }
	else if (globals.count("MINLENGTH") == 0){ cout << "WARNING: LX is not given in input file, used default value " << MINLENGTH << endl; }
	else { cout << "ERROR: Global variable MINLENGTH has multiple definition in input file, used default value MINLENGTH = " << MINLENGTH << endl; }

	if (globals.count("NEWEDGELENGTH") == 1){ NEWEDGELENGTH = atof(globals["NEWEDGELENGTH"].c_str()); }
	else if (globals.count("NEWEDGELENGTH") == 0){ cout << "WARNING: NEWEDGELENGTH is not given in input file, used default value " << NEWEDGELENGTH << endl; }
	else { cout << "ERROR: Global variable NEWEDGELENGTH has multiple definition in input file, used default value NEWEDGELENGTH = " << NEWEDGELENGTH << endl; }

	if (globals.count("GAMMA") == 1){ GAMMA = atof(globals["GAMMA"].c_str()); }
	else if (globals.count("GAMMA") == 0){ cout << "WARNING: GAMMA is not given in input file, used default value " << GAMMA << endl; }
	else { cout << "ERROR: Global variable GAMMA has multiple definition in input file, used default value GAMMA = " << GAMMA << endl; }

	if (globals.count("K") == 1){ K = atof(globals["K"].c_str()); }
	else if (globals.count("K") == 0){ cout << "WARNING: K is not given in input file, used default value " << K << endl; }
	else { cout << "ERROR: Global variable K has multiple definition in input file, used default value K = " << K << endl; }

	if (globals.count("SEED") == 1){ SEED = atof(globals["SEED"].c_str()); }
	else if (globals.count("SEED") == 0){ cout << "WARNING: SEED is not given in input file, used default value " << SEED << endl; }
	else { cout << "ERROR: Global variable SEED has multiple definition in input file, used default value SEED = " << SEED << endl; }

	if (globals.count("INITIALSTEPSIZE") == 1){ INITIALSTEPSIZE = atof(globals["INITIALSTEPSIZE"].c_str()); }
	else if (globals.count("INITIALSTEPSIZE") == 0){ cout << "WARNING: INITIALSTEPSIZE is not given in input file, used default value " << INITIALSTEPSIZE << endl; }
	else { cout << "ERROR: Global variable INITIALSTEPSIZE has multiple definition in input file, used default value INITIALSTEPSIZE = " << INITIALSTEPSIZE << endl; }

	if (globals.count("TIMEINTERVAL") == 1){ TIMEINTERVAL = atof(globals["TIMEINTERVAL"].c_str()); }
	else if (globals.count("TIMEINTERVAL") == 0){ cout << "WARNING: TIMEINTERVAL is not given in input file, used default value " << TIMEINTERVAL << endl; }
	else { cout << "ERROR: Global variable TIMEINTERVAL has multiple definition in input file, used default value TIMEINTERVAL = " << TIMEINTERVAL << endl; }

	if (globals.count("TOTALTIME") == 1){ TOTALTIME = atof(globals["TOTALTIME"].c_str()); }
	else if (globals.count("TOTALTIME") == 0){ cout << "WARNING: TOTALTIME is not given in input file, used default value " << TOTALTIME << endl; }
	else { cout << "ERROR: Global variable TOTALTIME has multiple definition in input file, used default value TOTALTIME = " << TOTALTIME << endl; }

	if (globals.count("CURRENTTIME") == 1){ CURRENTTIME = atof(globals["CURRENTTIME"].c_str()); }
	else if (globals.count("CURRENTTIME") == 0){ cout << "WARNING: CURRENTTIME is not given in input file, used default value " << CURRENTTIME << endl; }
	else { cout << "ERROR: Global variable CURRENTTIME has multiple definition in input file, used default value CURRENTTIME = " << CURRENTTIME << endl; }

	if (globals.count("error_tolerance_absolute_mechanical") == 1){ error_tolerance_absolute_mechanical = atof(globals["error_tolerance_absolute_mechanical"].c_str()); }
	else if (globals.count("error_tolerance_absolute_mechanical") == 0){ cout << "WARNING: error_tolerance_absolute_mechanical is not given in input file, used default value " << error_tolerance_absolute_mechanical << endl; }
	else { cout << "ERROR: Global variable error_tolerance_absolute_mechanical has multiple definition in input file, used default value error_tolerance_absolute_mechanical = " << error_tolerance_absolute_mechanical << endl; }

	if (globals.count("error_tolerance_absolute_protein") == 1){ error_tolerance_absolute_protein = atof(globals["error_tolerance_absolute_protein"].c_str()); }
	else if (globals.count("error_tolerance_absolute_protein") == 0){ cout << "WARNING: error_tolerance_absolute_protein is not given in input file, used default value " << error_tolerance_absolute_protein << endl; }
	else { cout << "ERROR: Global variable error_tolerance_absolute_protein has multiple definition in input file, used default value error_tolerance_absolute_protein = " << error_tolerance_absolute_protein << endl; }

	if (globals.count("ALPHAM") == 1){ alphaM = atof(globals["ALPHAM"].c_str()); }
	else if (globals.count("ALPHAM") == 0){ cout << "WARNING: ALPHAM is not given in input file, used default value " << alphaM << endl; }
	else { cout << "ERROR: Global variable ALPHAM has multiple definition in input file, used default value ALPHAM = " << alphaM << endl; }

	if (globals.count("INNERPRESSURE") == 1){ inner_pressure = atof(globals["INNERPRESSURE"].c_str()); }
	else if (globals.count("INNERPRESSURE") == 0){ cout << "WARNING: INNERPRESSURE is not given in input file, used default value " << inner_pressure << endl; }
	else { cout << "ERROR: Global variable INNERPRESSURE has multiple definition in input file, used default value INNERPRESSURE = " << inner_pressure << endl; }

	if (globals.count("OUTERPRESSURE") == 1){ outer_pressure = atof(globals["OUTERPRESSURE"].c_str()); }
	else if (globals.count("OUTERPRESSURE") == 0){ cout << "WARNING: OUTERPRESSURE is not given in input file, used default value " << outer_pressure << endl; }
	else { cout << "ERROR: Global variable OUTERPRESSURE has multiple definition in input file, used default value OUTERPRESSURE = " << outer_pressure << endl; }

	if (globals.count("EPSAB") == 1){ epsab = atof(globals["EPSAB"].c_str()); }
	else if (globals.count("EPSAB") == 0){ cout << "WARNING: EPSAB is not given in input file, used default value " << epsab << endl; }
	else { cout << "ERROR: Global variable EPSAB has multiple definition in input file, used default value EPSAB = " << epsab << endl; }

	if (globals.count("EPSNAUGHT") == 1){ epsnaught = atof(globals["EPSNAUGHT"].c_str()); }
	else if (globals.count("EPSNAUGHT") == 0){ cout << "WARNING: EPSNAUGHT is not given in input file, used default value " << epsnaught << endl; }
	else { cout << "ERROR: Global variable EPSNAUGHT has multiple definition in input file, used default value EPSNAUGHT = " << epsnaught << endl; }

	if (globals.count("EPSONE") == 1){ epsone = atof(globals["EPSONE"].c_str()); }
	else if (globals.count("EPSONE") == 0){ cout << "WARNING: EPSONE is not given in input file, used default value " << epsone << endl; }
	else { cout << "ERROR: Global variable EPSONE has multiple definition in input file, used default value EPSONE = " << epsone << endl; }

	if (globals.count("EPSTWO") == 1){ epstwo = atof(globals["EPSTWO"].c_str()); }
	else if (globals.count("EPSTWO") == 0){ cout << "WARNING: EPSTWO is not given in input file, used default value " << epstwo << endl; }
	else { cout << "ERROR: Global variable EPSTWO has multiple definition in input file, used default value EPSTWO = " << epstwo << endl; }

	if (globals.count("JONE") == 1){ jone = atof(globals["JONE"].c_str()); }
	else if (globals.count("JONE") == 0){ cout << "WARNING: JONE is not given in input file, used default value " << jone << endl; }
	else { cout << "ERROR: Global variable JONE has multiple definition in input file, used default value JONE = " << jone << endl; }

	if (globals.count("JTWO") == 1){ jtwo = atof(globals["JTWO"].c_str()); }
	else if (globals.count("JTWO") == 0){ cout << "WARNING: JTWO is not given in input file, used default value " << jtwo << endl; }
	else { cout << "ERROR: Global variable JTWO has multiple definition in input file, used default value JTWO = " << jtwo << endl; }

	if (globals.count("TE") == 1){ te = atof(globals["TE"].c_str()); }
	else if (globals.count("TE") == 0){ cout << "WARNING: TE is not given in input file, used default value " << te << endl; }
	else { cout << "ERROR: Global variable TE has multiple definition in input file, used default value TE = " << te << endl; }

	if (globals.count("KT1") == 1){ kT1 = atof(globals["KT1"].c_str()); }
	else if (globals.count("KT1") == 0){ cout << "WARNING: KT1 is not given in input file, used default value " << kT1 << endl; }
	else { cout << "ERROR: Global variable KT1 has multiple definition in input file, used default value KT1 = " << kT1 << endl; }

	if (globals.count("GAMMA0") == 1){ GAMMA0 = atof(globals["GAMMA0"].c_str()); }
	else if (globals.count("GAMMA0") == 0){ cout << "WARNING: GAMMA0 is not given in input file, used default value " << GAMMA0 << endl; }
	else { cout << "ERROR: Global variable GAMMA0 has multiple definition in input file, used default value GAMMA0 = " << GAMMA0 << endl; }

	if (globals.count("GAMMAC") == 1){ GAMMAC = atof(globals["GAMMAC"].c_str()); }
	else if (globals.count("GAMMAC") == 0){ cout << "WARNING: GAMMAC is not given in input file, used default value " << GAMMAC << endl; }
	else { cout << "ERROR: Global variable GAMMAC has multiple definition in input file, used default value GAMMAC = " << GAMMAC << endl; }

	if (globals.count("ALPHAP") == 1){ alphaP = atof(globals["ALPHAP"].c_str()); }
	else if (globals.count("ALPHAP") == 0){ cout << "WARNING: ALPHAP is not given in input file, used default value " << alphaP << endl; }
	else { cout << "ERROR: Global variable ALPHAP has multiple definition in input file, used default value ALPHAP = " << alphaP << endl; }

	if (globals.count("INCREMENT_FOR_PROPAGATION") == 1){ increment_for_propagation = atof(globals["INCREMENT_FOR_PROPAGATION"].c_str()); }
	else if (globals.count("INCREMENT_FOR_PROPAGATION") == 0){ cout << "WARNING: INCREMENT_FOR_PROPAGATION is not given in input file, used default value " << increment_for_propagation << endl; }
	else { cout << "ERROR: Global variable INCREMENT_FOR_PROPAGATION has multiple definition in input file, used default value INCREMENT_FOR_PROPAGATION = " << increment_for_propagation << endl; }

	if (globals.count("INCREMENT_FOR_RENORMALIZATION") == 1){ increment_for_renormalization = atof(globals["INCREMENT_FOR_RENORMALIZATION"].c_str()); }
	else if (globals.count("INCREMENT_FOR_RENORMALIZATION") == 0){ cout << "WARNING: INCREMENT_FOR_RENORMALIZATION is not given in input file, used default value " << increment_for_renormalization << endl; }
	else { cout << "ERROR: Global variable INCREMENT_FOR_RENORMALIZATION has multiple definition in input file, used default value INCREMENT_FOR_RENORMALIZATION = " << increment_for_renormalization << endl; }

	if (globals.count("MAX_STEP_SIZE") == 1){ max_step_size = atof(globals["MAX_STEP_SIZE"].c_str()); }
	else if (globals.count("MAX_STEP_SIZE") == 0){ cout << "WARNING: MAX_STEP_SIZE is not given in input file, used default value " << max_step_size << endl; }
	else { cout << "ERROR: Global variable MAX_STEP_SIZE has multiple definition in input file, used default value MAX_STEP_SIZE = " << max_step_size << endl; }


	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//
}

// Writes data into file specified by 'datafile'. If flag 'includeExtraInfo' is true, then additional information (such as tension, cell area, etc) is included.
// The extra information is non-essential in reconstructing the simulation state, and is just provided for convenience (e.g. for display in MATLAB).
void PCPMegaTissue::HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo){
	time_t rawtime;					//to deal with writing the date out
	struct tm * timeinfo;
	(*datafile).precision(10);
	(*datafile).setf(ios::scientific);
	(*datafile).setf(ios::showpoint);

	//First write the info about what you just ran and the time you ran it
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;
	*datafile << "Date:\t" << asctime(timeinfo) << endl;

	//Now write the global variables you used
	*datafile << "Global:" << endl;
	*datafile << "SEED\t" << SEED << endl;
	*datafile << "LX\t" << LX << endl;
	*datafile << "LY\t" << LY << endl;
	*datafile << "MINAREA\t" << MINAREA << endl;
	*datafile << "MINLENGTH\t" << MINLENGTH << endl;
	*datafile << "NEWEDGELENGTH\t" << NEWEDGELENGTH << endl;
	*datafile << "K\t" << K << endl;
	*datafile << "GAMMA\t" << GAMMA << endl;
	*datafile << "TIMEINTERVAL\t" << TIMEINTERVAL << endl;
	*datafile << "TOTALTIME\t" << TOTALTIME << endl;
	*datafile << "CURRENTTIME\t" << CURRENTTIME << endl;
	*datafile << "MAX_STEP_SIZE\t" << max_step_size << endl;
	*datafile << "INITIALSTEPSIZE\t" << INITIALSTEPSIZE << endl;
	*datafile << "error_tolerance_absolute_mechanical\t" << error_tolerance_absolute_mechanical << endl;
	*datafile << "error_tolerance_absolute_protein\t" << error_tolerance_absolute_protein << endl;
	*datafile << "INCREMENT_FOR_PROPAGATION\t" << increment_for_propagation << endl;
	*datafile << "INCREMENT_FOR_RENORMALIZATION\t" << increment_for_renormalization << endl;
	*datafile << "ALPHAP\t" << alphaP << endl;
	*datafile << "ALPHAM\t" << alphaM << endl;
	*datafile << "GAMMA0\t" << GAMMA0 << endl;
	*datafile << "GAMMAC\t" << GAMMAC << endl;
	*datafile << "EPSAB\t" << epsab << endl;
	*datafile << "EPSNAUGHT\t" << epsnaught << endl;
	*datafile << "EPSONE\t" << epsone << endl;
	*datafile << "EPSTWO\t" << epstwo << endl;
	*datafile << "TE\t" << te << endl;
	*datafile << "JONE\t" << jone << endl;
	*datafile << "JTWO\t" << jtwo << endl;
	*datafile << "KT1\t" << kT1 << endl;
	*datafile << "INNERPRESSURE\t" << inner_pressure << endl;
	*datafile << "OUTERPRESSURE\t" << outer_pressure << endl;
	//******ADD NEW OUTPUT HERE******//
	*datafile << "endGlobal" << endl;

	if (!includeExtraInfo)											// if flag is false, include only essential information to capture the state of the simulation
	{
		//print the vertex info
		*datafile << "Vertices" << endl;
		*datafile << "id\tx\ty" << endl;
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){			//loop through all of the vertices in the tissue
			*datafile << (*(*vit)).id << "\t" << (*(*vit)).r.x << "\t" << (*(*vit)).r.y << endl;						//write it's info to the file
		}

		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			*datafile << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << endl;			//write its info to the file
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;
		*datafile << "id\tA0\tEdgeNumb\tEdgeList\tType\tPartnerCell" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++){				//make comma-separated list of edges
				ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
			}

			PCPCell* pcp_cit = dynamic_cast<PCPCell*>(*cit);
			*datafile << (*pcp_cit).id << "\t" << (*pcp_cit).A0 << "\t" << (*pcp_cit).edges.size() << "\t" << ssEdgeList.str() << "\t" << (*pcp_cit).type << endl;
		}
	}

	else															// if flag is true, include only essential information to capture the state of the simulation
	{
		//print the vertex info with additional information
		*datafile << "Vertices" << endl;
		*datafile << "id\tx\ty\t_EdgesNumber" << endl;
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){			//loop through all of the vertices in the tissue
			*datafile << (*(*vit)).id << "\t" << (*(*vit)).r.x << "\t" << (*(*vit)).r.y << "\t" << (*(*vit)).edges.size() << endl;						//write it's info to the file
		}

		//Changes by Hayden for PCP: output PCP concentrations.
		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag\t_length\t_cell1\t_cell2\t_cell1AProtein\t_cell1BProtein\t_cell2AProtein\t_cell2BProtein" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			double L = (*(*eit)).GetLength();
			PCPEdge* pcpeit = dynamic_cast<PCPEdge*>(*eit);
			vector<double> planar_polarity_protein = (*pcpeit).planarpolarity;
			*datafile << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*pcpeit).yflag << "\t" << L
				<< "\t" << (*(*pcpeit).cell1).id << "\t" << (*(*pcpeit).cell2).id << "\t" << planar_polarity_protein[0] << "\t" << planar_polarity_protein[1] << "\t" << planar_polarity_protein[2] << "\t" << planar_polarity_protein[3] << endl;										//write its info to the file
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;
		*datafile << "id\tA0\tEdgeNumb\tType\t_xflag\t_yflag\t_Area\t_Pressure\t_Center_x\t_Center_y\tEdgeList\t_CWVertices\t_OrderedFlagsX\t_OrderedFlagsY\tPartnerCell" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++){				//make comma-separated list of edges
				ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
			}

			vector<Vertex*> CWVertices = (*(*cit)).GetClockWiseVertices();
			stringstream ssCWVertices;
			ssCWVertices.clear();
			for (i = 0; i < CWVertices.size(); i++){			//make comma-separated list of vertices ordered clockwise
				ssCWVertices << (*CWVertices[i]).id << ",";
			}

			//Changes by Hayden for PCP: Outputting this is useful for calculating shape polarization vectors in MATLAB.
			vector<vector<int>> CWVertexFlags = (*(*cit)).HelperOrderFlags(CWVertices);
			stringstream ssCWVerticesFlagsX;
			ssCWVerticesFlagsX.clear();
			stringstream ssCWVerticesFlagsY;
			ssCWVerticesFlagsY.clear();
			for (i = 0; i < CWVertices.size(); i++){			//make comma-separated list of vertices ordered clockwise
				ssCWVerticesFlagsX << CWVertexFlags[0][i] << ",";
			}
			for (i = 0; i < CWVertices.size(); i++){			//make comma-separated list of vertices ordered clockwise
				ssCWVerticesFlagsY << CWVertexFlags[1][i] << ",";
			}

			Vect centroid = (*(*cit)).GetCenter();
			PCPCell* pcp_cit = dynamic_cast<PCPCell*>(*cit);
			*datafile << (*pcp_cit).id << "\t" << (*pcp_cit).A0 << "\t" << (*pcp_cit).edges.size() << "\t" << (*pcp_cit).type << "\t" << (*pcp_cit).GetXFlag() << "\t" << (*pcp_cit).GetYFlag() << "\t"
				<< (*pcp_cit).GetArea() << "\t" << (*pcp_cit).GetPressure((*pcp_cit).GetArea()) << "\t" << centroid.x << "\t" << centroid.y << "\t" << ssEdgeList.str() << "\t" << ssCWVertices.str() << "\t" << ssCWVerticesFlagsX.str() << "\t" << ssCWVerticesFlagsY.str() << endl;
		}
	}


	*datafile << "End" << endl << endl;
}

//Locate trapped defects (often the result of short edge that has not undergone a T1).
vector<PCPCell*> PCPMegaTissue::locate_trapped_defect(vector<PCPCell*> next_differentiated_cells){
	//number of cells for differentiation
	int number_of_next_diff_cells = next_differentiated_cells.size();
	//Specify type 3 for cell that are next for PCP.
	for (int cell_index = 0; cell_index != number_of_next_diff_cells; cell_index++){
		(*next_differentiated_cells.at(cell_index)).type = pcp_cell_type;
	}

	//order cells by angle about the central megacell
	vector<PCPCell*> ordered_cells = (*this).order_cells_angle(next_differentiated_cells);
	//Vector for storing defects trapped by short edge.
	vector<int> located_defects;
	//Vector for storing cells that neighbor all cells that are neighbors to a given cell.
	vector<Cell*> neighboring_cells;
	vector<PCPCell*> defect_cells;
	//Count neighbors that are of type 0.
	int count_zero_type = 0;

	//Search for trapped defects by counting neighboring cells of type 0.
	for (int cell_index = 0; cell_index != number_of_next_diff_cells; cell_index++){
		count_zero_type = 0;
		//Find neighboring cells.
		neighboring_cells = (*ordered_cells.at(cell_index)).GetNeighboringCells();
		//Count number of neighboring cells.
		int number_of_neighboring_cells = neighboring_cells.size();

		//Loop through neighboring cells.
		for (int search_cell_index = 0; search_cell_index != number_of_neighboring_cells; search_cell_index++){
			//Increment counter if neighboring cell of type 0.
			if ((*neighboring_cells.at(search_cell_index)).type == 0)
				count_zero_type++;
		}
		//If all cells surrounding a given cell are of types other than 0, store the defect.
		if (count_zero_type == 0)
			located_defects.push_back(cell_index);
	}

	//Loop through located defects: push deactivated cells (one of the cells neighboring the defect)
	for (int cell_index = 0; cell_index != located_defects.size(); cell_index++){
		defect_cells.push_back(ordered_cells.at(located_defects.at(cell_index)));
	}

	//Reset next differentiated cell type to 0.
	for (int cell_index = 0; cell_index != number_of_next_diff_cells; cell_index++){
		(*next_differentiated_cells.at(cell_index)).type = 0;
	}

	return defect_cells;
}

//Locate trapped defects (often the result of short edge that has not undergone a T1).
vector<PCPCell*> PCPMegaTissue::locate_cells_for_deactivation(vector<PCPCell*> next_differentiated_cells){
	//number of cells for differentiation
	int number_of_next_diff_cells = next_differentiated_cells.size();
	//Specify type 3 for cell that are next for PCP.
	for (int cell_index = 0; cell_index != number_of_next_diff_cells; cell_index++){
		(*next_differentiated_cells.at(cell_index)).type = pcp_cell_type;
	}

	//order cells by angle about the central megacell
	vector<PCPCell*> ordered_cells = (*this).order_cells_angle(next_differentiated_cells);
	//Vector for storing defects trapped by short edge.
	vector<int> located_defects;
	//Vector for storing cells that neighbor all cells that are neighbors to a given cell.
	vector<Cell*> neighboring_cells;
	//Store cells that will be deactivated to avoid a trapped defect.
	vector<PCPCell*> deactivated_cells;
	//Count neighbors that are of type 0.
	int count_zero_type = 0;

	//Search for trapped defects by counting neighboring cells of type 0.
	for (int cell_index = 0; cell_index != number_of_next_diff_cells; cell_index++){
		count_zero_type = 0;
		//Find neighboring cells.
		neighboring_cells = (*ordered_cells.at(cell_index)).GetNeighboringCells();
		//Count number of neighboring cells.
		int number_of_neighboring_cells = neighboring_cells.size();

		//Loop through neighboring cells.
		for (int search_cell_index = 0; search_cell_index != number_of_neighboring_cells; search_cell_index++){
			//Increment counter if neighboring cell of type 0.
			if ((*neighboring_cells.at(search_cell_index)).type == 0)
				count_zero_type++;
		}
		//If all cells surrounding a given cell are of types other than 0, store the defect.
		if (count_zero_type == 0)
			located_defects.push_back(cell_index);
	}

	//Loop through located defects: push deactivated cells (one of the cells neighboring the defect)
	for (int cell_index = 0; cell_index != located_defects.size(); cell_index++){
		int previous_cell_index = (located_defects.at(cell_index) - 1) % number_of_next_diff_cells;
		int next_cell_index = (located_defects.at(cell_index) + 1) % number_of_next_diff_cells;
		//Perhaps, choose next, previous with 50-50 chance.
		deactivated_cells.push_back(ordered_cells.at(next_cell_index));
	}

	//Reset next differentiated cell type to 0.
	for (int cell_index = 0; cell_index != number_of_next_diff_cells; cell_index++){
		(*next_differentiated_cells.at(cell_index)).type = 0;
	}

	return deactivated_cells;
}

// Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids.
// Starts at current position of *datafile and reads until keyword 'End'.
void PCPMegaTissue::HelperReadInputData_CreateObjects(ifstream* datafile){
	string line, key, value, edg;					//temporary variables to store lines and words 
	string type;									//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline, ssedges;					//ssline is a string stream so we can comvert strings to stringstreams
	vector<string> keywords;						//vector of strings to hold the keywords
	vector<int> listedges;
	int i;
	map <string, string> objectInfo;				//map to hold the keywords and their data
	int newid, v1, v2, enumb, typenumb;				//temp varables to store the object info
	double A0;
	double x, y;
	bool xf, yf;


	//VERTICES
	//read the file untill you get to the word vertices
	while (line.compare("Vertices") != 0){
		getline(*datafile, line);
	}
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}
	//now recorde the data and make the vertex object untill you get to the word edges signifying the end of the data
	while (true){
		getline(*datafile, line);					//get the next line
		if ((line.compare("Edges") == 0)){ break; }
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')){				//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//make the object and as it to the linked list of vertices also record it's id in the master id list
		newid = atoi(objectInfo["id"].c_str());				//extract the data from the map
		x = atof(objectInfo["x"].c_str());
		y = atof(objectInfo["y"].c_str());
		Vertex* vptr = new Vertex;					// Create a vertes with an id and x,y pos and put it in the master list of all veritces
		(*vptr).id = newid;
		(*vptr).r.x = x;
		(*vptr).r.y = y;
		(vertices).push_back(vptr);
		vertexIDs.push_back(newid);					//and add its id to the list of ids
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*vertices.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}

	//EDGES
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}
	//now recorde the data and make the edge object untill you get to the word cells signifying the end of the data
	while (true){
		getline(*datafile, line);					//get the next line
		if (line.compare("Cells") == 0){ break; }
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')){				//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//make the object and as it to the linked list of vertices also record it's id in the master id list
		newid = atoi(objectInfo["id"].c_str());				//extract the data from the map
		v1 = atoi(objectInfo["vertex1"].c_str());
		v2 = atoi(objectInfo["vertex2"].c_str());
		xf = atoi(objectInfo["xflag"].c_str());
		yf = atoi(objectInfo["yflag"].c_str());
		PCPEdge* eptr = new PCPEdge;
		(*eptr).id = newid;
		(*eptr).IOVertexID[0] = v1;
		(*eptr).IOVertexID[1] = v2;
		(*eptr).xflag = xf;
		(*eptr).yflag = yf;
		(*eptr).cell1 = NULL;
		(*eptr).cell2 = NULL;
		//Initialize PCP concentrations to be equal to zero.
		(*eptr).planarpolarity.resize(4);
		fill((*eptr).planarpolarity.begin(), (*eptr).planarpolarity.end(), 0.0);
		(edges).push_back(eptr);					// Create a vertes with an id and x,y pos and put it in the master list of all veritces
		edgeIDs.push_back(newid);					//and add its id to the list of ids
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*edges.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}

	//CELLS
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}
	//now recorde the data and make the cell object untill you get to the word end signifying the end of the data
	while (true){
		getline(*datafile, line);					//get the next line
		if (line.compare("End") == 0){ break; }				//break if keyworkd 'End' is found
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')){				//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//extract the data from the map
		newid = atoi(objectInfo["id"].c_str());
		enumb = atoi(objectInfo["EdgeNumb"].c_str());
		if (objectInfo.count("Type")){
			typenumb = atoi(objectInfo["Type"].c_str());
		}
		else{ typenumb = 1; }
		if (objectInfo.count("A0")){
			A0 = atof(objectInfo["A0"].c_str());
		}
		ssedges.clear();
		ssedges << objectInfo["EdgeList"];
		//now worry about the edge list
		listedges.clear();
		while (getline(ssedges, edg, ',')){
			listedges.push_back(atoi(edg.c_str()));
		}
		//make the cell object
		PCPMegaCell* cptr = new PCPMegaCell;
		(*cptr).id = newid;
		(*cptr).type = typenumb;
		(*cptr).A0 = A0;
		(*cptr).IOEdgeID = listedges;
		(*cptr).time_of_pcp = 0.0;
		(*cptr).cell_generation = -1;
		cells.push_back(cptr);
		cellIDs.push_back(newid);
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}
}

//! This function is used to loop through the list of undifferentiated cells that are in contact with differentiated cells at a propagation time.
// \note Add a partial column after a couple of full columns.
void PCPMegaTissue::HelperEvolve_GeneratePCPProtein(bool delay, bool added_partial_column){
	//Find boundary cells (all cells of type 0 in contact with a cell of type 3)
	vector<PCPCell*> cells_for_differentiation = (*this).HelperEvolve_FindBoundaryCells();
	vector<PCPCell*> deactivated_cells = (*this).locate_cells_for_deactivation(cells_for_differentiation);
	vector<PCPCell*> deactivated_cells_with_cutoff;

	//Store number of cells for next round of differentiation.
	int number_of_cells_for_differentiation = cells_for_differentiation.size();

	//Order cells for differentiation by angle about the central megacell.
	vector<PCPCell*> ordered_cells_for_differentiation = (*this).order_cells_angle(cells_for_differentiation);
	//Store corresponding cell angles (this follows the convention of the ordered_cells_for_differentiation)
	vector<double> ordered_cell_angles = (*this).generate_cell_angle_vector(cells_for_differentiation);

	//Output array of contact for row definition.
	vector<PCPCell*> pcp_cell_along_boundary = (*this).HelperEvolve_FindCells_Contact_with_Boundary();
	vector<PCPCell*> ordered_pcp_cell_along_boundary = (*this).order_cells_angle(pcp_cell_along_boundary);

	//Define the sin_cut_off, any cell with sin(angle about central mega cell) should differentiate.
	//Effectively, this adds a partial column.
	double sin_cut_off = -0.5;

	//After you have added two full columns, add a partial column (based on sin(angle about the central megacell)).
	if ((differentiation_count == partial_column_insertion) && added_partial_column){
		if (delay){
			//Loop through all cell_for_differentiation.
			for (int cell_index = 0; cell_index != number_of_cells_for_differentiation; cell_index++){
				//Check sin of angle about central mega cell.
				if ((sin(ordered_cell_angles.at(cell_index)) >= sin_cut_off) && (!check_for_cell_pointer_within_list(ordered_cells_for_differentiation[cell_index], deactivated_cells))){
					//Initialize the PCP protein concentrations in the list of cells for differentiation.
					(*ordered_cells_for_differentiation[cell_index]).HelperEvolve_InitializePolarityProteins();
				}
				if ((check_for_cell_pointer_within_list(ordered_cells_for_differentiation[cell_index], deactivated_cells)) && (sin(ordered_cell_angles.at(cell_index)) >= sin_cut_off)){
					deactivated_cells_with_cutoff.push_back(ordered_cells_for_differentiation[cell_index]);
				}
			}
		}
		else{
			//Loop through all cell_for_differentiation.
			for (int cell_index = 0; cell_index != number_of_cells_for_differentiation; cell_index++){
				//Check sin of angle about central mega cell.
				if (sin(ordered_cell_angles.at(cell_index)) >= sin_cut_off){
					//Initialize the PCP protein concentrations in the list of cells for differentiation.
					(*ordered_cells_for_differentiation[cell_index]).HelperEvolve_InitializePolarityProteins();
				}
			}
		}
	}
	else{
		if (delay){
			//Loop through all cell_for_differentiation (don't do a sin cut off)
			for (int cell_index = 0; cell_index != number_of_cells_for_differentiation; cell_index++){
				if (!check_for_cell_pointer_within_list(ordered_cells_for_differentiation[cell_index], deactivated_cells))
					(*ordered_cells_for_differentiation[cell_index]).HelperEvolve_InitializePolarityProteins();
			}
		}
		else{
			//Loop through all cell_for_differentiation (don't do a sin cut off)
			for (int cell_index = 0; cell_index != number_of_cells_for_differentiation; cell_index++){
				(*ordered_cells_for_differentiation[cell_index]).HelperEvolve_InitializePolarityProteins();
			}
		}
	}
	if (delay){
		if ((differentiation_count == partial_column_insertion) && added_partial_column){
			fix_defects = deactivated_cells_with_cutoff;
		}
		else{
			fix_defects = deactivated_cells;
		}
	}

	//Increment the differentiation count.
	differentiation_count++;
}

//! This is a helper function for generating a checkerboard pattern in the next row: for a circular ring with a vertical stripe rather than for periodic boundary conditions.
void PCPMegaTissue::HelperCheckerboard_GeneratePattern(vector<PCPCell*> next_differentiated_cells, ofstream* datafile){
	//Order the differentiated cells based on angle: this helps later for flipping domains of spin.
	vector<PCPCell*> ordered_next_differentiated_cells = order_cells_angle(next_differentiated_cells);
	int number_of_next_differentiated_cells = ordered_next_differentiated_cells.size();
	//Generate random distribution of spins along the ring.
	vector<double> next_diff_cells_types = HelperCheckerboard_GenerateRandomSpins(number_of_next_differentiated_cells);

	//Build list of all differentiated cells connected to the next_differentiated_cells (This list includes next_differentiated_cells).
	vector<Cell*> connected_cells = HelperCheckerboard_BuildList_ConnectedCells(ordered_next_differentiated_cells);
	int number_of_connected_cells = connected_cells.size();

	//Build adjacency matrix that includes contact lengths between cells (zeros on the diagonal).
	gsl_matrix* contact_matrix = gsl_matrix_calloc(number_of_connected_cells, number_of_connected_cells);
	HelperCheckerboard_BuildAdjacencyMatrix(ordered_next_differentiated_cells, contact_matrix);

	//Vector for storing the different types of the connected cells.
	gsl_vector* cell_type_vectors = gsl_vector_calloc(number_of_connected_cells);
	//J matrix times cell_type_vectors.
	gsl_vector* coupling_vector = gsl_vector_calloc(number_of_connected_cells);

	//Build id2index map for the connected cells.
	unordered_map<long, int> connectedcells_id2index;
	connectedcells_id2index.reserve(number_of_connected_cells);								//reserve to contain vertices.size() elements

	//Store initial random configuration.
	vector<double> min_energy_configuration;
	//Store calculated energy with the given configuration.
	//calculated_energy -> energy calculated for a given configuration.
	//current_energy -> store energy of current configuration.
	//min_energy -> store energy of minimum energy (found during the entire search thus far).
	//delta_energy -> store change in energy between current configuration and proposed change.
	double calculated_energy, current_energy, min_energy, delta_energy, delta_min_energy;

	//Store a history of past energy in order to see if the system has equilibrized during single-spin-flips.
	vector<double> past_energy_vector;
	//randomly chosen index for spin flip.
	int random_spin_index;
	//define uniform integer distribution
	uniform_int_distribution<int> distribution(0, RAND_MAX);
	
	//calculate the characteristic length (coupling energies are based on contact lengths : only use the area filled by actual cells)
	double area_of_packing = LX * LY;
	list<Cell*>::iterator cit;
	for (cit = (*this).cells.begin(); cit != (*this).cells.end(); cit++){
		//Subtract out energy for inner and outer mega-cell (center_cell_type and outer_cell_type, respectively).
		if ((*(*cit)).type == center_cell_type || (*(*cit)).type == outer_cell_type){
			area_of_packing -= (*(*cit)).GetArea();
		}
	}
	//Store number of cells (subtract out both megacells).
	int number_of_cells = (*this).cells.size() - 2;
	double characteristic_length = sqrt(area_of_packing / number_of_cells / 3 / sqrt(3) * 2);
	double characteristic_energy_scale = characteristic_length;
	//variable for storing standard deviation of the past_energy_vector.
	double history_st_dev = characteristic_energy_scale;
	//store inverse temperature (~inverse of characteristic energy)
	double beta = 2. / characteristic_length;
	//length of history to be stored: sample time scale of approximately ~10 fluctuations from a metastable state.
	int length_of_history = (int)((1 + exp(-4. * beta * characteristic_energy_scale)) / exp(-4. * beta * characteristic_energy_scale));

	//Store map for id 2 index of the coupling matrix.
	int ind = 0;
	for (int connected_cell_index = 0; connected_cell_index != number_of_connected_cells; connected_cell_index++){
		connectedcells_id2index[(*connected_cells.at(connected_cell_index)).id] = ind;
		ind++;
	}

	//For cells that already have PCP that are connected to the next differentiated cells, store the corresponding spin (type 1 <-> -1.0, type 2 <-> 1.0)
	for (int cell_index = 0; cell_index != number_of_connected_cells; cell_index++){
		if ((*connected_cells.at(cell_index)).type == 1){
			gsl_vector_set(cell_type_vectors, cell_index, -1.0);
		}
		else if ((*connected_cells.at(cell_index)).type == 2){
			gsl_vector_set(cell_type_vectors, cell_index, 1.0);
		}
	}

	//Store the "spins" from randomly generated initial distribution.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
		gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(cell_index));
	}
	//Multiply the matrix of couplings by the spin vector (store result in coupling vector).
	gsl_blas_dgemv(CblasNoTrans, 1.0, contact_matrix, cell_type_vectors, 0.0, coupling_vector);
	//Do the dot product of coupling vector and spin vector -> store result in calculated energy.
	gsl_blas_ddot(cell_type_vectors, coupling_vector, &calculated_energy);
	//Store the initial energy as initial min energy, and initialize min_energy_configuration as initial "spin distribution"
	min_energy = calculated_energy;
	current_energy = min_energy;
	past_energy_vector.push_back(current_energy);
	min_energy_configuration = next_diff_cells_types;

	//Output initial type ("spin") distribution.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		if (min_energy_configuration.at(cell_index) == -1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 1;
		}
		else if (min_energy_configuration.at(cell_index) == 1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 2;
		}
	}
	(*this).WriteOutputDataAppend(datafile);

	//Keep track the number of attempts of spin flips (stop when you both satisfy the equilibration of spins and have at least gone length_of_history steps).
	int number_of_attempts = 0;
	do{
		//Select a random index for flipping the corresponding spin.
		random_spin_index = distribution(generator) % number_of_next_differentiated_cells;
		//Find absolute cell index -> in the convention of the coupling matrix.
		int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(random_spin_index)).id];
		//Flip that spin in the cell_type_vectors.
		gsl_vector_set(cell_type_vectors, absolute_cell_index, -1. * next_diff_cells_types.at(random_spin_index));
		//Multiply the matrix of couplings by the spin vector (store result in coupling vector).
		gsl_blas_dgemv(CblasNoTrans, 1.0, contact_matrix, cell_type_vectors, 0.0, coupling_vector);
		//Do the dot product of coupling vector and spin vector -> store result in calculated energy.
		gsl_blas_ddot(cell_type_vectors, coupling_vector, &calculated_energy);
		//Calculate the change in energy brought about by this spin flip (relative to the current configuration and the minimum energy configuration thus far).
		delta_energy = calculated_energy - current_energy;
		delta_min_energy = calculated_energy - min_energy;
		//Calculate the probability of the spin flip (based on delta_energy with Boltzmann weighting).
		double probability_flip = exp(-beta * delta_energy) / (1 + exp(-beta * delta_energy));
		//If this is currently the min energy state, store the configuration and the energy.
		if (delta_energy < 0.0){
			//Flip the corresponding spin.
			next_diff_cells_types.at(random_spin_index) = -1. * next_diff_cells_types.at(random_spin_index);
			//Store the calculated energy as the new minimum, and store the configuration.
			current_energy = calculated_energy;
		}
		//If the change in energy is positive, use a boltzmann weighting to judge if the change should proceed.
		else if (delta_energy > 0.0){
			bernoulli_distribution bern_distribution(probability_flip); 														//distribution for weighting chance of success
			//Use boltzmann weighting.
			if (bern_distribution(generator)){
				//Flip the corresponding spin.
				next_diff_cells_types.at(random_spin_index) = -1. * next_diff_cells_types.at(random_spin_index);
				//Store this energy as the current energy (but not the minimum).
				current_energy = calculated_energy;
			}
			else{
				//Change the spin in the cell type vectors back to the original state.
				gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(random_spin_index));
			}
		}
		//Increment the index for the number of attempts.
		number_of_attempts++;

		//If this configuration's energy is less than the currently stored minimum, update the minimum.
		if (delta_min_energy < 0.0){
			min_energy = calculated_energy;
			min_energy_configuration = next_diff_cells_types;
		}

		//If the number of attempts is less than the length of the history, push back the energy.
		if (number_of_attempts < length_of_history){
			past_energy_vector.push_back(current_energy);
		}
		//Else, update the energy vector on a rotating basis (dropping off oldest stored energy, adding new energy).
		else{
			vector<double> updated_energy_vector;
			for (int i = 1; i != length_of_history; i++) updated_energy_vector.push_back(past_energy_vector.at(i));
			updated_energy_vector.push_back(current_energy);
			past_energy_vector = updated_energy_vector;
			//Store the standard deviation of the history.
			history_st_dev = sqrt(variance_history(past_energy_vector));
		}
	} while (history_st_dev > characteristic_energy_scale || number_of_attempts < length_of_history + 1);

	//Output the configuration after this individual spin flip step.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		if (min_energy_configuration.at(cell_index) == -1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 1;
		}
		else if (min_energy_configuration.at(cell_index) == 1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 2;
		}
	}
	(*this).WriteOutputDataAppend(datafile);

	//Find domains of spin, attempt to flip domains: start with min energy configuration.
	next_diff_cells_types = min_energy_configuration;
	current_energy = min_energy;
	//Store this configuration in the cell_type_vectors.
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
		gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(cell_index));
	}

	//Find mismatches in cell types in the chain of cells.
	vector<int> mismatched_cells = find_domains_of_spin(next_diff_cells_types);
	//Only attempt the number of domain flips consistent with the number of mismatches.
	int number_of_domain_flip_attempts = mismatched_cells.size();

	//Keep track the number of attempts of spin domain flips (stop when you reach number_of_domain_flip_attempts).
	int spin_flip_attempt = 0;
	while (spin_flip_attempt != number_of_domain_flip_attempts && mismatched_cells.size() != 0){
		//Find mismatches in cell types in the chain of cells.
		mismatched_cells = find_domains_of_spin(next_diff_cells_types);
		//Select a domain randomly by randomly selecting the corresponding mismatch.
		int random_mismatch_index = distribution(generator) % mismatched_cells.size();

		//Start at the domain mismatch.
		int cell_index = mismatched_cells.at(random_mismatch_index);
		int mod_cell_index;
		//Loop through flipping the full domain of spins until the next mismatch.
		do{
			//As you increment the cell index, use mod to wrap through the list of cells (the cells are in a ring).
			mod_cell_index = cell_index % number_of_next_differentiated_cells;
			//Find absolute cell index -> in the convention of the coupling matrix.
			int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(mod_cell_index)).id];
			//Flip the spin in the cell_type_vectors.
			gsl_vector_set(cell_type_vectors, absolute_cell_index, -1 * next_diff_cells_types.at(mod_cell_index));
			//Increment cell index.
			cell_index++;
			mod_cell_index = cell_index % number_of_next_differentiated_cells;
		} while (mod_cell_index != mismatched_cells.at((random_mismatch_index + 1) % mismatched_cells.size()));
		gsl_blas_dgemv(CblasNoTrans, 1.0, contact_matrix, cell_type_vectors, 0.0, coupling_vector);
		//Do the dot product of coupling vector and spin vector -> store result in calculated energy.
		gsl_blas_ddot(cell_type_vectors, coupling_vector, &calculated_energy);
		//Calculate the change in energy brought about by this spin flip (relative to the current configuration and the minimum energy configuration thus far).
		delta_energy = calculated_energy - current_energy;
		delta_min_energy = calculated_energy - min_energy;
		//Calculate the probability of the spin domain flip (based on delta_energy with Boltzmann weighting).
		double probability_flip = exp(-beta * delta_energy) / (1 + exp(-beta * delta_energy));
		//Increment spin_flip_attempt index.
		spin_flip_attempt++;

		//If this is currently the min energy state, store the configuration and the energy.
		if (delta_energy < 0.0){
			//Update the vector of next_diff_cells.
			for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
				int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
				next_diff_cells_types.at(cell_index) = gsl_vector_get(cell_type_vectors, absolute_cell_index);
			}
			//Store this energy as the current energy (but not the minimum).
			current_energy = calculated_energy;
		}
		//If the change in energy is positive, use a boltzmann weighting to judge if the change should proceed.
		else if (delta_energy > 0.0){
			bernoulli_distribution bern_distribution(probability_flip); 																//distribution for weighting chance of success
			//Use boltzmann weighting.
			if (bern_distribution(generator)){
				//Flip the corresponding domain of spins.
				for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
					int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(cell_index)).id];
					next_diff_cells_types.at(cell_index) = gsl_vector_get(cell_type_vectors, absolute_cell_index);
				}
				//Update current energy.
				current_energy = calculated_energy;
			}
			else{
				//If not, reset the cell_type_vectors.
				cell_index = mismatched_cells.at(random_mismatch_index);
				do{
					mod_cell_index = cell_index % number_of_next_differentiated_cells;
					int absolute_cell_index = connectedcells_id2index[(*ordered_next_differentiated_cells.at(mod_cell_index)).id];
					gsl_vector_set(cell_type_vectors, absolute_cell_index, next_diff_cells_types.at(mod_cell_index));
					cell_index++;
					mod_cell_index = cell_index % number_of_next_differentiated_cells;
				} while (mod_cell_index != mismatched_cells.at((random_mismatch_index + 1) % mismatched_cells.size()));
			}
		}

		//If this configuration's energy is less than the currently stored minimum, update the minimum.
		if (delta_min_energy < 0.0){
			min_energy = calculated_energy;
			min_energy_configuration = next_diff_cells_types;
		}
		mismatched_cells = find_domains_of_spin(next_diff_cells_types);
	}

	//Now differentiate the cells (based on minimum energy configuration).
	for (int cell_index = 0; cell_index != ordered_next_differentiated_cells.size(); cell_index++){
		if (min_energy_configuration.at(cell_index) == -1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 1;
		}
		else if (min_energy_configuration.at(cell_index) == 1.){
			(*ordered_next_differentiated_cells[cell_index]).type = 2;
		}
	}
	(*this).WriteOutputDataAppend(datafile);

	//Free the gsl allocations.
	gsl_matrix_free(contact_matrix);
	gsl_vector_free(coupling_vector);
	gsl_vector_free(cell_type_vectors);
}

//! Order cell pointer by the angle that the cell centroid about the central mega cell.
vector<PCPCell*> PCPMegaTissue::order_cells_angle(const vector<PCPCell*> &unordered_cells){
	//Store the number of cells in the list.
	int number_of_cells = unordered_cells.size();
	//Build a vector for storing angles of cells relative to central mega cell.
	vector<double> angle_vector;
	angle_vector.resize(number_of_cells);
	//Build a vector for storing ordered cells (pointers).
	vector<PCPCell*> ordered_cells;
	ordered_cells.resize(number_of_cells);
	//Store a list of angle -> cell.
	unordered_map<double, PCPCell*> cell_angle2pointer;
	cell_angle2pointer.reserve(number_of_cells);							//reserve to contain vertices.size() elements

	//Store the origin of the central cell.
	Vect origin;
	//Initialize all cell types to zero (verify).
	list<Cell*>::iterator cit;
	for (cit = (*this).cells.begin(); cit != (*this).cells.end(); cit++){
		if ((*(*cit)).type == center_cell_type){
			origin = (*(*cit)).GetCenter();
		}
	}

	//Store the angle for each cell in the unordered_cells list.
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		Vect vect_to_origin = (*unordered_cells.at(cell_index)).GetCenter() - origin;
		angle_vector.at(cell_index) = vect_to_origin.Angle();
	}

	//Build map of angle->cell.
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		cell_angle2pointer[angle_vector.at(cell_index)] = unordered_cells.at(cell_index);
	}

	//Sort the pointers by angle.
	std::sort(angle_vector.begin(), angle_vector.end());

	//Store the ordered cells (by angle).
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		ordered_cells.at(cell_index) = cell_angle2pointer[angle_vector.at(cell_index)];
	}

	return ordered_cells;
}

//! Output vector of angles about the central mega cell.
// \note The vector obeys the conventions of order_cell_angles.
vector<double> PCPMegaTissue::generate_cell_angle_vector(const vector<PCPCell*> &unordered_cells){
	//Store the number of cells in the list.
	int number_of_cells = unordered_cells.size();
	//Build a vector for storing angles of cells relative to central mega cell.
	vector<double> angle_vector;
	angle_vector.resize(number_of_cells);
	//Build a vector for storing ordered cells (pointers).
	vector<PCPCell*> ordered_cells;
	ordered_cells.resize(number_of_cells);
	//Store a list of angle -> cell.
	unordered_map<double, PCPCell*> cell_angle2pointer;
	cell_angle2pointer.reserve(number_of_cells);							//reserve to contain vertices.size() elements

	//Store the angle for each cell in the unordered_cells list.
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		PCPCell* pcp_cell_ptr = dynamic_cast<PCPCell*>(unordered_cells.at(cell_index));
		angle_vector.at(cell_index) = (*this).Helper_Find_Cell_Angle(pcp_cell_ptr);
	}

	//Build map of angle->cell.
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		cell_angle2pointer[angle_vector.at(cell_index)] = unordered_cells.at(cell_index);
	}

	//Sort the pointers by angle.
	std::sort(angle_vector.begin(), angle_vector.end());

	//Store the ordered cells (by angle).
	for (int cell_index = 0; cell_index != number_of_cells; cell_index++){
		ordered_cells.at(cell_index) = cell_angle2pointer[angle_vector.at(cell_index)];
	}

	return angle_vector;
}

//! Find angle of cell about central megacell.
double PCPMegaTissue::Helper_Find_Cell_Angle(PCPCell* pcp_cell_ptr){
	double angle_calc;
	//Store the origin of the central cell.
	Vect origin;
	//Initialize all cell types to zero (verify).
	list<Cell*>::iterator cit;
	for (cit = (*this).cells.begin(); cit != (*this).cells.end(); cit++){
		if ((*(*cit)).type == center_cell_type){
			origin = (*(*cit)).GetCenter();
		}
	}

	Vect vect_to_origin = (*pcp_cell_ptr).GetCenter() - origin;
	angle_calc = vect_to_origin.Angle();
	return angle_calc;
}

//! Make sure that you propogate PCP outward from origin_mega_cell.
Edge* PCPMegaTissue::edge_candidate_outward(const vector<Edge*> &edge_candidates){
	//Store all edges in edge_candidates that form the right boundary for PCP propogation (make sure that front propogation is consistently in one direction).
	vector<Edge*> edges_to_right;

	//Loop through edge_candidates.
	for (int edge_index = 0; edge_index != edge_candidates.size(); edge_index++){
		//Define pointer to cell 1 and cell2.
		Cell* cell1ptr = (*edge_candidates.at(edge_index)).cell1;
		Cell* cell2ptr = (*edge_candidates.at(edge_index)).cell2;
		Cell* type0_cellptr;

		//Figure out which cell is of type 0.
		if ((*cell1ptr).type == 0){
			type0_cellptr = cell1ptr;
		}
		else if ((*cell2ptr).type == 0){
			type0_cellptr = cell2ptr;
		}

		//Returns a vector of pointers to all vertices in the cell in clockwise (CW) order: for the cell of type 0.
		vector<Vertex*> cw_verts = (*type0_cellptr).GetClockWiseVertices();
		//Store corresponding ordered vert positions (respecting PBC).
		vector<vector<double>> ordered_vert_positions = (*type0_cellptr).HelperOrderVerticesPos(cw_verts);
		//Define the geocenter of the ordered_vert_positions.
		Vect position = GetGeoCenter(ordered_vert_positions[0], ordered_vert_positions[1]);
		//Define vector from origin_mega_cell to position.
		Vect vector_from_center = position - origin_mega_cell;
		//Store a vector of indices for edges that a ray pointing to the right intersects.
		vector<int> intersection_indices = which_Edge_Intersection(position, ordered_vert_positions, vector_from_center.Angle());
		//For a convex polygon, there should only one intersection.
		assert((intersection_indices.size() == 1) && "ERROR: COULD NOT DETERMINE DIRECTION FOR FRONT PROPAGATION");
		int intersection_index = intersection_indices.at(0);
		//Store the number of vertices in this cell (type0_cell).
		int number_of_vertices = cw_verts.size();
		//Find next index in the list of vertices (wrapping to beginning of vector by mod).
		int next_index = (intersection_index + 1) % number_of_vertices;
		//Find the edge between the two vertices (the intersected edge).
		Edge* edge_between = (*cw_verts.at(intersection_index)).GetEdgeBetween(cw_verts.at(next_index));
		//Find the cell on the other side of this edge.
		Cell* other_cell = (*edge_between).GetOtherCell(type0_cellptr);
		//If the other cell is of type 0, edge_candidates.at(edge_index) can serve for propogation of PCP to the right.
		if ((*other_cell).type == 0){
			edges_to_right.push_back(edge_candidates.at(edge_index));
		}
	}

	//Assert that you have found at least one edge candidate for propogating PCP to the right.
	assert((edges_to_right.size() != 0) && "ERROR: COULD NOT DETERMINE DIRECTION FOR FRONT PROPAGATION");

	//Return only one of the edges.
	Edge* chosen_edge_candidate = edges_to_right.at(0);

	return chosen_edge_candidate;
}

//Check fourfold vertex for stability (make sure that this function creates appropriate types of objects -> PCPEdge and PCPMegaCell)
double PCPMegaTissue::HelperResolveVertices_CheckStability(Vertex* vert, Edge* edge1, Edge* edge2, Edge* edge3, Edge* edge4, Cell* cell1, Cell* cell2, Cell* cell3, Cell* cell4, Vect *Force){
	//First make new objects which are coppies of the objects in the linked lists we will minipulate these objecs and they will be destroyed when the function ends
	PCPMegaCell* pcp_c1 = dynamic_cast<PCPMegaCell*>(cell1);
	PCPMegaCell* pcp_c2 = dynamic_cast<PCPMegaCell*>(cell2);
	PCPMegaCell* pcp_c3 = dynamic_cast<PCPMegaCell*>(cell3);
	PCPMegaCell* pcp_c4 = dynamic_cast<PCPMegaCell*>(cell4);
	PCPMegaCell c1 = (*pcp_c1);
	PCPMegaCell c2 = (*pcp_c2);
	PCPMegaCell c3 = (*pcp_c3);
	PCPMegaCell c4 = (*pcp_c4);
	PCPEdge* pcp_edge1 = dynamic_cast<PCPEdge*>(edge1);
	PCPEdge* pcp_edge2 = dynamic_cast<PCPEdge*>(edge2);
	PCPEdge* pcp_edge3 = dynamic_cast<PCPEdge*>(edge3);
	PCPEdge* pcp_edge4 = dynamic_cast<PCPEdge*>(edge4);
	PCPEdge e1 = (*pcp_edge1);
	PCPEdge e2 = (*pcp_edge2);
	PCPEdge e3 = (*pcp_edge3);
	PCPEdge e4 = (*pcp_edge4);
	int stability;
	//Changes by Hayden for PCP.
	int largestID = (*max_element(edgeIDs.begin(), edgeIDs.end()));
	double a1, a2, a3, a4;

	//Make sure that these edges have distinct ids
	Vertex v12; //the vertex on edges 1 and 2
	Vertex v34; //the vertex on edges 3 and 4

	//Note: important change, initilizing all at once (also assigning an actual id rather than having garbage handed to it).
	//Make the new edge and vertices resulting from the split into 2 3fold vertices
	//Please note: enew has no polarity protein.
	PCPEdge enew = PCPEdge::PCPEdge(largestID + 1, 0, 0, &v12, &v34, &c4, &c2);

	//Now update the info on the coppies of the objects so they are 2 3fold vertices

	//For the new vertices
	v12.r = (*vert).r;							//put the two new verts where the old fourfold was
	v12.edges.push_back(&e1);
	v12.edges.push_back(&e2);
	v12.edges.push_back(&enew);
	v34.r = (*vert).r;							//put the two new verts where the old fourfold was
	v34.edges.push_back(&e3);
	v34.edges.push_back(&e4);
	v34.edges.push_back(&enew);

	//For the old edges
	e1.UpdateVertex(vert, &v12);
	e2.UpdateVertex(vert, &v12);
	e3.UpdateVertex(vert, &v34);
	e4.UpdateVertex(vert, &v34);

	//Now we have at the position of the fourfold vertex two threefold vertices in which everyone knows all of their neighbors correclty, we can now calculate the forces based on what the tensiions and pressures are in this configuration

	//pressure differences across cells
	a1 = c1.GetArea();
	a2 = c2.GetArea();
	a3 = c3.GetArea();
	a4 = c4.GetArea();
	double p2 = (c2.GetPressure(a2) - c1.GetPressure(a1)); 			//pressure difference acorss edge2 between cell 1 and cell 2 taken CCW
	double p3 = (c3.GetPressure(a3) - c2.GetPressure(a2));
	double p4 = (c4.GetPressure(a4) - c3.GetPressure(a3));
	double p1 = (c1.GetPressure(a1) - c4.GetPressure(a4));

	//getting z X edge (both the diredction and magnitude are correct
	Vect zXe1 = (e1).GetVectorDifference(&v12, v12.GetOtherVertex(&e1)).Rotate(PI / 2.0);
	Vect zXe2 = (e2).GetVectorDifference(&v12, v12.GetOtherVertex(&e2)).Rotate(PI / 2.0);
	Vect zXe3 = (e3).GetVectorDifference(&v34, v34.GetOtherVertex(&e3)).Rotate(PI / 2.0);
	Vect zXe4 = (e4).GetVectorDifference(&v34, v34.GetOtherVertex(&e4)).Rotate(PI / 2.0);

	//Get the forces
	Vect ForceT1 = (e1).GetTension_ForceOnVertex(&v12);
	Vect ForceT2 = (e2).GetTension_ForceOnVertex(&v12);
	Vect ForceT3 = (e3).GetTension_ForceOnVertex(&v34);
	Vect ForceT4 = (e4).GetTension_ForceOnVertex(&v34);
	Vect ForceP1 = zXe1 * (p1 * 0.5);
	Vect ForceP2 = zXe2 * (p2 * 0.5);
	Vect ForceP3 = zXe3 * (p3 * 0.5);
	Vect ForceP4 = zXe4 * (p4 * 0.5);
	*Force = ((ForceT1 + ForceT2 + ForceP1 + ForceP2) - (ForceT3 + ForceT4 + ForceP3 + ForceP4))* 0.5;

	Vect ForceNewEdge;
	//c4 -> cell1 of enew
	//c2 -> cell2 of enew
	//Check which edges belong to c4 and which belong to c2.
	//e1 and e2 are paired, and e3 and e4 are paired.

	//Calculate TENSION ON VIRTUAL EDGE.
	if (((c4 == *e1.cell1) || (c4 == *e1.cell2)) && ((c4 == *e3.cell1) || (c4 == *e3.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e3, &e2, &e4);
	}
	else if (((c4 == *e1.cell1) || (c4 == *e1.cell2)) && ((c4 == *e4.cell1) || (c4 == *e4.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e4, &e2, &e3);
	}
	else if (((c2 == *e1.cell1) || (c2 == *e1.cell2)) && ((c2 == *e3.cell1) || (c2 == *e3.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e3, &e2, &e4);
	}
	else if (((c2 == *e1.cell1) || (c2 == *e1.cell2)) && ((c2 == *e4.cell1) || (c2 == *e4.cell2))){
		ForceNewEdge = enew.GetTension_ForceOnVertex_Virtual(&v12, &e1, &e4, &e2, &e3);
	}

	//Check the stability
	stability = (*Force).Magnitude() - ForceNewEdge.Magnitude();

	return stability;
}