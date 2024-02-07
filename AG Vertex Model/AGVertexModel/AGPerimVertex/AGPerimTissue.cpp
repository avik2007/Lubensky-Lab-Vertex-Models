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
#include <array>
#include <typeinfo>
#include <numeric>
#include <stddef.h>
#include <iomanip>
#include <cmath>
#include <complex> 
#include <Eigen/Dense>
#include <unordered_map>

#include "Constants.h"
#include "AGConstants.h"
#include "TissueClass.h"
#include "AGTissue.h"
#include "AGCell.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "AGPerimTissue.h"
#include "AGPerimConstants.h"
#include "AGPerimCell.h"
#include "AGPerimEdge.h"

using namespace Eigen;

/*
* void HelperReadInputData_SetGlobalValues(ifstream* datafile)
*
* Reads the data and sets all of the global variables.
* Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
*/
void AGPerimTissue::HelperReadInputData_SetGlobalValues(ifstream* datafile){

	string line, word, key, value;		//temporary variables to store lines and words 
	string type;						//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline;				//ssline is a string stream so we can comvert strings to stringstreams
	map <string, string> globals;



	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//

	//Read until you get to the word Global signifying the begining of the global variables
	while (line.compare("Global:") != 0){
		getline(*datafile, line);
		cout << line << endl;
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

	if (globals.count("ALPHAM") == 1){ alphaM = atof(globals["ALPHAM"].c_str()); }
	else if (globals.count("ALPHAM") == 0){ cout << "WARNING: ALPHAM is not given in input file, used default value " << alphaM << endl; }
	else { cout << "ERROR: Global variable ALPHAM has multiple definition in input file, used default value ALPHAM = " << alphaM << endl; }

	if (globals.count("TISSUE_ALPHA") == 1){ TISSUE_ALPHA = atof(globals["TISSUE_ALPHA"].c_str()); }
	else if (globals.count("TISSUE_ALPHA") == 0){ cout << "WARNING: TISSUE_ALPHA is not given in input file, used default value " << TISSUE_ALPHA << endl; }
	else { cout << "ERROR: Global variable TISSUE_ALPHA has multiple definition in input file, used default value TISSUE_ALPHA = " << TISSUE_ALPHA << endl; }

	if (globals.count("TISSUE_EXT_STRESS") == 1){ TISSUE_EXT_STRESS = atof(globals["TISSUE_EXT_STRESS"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS") == 0){ cout << "WARNING: TISSUE_EXT_STRESS is not given in input file, used default value " << TISSUE_EXT_STRESS << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS has multiple definition in input file, used default value TISSUE_EXT_STRESS = " << TISSUE_EXT_STRESS << endl; }

	if (globals.count("TISSUE_EXT_STRESS_X") == 1){ TISSUE_EXT_STRESS_X = atof(globals["TISSUE_EXT_STRESS_X"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_X") == 0){ cout << "WARNING: TISSUE_EXT_STRESS_X is not given in input file, used default value " << TISSUE_EXT_STRESS_X << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_X has multiple definition in input file, used default value TISSUE_EXT_STRESS_X = " << TISSUE_EXT_STRESS_X << endl; }

	if (globals.count("TISSUE_EXT_STRESS_Y") == 1){ TISSUE_EXT_STRESS_Y = atof(globals["TISSUE_EXT_STRESS_Y"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_Y") == 0){ cout << "WARNING: TISSUE_EXT_STRESS_Y is not given in input file, used default value " << TISSUE_EXT_STRESS_Y << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_Y has multiple definition in input file, used default value TISSUE_EXT_STRESS_Y = " << TISSUE_EXT_STRESS_Y << endl; }

	if (globals.count("MITOSISON") == 1){ MITOSISON = (atof(globals["MITOSISON"].c_str()) > 0); }
	else if (globals.count("MITOSISON") == 0){ cout << "WARNING: MITOSISON is not given in input file, used default value " << MITOSISON << endl; }
	else { cout << "ERROR: Global variable MITOSISON has multiple definition in input file, used default value MITOSISON = " << MITOSISON << endl; }

	if (globals.count("T1ON") == 1){ T1ON = (atof(globals["T1ON"].c_str()) > 0); }
	else if (globals.count("T1ON") == 0){ cout << "WARNING: T1ON is not given in input file, used default value " << T1ON << endl; }
	else { cout << "ERROR: Global variable T1ON has multiple definition in input file, used default value T1ON = " << T1ON << endl; }

	if (globals.count("APOPTOSISON") == 1){ APOPTOSISON = (atof(globals["APOPTOSISON"].c_str()) > 0); }
	else if (globals.count("APOPTOSISON") == 0){ cout << "WARNING: APOPTOSISON is not given in input file, used default value " << APOPTOSISON << endl; }
	else { cout << "ERROR: Global variable APOPTOSISON has multiple definition in input file, used default value APOPTOSISON = " << APOPTOSISON << endl; }

	if (globals.count("WRITEMITOSIS") == 1){ WRITEMITOSIS = atof(globals["WRITEMITOSIS"].c_str()); }
	else if (globals.count("WRITEMITOSIS") == 0){ cout << "WARNING: WRITEMITOSIS is not given in input file, used default value " << WRITEMITOSIS << endl; }
	else { cout << "ERROR: Global variable WRITEMITOSIS has multiple definition in input file, used default value WRITEMITSOSIS = " << WRITEMITOSIS << endl; }

	if (globals.count("VERBOSEON") == 1){ VERBOSEON = (atof(globals["VERBOSEON"].c_str()) > 0); }
	else if (globals.count("VERBOSEON") == 0){ cout << "WARNING: VERBOSEON is not given in input file, used default value " << VERBOSEON << endl; }
	else { cout << "ERROR: Global variable VERBOSEON has multiple definition in input file, used default value APOPTOSISON = " << VERBOSEON << endl; }

	if (globals.count("EQUILIBRIUMMITOSIS") == 1){ EQUILIBRIUMMITOSIS = (atof(globals["EQUILIBRIUMMITOSIS"].c_str()) > 0); }
	else if (globals.count("EQUILIBRIUMMITOSIS") == 0){ cout << "WARNING: EQUILIBRIUMMITOSIS is not given in input file, used default value " << EQUILIBRIUMMITOSIS << endl; }
	else { cout << "ERROR: Global variable EQUILIBRIUMMITOSIS has multiple definition in input file, used default value EQUILIBRIUMMITOSIS = " << EQUILIBRIUMMITOSIS << endl; }

	if (globals.count("P0") == 1){ P0 = (atof(globals["P0"].c_str())); }
	else if (globals.count("P0") == 0){ cout << "WARNING: P0 is not given in input file, used default value " << P0 << endl; }
	else { cout << "ERROR: Global variable P0 has multiple definition in input file, used default value P0 = " << P0 << endl; }

	if (globals.count("LAMBDA") == 1){ LAMBDA = (atof(globals["LAMBDA"].c_str())); }
	else if (globals.count("LAMBDA") == 0){ cout << "WARNING: LAMBDA is not given in input file, used default value " << LAMBDA << endl; }
	else { cout << "ERROR: Global variable LAMBDA has multiple definition in input file, used default value LAMBDA = " << LAMBDA << endl; }

	if (globals.count("A0") == 1){ A0 = (atof(globals["A0"].c_str())); }
	else if (globals.count("A0") == 0){ cout << "WARNING: A0 is not given in input file, used default value " << A0 << endl; }
	else { cout << "ERROR: Global variable A0 has multiple definition in input file, used default value A0 = " << A0 << endl; }

	if (globals.count("energyListSize") == 1){ energyListSize = (atof(globals["energyListSize"].c_str())); }
	else if (globals.count("energyListSize") == 0){ cout << "WARNING: energyListSize is not given in input file, used default value " << energyListSize << endl; }
	else { cout << "ERROR: Global variable energyListSize has multiple definition in input file, used default value energyListSize = " << energyListSize << endl; }

	if (globals.count("equilibriumtolerance") == 1){ equilibriumtolerance = (atof(globals["equilibriumtolerance"].c_str())); }
	else if (globals.count("equilibriumtolerance") == 0){ cout << "WARNING: equilibriumtolerance is not given in input file, used default value " << equilibriumtolerance << endl; }
	else { cout << "ERROR: Global variable equilibriumtolerance has multiple definition in input file, used default value equilibriumtolerance = " << equilibriumtolerance << endl; }

	if (globals.count("C4") == 1) { C4 = atof(globals["C4"].c_str()); }
	else if (globals.count("C4") == 0) { cout << "WARNING: C4 is not given in input file, used default value " << C4 << endl; }
	else { cout << "ERROR: Global variable C4 has multiple definition in input file, used default value C4 = " << C4 << endl; }
	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//
	int i = 1;
	if (energyListSize < 3) { energyListSize = 3; } // minimal size for this has to be 2, I'm setting it to 3 here since if it is equal to 2 there would be no coarse grained slope analysis
	// setting set the initial box dimensions here
	LX0 = LX;
	LY0 = LY;

	// setting the initial stepsize upperbound here
	//(*this).stepsize_upperbound = 1; 
	// nominally, there's definitely a better place to set this - consider upgrading this at some pointd
}

/*
* void HelperReadInputData_CreateObjects
*
* Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's
* list of ids.
* Starts at current position of *datafile and reads until keyword 'End'.
*/
void AGPerimTissue::HelperReadInputData_CreateObjects(ifstream* datafile){

	string line, key, value, edg;			//temporary variables to store lines and words 
	string type;							//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline, ssedges;			//ssline is a string stream so we can comvert strings to stringstreams
	vector<string> keywords;				//vector of strings to hold the keywords
	vector<int> listedges;
	int i, j;
	map <string, string> objectInfo;		//map to hold the keywords and their data
	int newid, v1, v2, enumb, typenumb;					//temp varables to store the object info
	double x, y;
	bool xf, yf;
	double A0;
	string potential, cellpressuretype;
	double tcjstrength = 0;
	double thetaan = 0;
	double anstrength = 0;

	cout << "Making a base class tissue" << endl;

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
		getline(*datafile, line);							//get the next line
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
		newid = atoi(objectInfo["id"].c_str());						//extract the data from the map
		x = atof(objectInfo["x"].c_str());
		y = atof(objectInfo["y"].c_str());
		Vertex* vptr = new Vertex;									// Create a vertes with an id and x,y pos and put it in the master list of all veritces
		(*vptr).id = newid;
		(*vptr).r.x = x;
		(*vptr).r.y = y;
		(vertices).push_back(vptr);
		vertexIDs.push_back(newid);								//and add its id to the list of ids
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
		getline(*datafile, line);							//get the next line
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
		newid = atoi(objectInfo["id"].c_str());						//extract the data from the map
		v1 = atoi(objectInfo["vertex1"].c_str());
		v2 = atoi(objectInfo["vertex2"].c_str());
		xf = atoi(objectInfo["xflag"].c_str());
		yf = atoi(objectInfo["yflag"].c_str());
		AGPerimEdge* eptr = new AGPerimEdge;
		(*eptr).id = newid;
		(*eptr).IOVertexID[0] = v1;
		(*eptr).IOVertexID[1] = v2;
		(*eptr).xflag = xf;
		(*eptr).yflag = yf;
		(*eptr).cell1 = NULL;
		(*eptr).cell2 = NULL;
		(edges).push_back(eptr);		// Create a vertes with an id and x,y pos and put it in the master list of all veritces
		edgeIDs.push_back(newid);								//and add its id to the list of ids
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
		getline(*datafile, line);							//get the next line
		if (line.compare("End") == 0 || line.compare("Fibers") == 0){ break; }				//break if keyworkd 'End' is found
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
		else{ typenumb = 0; }
		if (objectInfo.count("A0")){
			A0 = atof(objectInfo["A0"].c_str());
		}
		// you only want atoi if you have a number you want to store as a string. but we should probably have a string


		ssedges.clear();
		ssedges << objectInfo["EdgeList"];
		//now worry about the edge list
		listedges.clear();
		while (getline(ssedges, edg, ',')){
			listedges.push_back(atoi(edg.c_str()));
		}

		if (objectInfo.count("Potential"))
		{
			potential = objectInfo["Potential"];
		}
		else
		{
			potential = "uniform";
		}

		if (objectInfo.count("TCJStrength"))
		{
			tcjstrength = atof(objectInfo["TCJStrength"].c_str());
		}

		if (objectInfo.count("ThetaAn"))
		{
			thetaan = atof(objectInfo["ThetaAn"].c_str());
		}

		if (objectInfo.count("AnStrength"))
		{
			anstrength = atof(objectInfo["AnStrength"].c_str());
		}
		if (objectInfo.count("CellPressureType"))
		{
			cellpressuretype = objectInfo["CellPressureType"];
		}
		else
		{
			cellpressuretype = "QUADRATIC";
		}
		//make the cell object
		AGPerimCell* cptr = new AGPerimCell;
		(*cptr).id = newid;
		(*cptr).type = typenumb;
		(*cptr).A0 = A0;
		(*cptr).IOEdgeID = listedges;
		(*cptr).DIVISIONPOTENTIAL = potential;
		(*cptr).P0 = P0;
		(*cptr).currentEnergy = 0;
		(*cptr).pastEnergy = 0;
		(*this).A0 = A0;
		(*this).P0 = P0;
		(*cptr).pressureType = cellpressuretype;

		if (tcjstrength != 0)
		{
			(*cptr).TCJSTRENGTH = tcjstrength;
			// if this has been included in the input file, assign the input file values
		}
		else
		{
			(*cptr).TCJSTRENGTH = TCJSTRENGTH;
		}

		if (anstrength != 0 && thetaan != 0)
		{
			(*cptr).THETAAN = thetaan;
			(*cptr).ANSTRENGTH = anstrength;
			// if these have been included in the input file, then assign input file values 
		}
		else
		{
			(*cptr).THETAAN = THETAAN;
			(*cptr).ANSTRENGTH = ANSTRENGTH;
			// otherwise, use the default values assigned in the Main function (which have been included as failsafe)
		}


		cells.push_back(cptr);
		cellIDs.push_back(newid);


		
		

		// instantiate the free energy to be larger than the actual energy of the packing
		// just so that this doesn't mess with an adaptive stepsize changer (for finding
		// equilibrium)


		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
 
	}
	(*this).currentFreeEnergy = A0*A0*A0*(K * A0 + LAMBDA * P0 * P0) * cells.size(); // this just needs to start out as a bigger number than the actual tissue energy
	(*this).pastFreeEnergy = (*this).currentFreeEnergy;
}


/*
* bool CalculateFreeEnergy()
*
* Calculates the energy of the packing using the Hamiltonian
* in this model:
* H = (sum over cells) [ 0.5 K(A_cells - A_0) + 0.5 LAMBDA (P_cell - P_0)^2]
*/
double AGPerimTissue::CalculateFreeEnergy()
 {
	double carea, cperim;
	double a0, p0;
	double freeenergy = 0;
	
	double perim_difference_sqr = 0;
	// (P_c - P_0)^2
	double area_difference_sqr = 0;
	// (A_c - A_0)^2

	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		carea = ((*agit)).GetArea();
		cperim = ((*agit)).GetPerimeter();
		p0 = (*agit).P0;
		a0 = ((*agit)).A0;
		if ((*agit).pressureType == "LOG")
		{
			area_difference_sqr += -2 * log(1 + carea/a0);
			// where the energy would look like Ep = - K A0^2 Log(1 + A/A0) ([K] = pressure)
		}
		else if ((*agit).pressureType == "EXPONENTIAL")
		{
			area_difference_sqr += 2 *exp(-carea / a0); 
			// where the energy would look like Ep = K A0^2 exp((A/A0-1)) ([K] = pressure)
		}
		else if ((*agit).pressureType == "QUARTIC")
		{
			area_difference_sqr += (pow(carea / a0 - 1, 2) + C4*pow(carea / a0 - 1, 4));
		}
		else if ((*agit).pressureType == "QUARTIC2.0")
		{
			area_difference_sqr +=  (pow(carea / a0 - 1, 2)*(1 - 6 * C4 ) + C4 * pow(carea / a0 - 1, 4));
		}
		else
		{
			area_difference_sqr += a0*pow(carea/a0 - 1, 2);
			// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
		}
		 
		
		perim_difference_sqr += pow( cperim - p0, 2); 
		
		// extra factor of 1/A0 is because K is written in units of pressure instead of E/length^4
	}
	freeenergy =0.5*LAMBDA*perim_difference_sqr + 0.5*K*a0*area_difference_sqr - (*this).CalculateBoxEnergy();


	return freeenergy;
}

/*
* Edge* CreateEdge
*
* Creates a new AGPerimEdge and adds it to the list of edges and returns a pointer
*/
Edge* AGPerimTissue::CreateEdge(Edge* subclass, bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2){
	int largestID, newid;
	Edge* newedge{};
	//Get the id
	largestID = (*max_element(edgeIDs.begin(), edgeIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;								//Now make a new ID that is unique

	if (newid > 1000000000) {
		cout << "ERROR to many edges" << endl;
		return NULL;
	}
	//create the vertex and add the vertex to the lists
	AGPerimEdge* eptr = new AGPerimEdge;
	(*eptr).id = newid;
	(*eptr).xflag = xflg;
	(*eptr).yflag = yflg;
	(*eptr).vertex1 = v1;
	(*eptr).vertex2 = v2;
	(*eptr).cell1 = c1;
	(*eptr).cell2 = c2;
	(edges).push_back(eptr);
	edgeIDs.push_back(newid);

	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){
		if ((*(*eit)).id == newid){
			newedge = &(*(*eit));
			break;
		}
	}
	return newedge;
}
/*
* void Evolve( double interval, double *stepsize, ofstream* errorfile)
*
* Evolves the tissue one time interval taking steps no larger than dt and updates the vertex positions
*/
void AGPerimTissue::Evolve(double interval, double *stepsize, ofstream* errorfile)
{

	//find the end point of this evolve
	double endtime = CURRENTTIME + interval;
	bool reinitialize = true;
	bool  flag1, flag2, flag3, flag4;
	int count = 0;
	double previoustime = 0;
	bool energyDidIncrease = 0;
	int index = 0; 
	vector<Vect> initial_positions;
	int dimension_mechanical = 0; 
	bool breakVar = false;
	double quadraticPressureEnergy = 0;
	double expPressureEnergy = 0;
	double logPressureEnergy = 0;
	double perimEnergy = 0;
	double energy_wo_quartic = 0;
	double cell_type_sum = 0; // 1/24/23, useful for checking if any cells need to divide still
	//double holding = 0;
	

	while (CURRENTTIME < endtime)
	{
		count++;
		//STEP 1: If the topology has changed get a list of the clockwise vertices in each cell
		if (reinitialize){																	//If their has been any topological change OR vertices have crossed a boundary 
			for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
				(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
				cell_type_sum += (*(*cit)).type;
				reinitialize = false;
			}
		}

		//STEP 2: Integrate 
		//Loop through the vertices and record the current position in rold 
		// saving initial conditions
		dimension_mechanical = vertices.size() + 1;
		initial_positions.resize(dimension_mechanical);
	
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++) { (*(*vit)).rold = (*(*vit)).r; }
		(*this).Store_Initial_Condition(initial_positions);
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++)
		{
			AGPerimEdge * agpeit = dynamic_cast<AGPerimEdge*>(*eit);
			(*agpeit).xflagold = (*(*eit)).xflag;
			(*agpeit).yflagold = (*(*eit)).yflag;
		}
		LXold = LX;				//Record the old box size
		LYold = LY;

		previoustime = CURRENTTIME;


		// enter do while loop - this makes sure that the integrator doesn't take too large a step and keeps resulting 
		// in decreases in energy
		do
		{
			//(*this).WriteOutputDataAppend(errorfile);
			//******CHOOSE YOUR INTEGRATOR******
			if (INTEGRATORNAME == "RK4_Adaptive"){
				//4th order Runge-Kutta with adaptive stepsize
				(*this).Integrate_RK4_Adaptive(stepsize, error_tolerance_absolute_mechanical, &CURRENTTIME, endtime, errorfile); //Evolve the vertices and box size dimensions with the given error tolerance (adaptive step size). time is incremented within this function.
			}
			else if (INTEGRATORNAME == "RK4"){
				//4th order Runge-Kutta with fized stepsize
				(*this).Integrate_RK4(*stepsize, errorfile);
				CURRENTTIME += *stepsize;	//Update the current time
			}
			else if (INTEGRATORNAME == "Euiler"){
				//Simple Euiler with Fixed stepsize
				(*this).Integrate_Euiler(*stepsize, errorfile);
				CURRENTTIME += *stepsize;	//Update the current time
			}
			else
			{
				cout << "ERROR invalid integrator name! Using Euiler!";
				(*this).Integrate_Euiler(*stepsize, errorfile);
				CURRENTTIME += *stepsize;	//Update the current time
			}


			//Step 3: Update box and boundary flags
			//(*this).Store_New_Positions(trial_integration_step);

			LX = LXnew;				//Update the current box size 
			LY = LYnew;
			(*this).HelperEvolve_UpdatePosition();  //Update the vertex positions (this automaticaly updates edge flags if a vertex has crossed the boundary)


			// STEP 4: Update Free Energy

			// I want to check free energy here! - this way, I can check if we need to backtrack
			pastFreeEnergy = currentFreeEnergy;
			currentFreeEnergy = (*this).CalculateFreeEnergy();
			energyDidIncrease = (currentFreeEnergy - pastFreeEnergy > 0);


			// STEP 5: Check if we're overshooting equilibrium. If so, backtrack, decrease stepsize and change adaptive 
			// integrator to rk4
			
			if (energyDidIncrease)
			{
				//if (*stepsize < 0.0000001) // this could be turned into a parameter, the minumum stepsize
				if (*stepsize <   0.000000000001)
				{
					breakVar = true;
					//*errorfile << "Stepsize got too small. Break" << endl;
					break;
					// This is a step I take so that the adaptive time stepper doesn't keep pushing the stepsize all the way down the machine precision 
				}
				//*errorfile << "Current FE \t Stepsize\t Energy Change \t Polydispersity \n";
				//*errorfile << std::setprecision(outputPrecision) << currentFreeEnergy  << "\t" << *stepsize << "\t" << currentFreeEnergy - pastFreeEnergy << "\t" << (*this).calculatePolydispersity() << endl;
				//*errorfile << "Energy Increase - change timestep and redo step\n";
				*errorfile << "\n";
				//*errorfile << "--------------------------------------------------------------------------\n";
				
				//*errorfile << "--------------------------------------------------------------------------\n";

				// BACKTRACK STEP 1 ----- update to initial positions
				(*this).Update_Old_Mechanical_Variables(initial_positions);

				//BACKTRACK STEP 2 ------ update to inital flags (seems not to throw errors)
				for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++)
				{
					// reverts flags back to previous state - not sure if this 
					AGPerimEdge * agpeit = dynamic_cast<AGPerimEdge*>(*eit);
					(*(*eit)).xflag = (*agpeit).xflagold;
					(*(*eit)).yflag = (*agpeit).yflagold;
				}
				// BACKTRACK STEP 3 - revert the time
				CURRENTTIME = previoustime;
				// BACTRACK STEP 4 - revert the energy - this is error prone due to positions being wrong
				//holding = (*this).CalculateFreeEnergy();//pastFreeEnergy; // revert current free energy back to where it was before the timestel
				//*errorfile << "\n Does the energy match the pre-backtracker energy? " << holding << " vs the old " << pastFreeEnergy << " and " << (holding - pastFreeEnergy) <<" \n";
				currentFreeEnergy = pastFreeEnergy;//holding;

				//BACKTRACK STEP 5: Decrease stepsize AND change integrator
				*stepsize = *stepsize * 0.75;
				INTEGRATORNAME = "RK4";

				//*errorfile << "--------------------------------------------------------------------------\n";
				//(*this).WriteOutputDataAppend(errorfile);
				//*errorfile << " Stepsize is now:" << *stepsize << endl;
				//*errorfile << "--------------------------------------------------------------------------\n";
			}

		} while (energyDidIncrease);

		double energySlope = getLinearRegressionSlope(timeList, energyList);
	
		//*errorfile << std::setprecision(outputPrecision) << CURRENTTIME << "\t" << *stepsize << "\t" << currentFreeEnergy << "\t" << energySlope << "\t" << (currentFreeEnergy - pastFreeEnergy) / (cells.size() * K * A0) << "\n";

		//quadraticPressureEnergy = (*this).CalculateQuadraticPressureEnergy();
		perimEnergy = (*this).CalculatePerimeterEnergy();
		//energy_wo_quartic = perimEnergy + quadraticPressureEnergy;
		
		//*errorfile << std::setprecision(outputPrecision) << CURRENTTIME << "\t" << currentFreeEnergy << "\t" << energySlope << "\t"<< (*this).calculatePolydispersity() << endl;

		//if ((abs((currentFreeEnergy - pastFreeEnergy) / ((CURRENTTIME - previoustime)*(cells.size() * K * A0))) < equilibriumtolerance) && currentFreeEnergy - pastFreeEnergy < 0)

		//STEP 6: look for topological changes and set flags if any have occured 
		flag1 = false;
		flag2 = false;	
		flag3 = false;
		flag4 = false;
		if (APOPTOSISON) flag1 = (*this).Apoptosis(errorfile);


		if (T1ON)
		{
			//*errorfile << "Free Energy BEFORE T1: " << (*this).CalculateFreeEnergy() << endl;
			flag2 = (*this).T1(errorfile);
		}
		if (flag2) { *errorfile << std::setprecision(outputPrecision) << "Free Energy AFTER T1: " << (*this).CalculateFreeEnergy() << endl; }
		flag3 = (*this).ResolveVertices(errorfile);
		

		
		double timeratio = CURRENTTIME / TAU_STRESS;
		double floortr = floor(CURRENTTIME / TAU_STRESS);
		// This needs a proper calculate equilibrium function
		if (MITOSISON)
		{
			// only does mitosis if packing is equilibrated 
			if (EQUILIBRIUMMITOSIS)
			{
				if ((*this).IsEquilibrated())
				{
					flag4 = (*this).Mitosis();
				}
			}
			else //if (timeratio > 1 && timeratio > floortr && timeratio < 2 * INITIALSTEPSIZE + floortr)
				// Gives two timestep windows for mitosis
			{
				flag4 = (*this).Mitosis();
				cell_type_sum = 0;
			}
			// THIS ACTIVATES MITOSIS!!!			
			//if (WRITEMITOSIS && flag4) (*this).WriteOutputDataAppend(errorfile);
		}
		if (flag1  || flag2 || flag3 || flag4) reinitialize = true;  //Determine if we need to reinitialize (recalculate the cell topology)

		

		

		// STEP 7: Update equilibrium checker. Make sure list of energies being analyzed doesn't get bigger
		// than desired limit. 
		if (energyList.size() > energyListSize)
		{ 
			energyList.pop_front();
			timeList.pop_front();
		}

		energyList.push_back(currentFreeEnergy / (cells.size() * K * A0));
		timeList.push_back(CURRENTTIME);

		
		/*
		double total_work = 0;
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++) 
		{ 
			total_work += (*(*vit)).SavedInfoForEvolve_TotalForce.DotProduct( ((*(*vit)).rnew - (*(*vit)).rold) );
		}
		*/
		

		// STEP 8: Check for whether we have hit equilibrium or too small a timestep
		if (energyList.size() > energyListSize)
		{
			(*this).equilibrated = (*this).IsEquilibrated();
			if ((*this).equilibrated)
			{
				//(*this).equilibrated = (*this).IsEquilibrated();
				*errorfile << "Equilibrium condition was met." << endl;
				cout << "Equilibrium condition was met" << endl;
				breakVar = true;
			}
		}
 
		if (breakVar)
		{
			(*this).equilibrated = true;
			if (cell_type_sum < 1)
				break;
			else
				breakVar = false;
				cell_type_sum = 0;
		}
	}

	cout << "TISSUE EVOLVED TO TIME " << CURRENTTIME << " in " << count << " steps" << endl;
	cout << "LX = " << LX << "  LY = " << LY << endl;

}
/*
* bool IsEquilibrated()
*
* checks if packing is equilibrated (this function needs to be formalized still - once we figure out an
* appropriate energy threshold). In AGPerimTissue, this means that the stepsize of the simulation is too
* small or the coarse grained derivative of the energy v time plot is sufficiently small 
*/
bool AGPerimTissue::IsEquilibrated()
{
	//bool energyratio = (abs((currentFreeEnergy - pastFreeEnergy) / (A0 * K * cells.size())) < equilibriumtolerance);
	double energySlope = getLinearRegressionSlope(timeList, energyList);
	bool equilibrated= (abs(energySlope) < equilibriumtolerance) || (*this).equilibrated;
	return equilibrated;
}


// Looks through all of the edges and find any that are to small and have them undergo T1 transitions, returns TRUE if successful.
/*!
\par Pseudocode
\code{.cpp}
startfcn:
for each edge in the tissue, do:
if length < minlength, then:
record the edge
set the flag = TRUE
endif
endfor

if flag is FALSE, then:
return FALSE
else
for each recorded edge, do:
call HelperT1_RemoveEdge()
endfor

if the T1 is successful
return TRUE
else
return FALSE
endif
endif

endfcn
\endcode
*/
//Note: to debug, take as input ofstream* datafile, then use (*this).WriteOutputDataAppend(datafile); before and after HelperT1_RemoveEdge

bool AGPerimTissue::T1(ofstream* errorfile){
	int flag = 0;
	int i;
	vector<Edge*> t1edges;								//list to store the edges to t1 in (you must put them in a seprate list and then t1 them or the number of edges will change and the itterator will not function properly
	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){
		if ((*(*eit)).GetLength() < MINLENGTH){
			t1edges.push_back(&(*(*eit)));
			*errorfile << "\n" << "t1 " << "\n" << endl;
			flag++;
		}
	}
	if (flag == 0){
		return false;
	}
	else{
		//Output state of system before and after T1.
		//*errorfile << "Free Energy BEFORE T1: " << (*this).CalculateFreeEnergy() << endl;
		//(*this).WriteOutputDataAppend(errorfile);
		for (i = 0; i<t1edges.size(); i++){
			cout << "T1 of Edge " << (*t1edges[i]).id << "   at location   (" << (*t1edges[i]).GetCenter().x << " , " << (*t1edges[i]).GetCenter().y << ")" << endl;
			int nT1;
			Vertex* v = (*this).HelperT1_RemoveEdge(t1edges[i], errorfile, &nT1);
			// (*this).WriteOutputDataAppend(errorfile);
			HelperResolveVertices_ResolveFourfoldVertex(v);
		}
		 //(*this).WriteOutputDataAppend(errorfile);
		//*errorfile << "Free Energy AFTER T1: " << (*this).CalculateFreeEnergy() << endl;
		return true;
	}
}

/*
* double CalculatePressureEnergy()
*
* Calculates the area energy of the packing using the Hamiltonian
* in this model:
* H = (sum over cells) [ 0.5 K(A_cells - A_0) ]
*/
double AGPerimTissue::CalculatePressureEnergy()
{
	double freeenergy = 0;
	
	double area_difference_sqr = 0;
	// (A_c - A_0)^2
	double carea;
	double a0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		carea = (*agit).GetArea();
		a0 = (*agit).A0;
		if ((*agit).pressureType == "LOG")
		{
			area_difference_sqr += -2 * log(1 + carea / a0);
			// where the energy would look like Ep = - 1 / 2 K A0^2 Log(1 + A/A0) ([K] = pressure)
		}
		else if ((*agit).pressureType == "EXPONENTIAL")
		{
			area_difference_sqr += 2  * exp(- carea / a0 );
			// where the energy would look like Ep = 1/2 K A0^2 exp(-A/A0) ([K] = pressure)
		}
		else if ((*agit).pressureType == "QUARTIC")
		{
			area_difference_sqr += (pow(carea / a0 - 1, 2) + C4*pow(carea/a0 - 1, 4) );
		}
		else if ((*agit).pressureType == "QUARTIC2.0")
		{
			area_difference_sqr += (pow(carea / a0 - 1, 2) * ( 1 - 6 * C4 ) + C4 * pow(carea / a0 - 1, 4));
		}
		else
		{
			area_difference_sqr += a0 * pow(carea / a0 - 1, 2);
			// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
		}
	}
	freeenergy = 0.5*K*a0*area_difference_sqr;
	return freeenergy;
}

/*
* bool CalculatePerimeterEnergy()
*
* Calculates the perimeter energy of the packing using the Hamiltonian
* in this model:
* H = (sum over cells) [  0.5 LAMBDA (P_cell - P_0)^2]
*/
double AGPerimTissue::CalculatePerimeterEnergy()
{
	double freeenergy = 0;

	double perim_difference_sqr = 0;
	// (P_c - P_0)^2


	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		perim_difference_sqr += pow(((*agit)).GetPerimeter() - (*agit).P0, 2);
		
	}
	freeenergy = 0.5*LAMBDA*perim_difference_sqr;
	return freeenergy;
}

/*
* double CalculateQuadraticPressureEnergy
* 
* Calculates the pressure energy of the packing as if it were quadratic
*/
double AGPerimTissue::CalculateQuadraticPressureEnergy()
{
	double freeenergy = 0;

	double area_difference_sqr = 0;
	// (A_c - A_0)^2
	double carea;
	double a0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		carea = (*agit).GetArea();
		a0 = (*agit).A0;
		area_difference_sqr += a0 * pow(carea / a0 - 1, 2);
		// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
	}
	freeenergy = 0.5 * K * area_difference_sqr;
	return freeenergy;
}

/*
* double CalculateLogPressureEnergy
* 
* Calculates the pressure energy of the packing as if it were logarithmic
*/
double AGPerimTissue::CalculateLogPressureEnergy()
{
	double freeenergy = 0;

	double area_difference_sqr = 0;
	// (A_c - A_0)^2
	double carea;
	double a0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		carea = (*agit).GetArea();
		a0 = (*agit).A0;
		area_difference_sqr += -2 * a0 * log(1 + carea / a0);
		// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
	}
	freeenergy = 0.5 * K * area_difference_sqr;
	return freeenergy;
}

/*
* double CalculateExpPressureEnergy
* 
* Calculates the pressure energy of the packing as if it were exponential
* 
*/
double AGPerimTissue::CalculateExpPressureEnergy()
{
	double freeenergy = 0;

	double area_difference_sqr = 0;
	// (A_c - A_0)^2
	double carea;
	double a0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		carea = (*agit).GetArea();
		a0 = (*agit).A0;
		area_difference_sqr += 2 * a0 * exp(-carea / a0);
		// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
	}
	freeenergy = 0.5 * K * area_difference_sqr;
	return freeenergy;
}

// Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration.Splits vertices so edges(e1, e2) and(e3, e4) are paired and cells c2 and c4 are neighbors retruns false if the vertex is unstable in that topology
/**
* HelperResolveVertices_CheckStability( see Meryl's thesis to see what all goes into here/ check out Tissue class )
*
* An updated stability checker for 4-fold vertices in the perimeter squared model
*
*/
double AGPerimTissue::HelperResolveVertices_CheckStability(Vertex* vert, Edge* edge1, Edge* edge2, Edge* edge3, Edge* edge4, Cell* cell1, Cell* cell2, Cell* cell3, Cell* cell4, Vect* Force){
	//First make new objects which are coppies of the objects in the linked lists we will minipulate these objecs and they will be destroyed when the function ends
	Cell c1 = (*cell1);
	Cell c2 = (*cell2);
	Cell c3 = (*cell3);
	Cell c4 = (*cell4);
	AGPerimEdge* age1 = dynamic_cast<AGPerimEdge*>(edge1);
	AGPerimEdge e1 = (*age1);
	AGPerimEdge* age2 = dynamic_cast<AGPerimEdge*>(edge2);
	AGPerimEdge e2 = (*age2);
	AGPerimEdge* age3 = dynamic_cast<AGPerimEdge*>(edge3);
	AGPerimEdge e3 = (*age3);
	AGPerimEdge* age4 = dynamic_cast<AGPerimEdge*>(edge4);
	AGPerimEdge e4 = (*age4);
	double stability;
	double a1, a2, a3, a4;

	//Make the new edge and vertices resulting from the split into 2 3fold vertices
	Edge enew;
	Vertex v12; //the vertex on edges 1 and 2
	Vertex v34; //the vertex on edges 3 and 4

	//Now update the info on the coppies of the objects so they are 2 3fold vertices

	//For the new vertices
	v12.r = (*vert).r;							//put the two new verts where the olf 4ourfold was
	v12.edges.push_back(&e1);
	v12.edges.push_back(&e2);
	v12.edges.push_back(&enew);
	v34.r = (*vert).r;							//put the two new verts where the olf 4ourfold was
	v34.edges.push_back(&e3);
	v34.edges.push_back(&e4);
	v34.edges.push_back(&enew);

	//For the new edge
	enew.yflag = 0;		//the new edge has no length and therefor cann not have a flag
	enew.xflag = 0;
	enew.vertex1 = &v12;
	enew.vertex2 = &v34;
	enew.cell1 = &c4;
	enew.cell2 = &c2;

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
	double p2 = (c2.GetPressure(a2) - c1.GetPressure(a1)); //pressure difference acorss edge2 between cell 1 and cell 2 taken CCW
	double p3 = (c3.GetPressure(a3) - c2.GetPressure(a2));
	double p4 = (c4.GetPressure(a4) - c3.GetPressure(a3));
	double p1 = (c1.GetPressure(a1) - c4.GetPressure(a4));

	//getting z X edge (both the diredction and magnitude are correctHelperResolveVertices_ResolveFourfoldVertex
	Vect zXe1 = (e1).GetVectorDifference(&v12, v12.GetOtherVertex(&e1)).Rotate(PI / 2.0);
	Vect zXe2 = (e2).GetVectorDifference(&v12, v12.GetOtherVertex(&e2)).Rotate(PI / 2.0);
	Vect zXe3 = (e3).GetVectorDifference(&v34, v34.GetOtherVertex(&e3)).Rotate(PI / 2.0);
	Vect zXe4 = (e4).GetVectorDifference(&v34, v34.GetOtherVertex(&e4)).Rotate(PI / 2.0);

	//Get the forces (tension first)
	Vect ForceT1 = (e1).GetTension_ForceOnVertex(&v12);
	Vect ForceT2 = (e2).GetTension_ForceOnVertex(&v12);
	Vect ForceT3 = (e3).GetTension_ForceOnVertex(&v34);
	Vect ForceT4 = (e4).GetTension_ForceOnVertex(&v34);
	//Generate unit vectors corresponding to edges
	//Vect l1 = ForceT1 / ForceT1.Magnitude();
	//Vect l2 = ForceT2 / ForceT2.Magnitude();
	//Vect l3 = ForceT3 / ForceT3.Magnitude();
	//Vect l4 = ForceT4 / ForceT4.Magnitude();

	//Get the forces (due to pressure)
	Vect ForceP1 = zXe1 * (p1 * 0.5);
	Vect ForceP2 = zXe2 * (p2 * 0.5);
	Vect ForceP3 = zXe3 * (p3 * 0.5);
	Vect ForceP4 = zXe4 * (p4 * 0.5);
	*Force = ((ForceT1 + ForceT2 + ForceP1 + ForceP2) - (ForceT3 + ForceT4 + ForceP3 + ForceP4))*0.5;

	

	//Calculate the tension on the new edge
	Vect ForceNewEdge_zerolength = enew.GetTension_ForceOnVertex(&v12);
	

	//Check the stability
	stability = (*Force).Magnitude() - ForceNewEdge_zerolength.Magnitude();
	// You changed the location of the factors of two that were in the previous verions of this code
	// Consider putting them back if necessary
	return stability;
}

/*
* void HelperWriteOutputData_Write2File
*
* Writes data into file specified by 'datafile'. If flag 'includeExtraInfo' is true, then additional information (such as tension, cell area, etc) is included.
*The extra information is non-essential in reconstructing the simulation state, and is just provided for convenience (e.g. for display in MATLAB).
* (should be helpful for incremental coding erro debug sessions)
*
*/

void AGPerimTissue::HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo){
	time_t rawtime;					//to deal with writing the date out
	struct tm * timeinfo;
	*datafile << "Global:" << endl;
	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;
	//Now write the global variables you used
	*datafile << std::setprecision(outputPrecision) << "CURRENTTIME\t" << CURRENTTIME << endl;
	*datafile << "SEED\t" << SEED << endl;
	*datafile << std::setprecision(outputPrecision) << "TIMEINTERVAL\t" << TIMEINTERVAL << endl;
	*datafile << std::setprecision(outputPrecision) << "LX\t" << LX << endl;
	*datafile << std::setprecision(outputPrecision) << "LY\t" << LY << endl;
	*datafile << std::setprecision(outputPrecision) << "Abox\t" << LX*LY << endl;
	*datafile << "MINAREA\t" << MINAREA << endl;
	*datafile << "MINLENGTH\t" << MINLENGTH << endl;
	*datafile << "NEWEDGELENGTH\t" << NEWEDGELENGTH << endl;
	*datafile << std::setprecision(outputPrecision) << "K\t" << K << endl;
	*datafile << std::setprecision(outputPrecision) << "GAMMA\t" << GAMMA << endl;
	*datafile << std::setprecision(outputPrecision) << "TOTALTIME\t" << TOTALTIME << endl;
	*datafile << std::setprecision(outputPrecision) << "INITIALSTEPSIZE\t" << INITIALSTEPSIZE << endl;
	*datafile << "error_tolerance_absolute_mechanical\t" << error_tolerance_absolute_mechanical << endl;
	*datafile << std::setprecision(outputPrecision) << "TISSUE_ALPHA\t" << TISSUE_ALPHA << endl;
	*datafile << std::setprecision(outputPrecision) << "ALPHAM\t" << alphaM << endl;
	*datafile << std::setprecision(outputPrecision) << "TISSUE_EXT_STRESS_X\t" << TISSUE_EXT_STRESS_X << endl;
	*datafile << std::setprecision(outputPrecision) << "TISSUE_EXT_STRESS_Y\t" << TISSUE_EXT_STRESS_Y << endl;
	*datafile << "N\t" << cells.size() << endl;
	*datafile << std::setprecision(outputPrecision) << "A0\t" << (*this).A0 << endl;
	*datafile << "APOPTOSISON\t" << APOPTOSISON << endl;
	*datafile << "T1ON\t" << T1ON << endl;
	*datafile << "MITOSISON\t" << MITOSISON << endl;
	*datafile << "WRITEMITOSIS\t" << WRITEMITOSIS << endl;
	*datafile << "VERBOSEON\t" << VERBOSEON << endl;
	*datafile << "EQUILIBRIUMMITOSIS\t" << EQUILIBRIUMMITOSIS << endl;
	*datafile << std::setprecision(outputPrecision) << "LAMBDA\t" << LAMBDA << endl;
	*datafile << "energyListSize\t" << energyListSize << endl;
	*datafile << "C4\t" << C4 << endl;
	*datafile << std::setprecision(outputPrecision) << "equilibriumtolerance\t" << equilibriumtolerance << endl;
	*datafile << std::setprecision(outputPrecision) << "P0\t" << P0 << endl;
		//******ADD NEW OUTPUT HERE******//
	*datafile << "endGlobal" << endl;

	if (!includeExtraInfo)											// if flag is false, include only essential information to capture the state of the simulation
	{
		//print the vertex info
		*datafile << "Vertices" << endl;
		*datafile << "id\tx\ty" << endl;
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){			//loop through all of the vertices in the tissue
			*datafile << std::setprecision(outputPrecision) << (*(*vit)).id << "\t" << (*(*vit)).r.x << "\t" << (*(*vit)).r.y << endl;						//write it's info to the file
		}

		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			*datafile 
				<< (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << endl;			//write its info to the file
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;
		*datafile << "id\tA0\tEdgeNumb\tEdgeList\tType" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++){				//make comma-separated list of edges
				ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
			}

			*datafile << std::setprecision(outputPrecision) << (*(*cit)).id << "\t" << (*(*cit)).A0 << "\t" << (*(*cit)).edges.size() << "\t" << ssEdgeList.str() << "\t" << (*(*cit)).type << endl;
		}
	}

	else															// if flag is true, include only essential information to capture the state of the simulation
	{
		//print the vertex info with additional information
		*datafile << "Vertices" << endl;
		*datafile << "id\tx\ty\t_EdgesNumber\tForceX\tForceY" << endl;
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
		{			
			//loop through all of the vertices in the tissue
			*datafile << std::setprecision(outputPrecision) << (*(*vit)).id << "\t" << (*(*vit)).r.x << "\t" << (*(*vit)).r.y << "\t" << (*(*vit)).edges.size() << "\t" << (*(*vit)).SavedInfoForEvolve_TotalForce.x << "\t" << (*(*vit)).SavedInfoForEvolve_TotalForce.y << endl;						//write it's info to the file
		}

		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\tlength\txflag\tyflag\ttension\t_cell1\t_cell2" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++)
		{			
			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			double L = (*(*eit)).GetLength();
			*datafile << std::setprecision(outputPrecision) << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t"<< L << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << "\t" << (*(*eit)).GetTension_Magnitude() << "\t" << (*(*(*eit)).cell1).id << "\t" << (*(*(*eit)).cell2).id << endl;
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;

		/*
		// debugging GetBaryCentricVectors
		*datafile << "Cenx\tCeny\t\tBVx\t\t\t\t\t\t\t\t\tBVy" << endl;
		for (list<Cell*>::iterator cvit = cells.begin(); cvit != cells.end(); cvit++)
		{
		stringstream bVectorListx;
		stringstream bVectorListy;
		bVectorListx.clear();
		bVectorListy.clear();
		AGCell* agit = dynamic_cast<AGCell*>((*cvit));
		Vect cent = (*agit).GetCenter();
		list<Vect> tempVList = (*agit).GetBarycentricUnitVectors();
		// include GETCenter in here too
		for (list<Vect>::iterator tvit = tempVList.begin(); tvit != tempVList.end(); tvit++){
		bVectorListx << (*tvit).x << ",";
		bVectorListy << (*tvit).y << ",";
		}
		*datafile << cent.x << "\t" << cent.y << "\t \t" << bVectorListx.str() << "\t \t \t" << bVectorListy.str() << endl;
		*/

		/* FOR NORMAL RUNS - - -
		*datafile << "id\tA0\tEdgeNumb\tType\txflag\tyflag\tArea\tPressure\tCenterx\tCentery\tEdgeList\tCWVertices\tThetaTCJ\tPotential\tTCJSTRENGTH\tTHETAAN\tANSTRENGTH" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
		{			//loop through all of the cells in the tissue
		stringstream ssEdgeList;
		ssEdgeList.clear();
		for (i = 0; i < (*(*cit)).edges.size(); i++)
		{				//make comma-separated list of edges
		ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
		}

		vector<Vertex*> CWVertices = (*(*cit)).GetClockWiseVertices();
		stringstream ssCWVertices;
		ssCWVertices.clear();
		for (i = 0; i < CWVertices.size(); i++)
		{			//make comma-separated list of vertices ordered clockwise
		ssCWVertices << (*CWVertices[i]).id << ",";
		}
		Vect centroid = (*(*cit)).GetCenter();
		AGCell * agit = dynamic_cast<AGCell*>(*cit);
		*datafile << (*(*cit)).id << "\t" << (*(*cit)).A0 << "\t" << (*(*cit)).edges.size() << "\t\t" << (*(*cit)).type << "\t" << (*(*cit)).GetXFlag() << "\t" << (*(*cit)).GetYFlag() << "\t"
		<< (*(*cit)).GetArea() << "\t" << (*(*cit)).GetPressure((*(*cit)).GetArea()) << "\t" << centroid.x << "\t" << centroid.y << "\t" << ssEdgeList.str() << "\t" << ssCWVertices.str() << "\t"<< (*agit).GetThetaTCJ() << "\t" << (*agit).DIVISIONPOTENTIAL
		<< "\t" << (*agit).TCJSTRENGTH << "\t" << (*agit).THETAAN << "\t" << (*agit).ANSTRENGTH << endl;
		}
		}
		*/
		// FOR DEBUGGING Manning model

		*datafile << std::setprecision(outputPrecision) << "id\tA0\tEdgeNumb\tType\txflag\tyflag\tArea\tPerimeter\tPressure\tCenterx\tCentery\tEdgeList\tCWVertices\tCellPressureType" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
		{			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++)
			{				//make comma-separated list of edges
				ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
			}

			vector<Vertex*> CWVertices = (*(*cit)).GetClockWiseVertices();
			stringstream ssCWVertices;
			ssCWVertices.clear();
			for (i = 0; i < CWVertices.size(); i++)
			{			//make comma-separated list of vertices ordered clockwise
				ssCWVertices << (*CWVertices[i]).id << ",";
			}
			Vect centroid = (*(*cit)).GetCenter();
			AGCell * agit = dynamic_cast<AGCell*>(*cit);
			*datafile << std::setprecision(outputPrecision) << (*(*cit)).id << "\t" << (*agit).A0 << "\t" << (*(*cit)).edges.size() << "\t" << (*(*cit)).type << "\t" << (*(*cit)).GetXFlag() << "\t" << (*(*cit)).GetYFlag() << "\t"
				<< (*agit).GetArea() << "\t" << (*agit).GetPerimeter() << "\t" << (*agit).GetPressure((*(*cit)).GetArea()) << "\t" << centroid.x << "\t" << centroid.y << "\t" << ssEdgeList.str() << "\t" << ssCWVertices.str() << "\t" << (*agit).pressureType << endl;
		}
	}

	*datafile << "End" << endl << endl;
}

/*
* bool Mitosis (overload)
*
* Do mitosis on cells. For AGPerimTissue (1/24/2023), I want to make it so that it will divide cells by type
*/
bool AGPerimTissue::Mitosis() {
	vector<Cell*> mitcells;
	AGPerimCell* ag1;
	int flag = 0;
	int i;
	bool convex_condition;

	double f = 0; //feedback strength - between 0 and 1
	// double alpha = p0 / (K*(A0 - f*A0)); //this defines feedback so that if a cell is smaller than f times the preferred area of cells, it can't divide
	double alpha = 0;
	double thisarea = 0;
	double area_ratio = 0;
	double p = 0;
	double angle;
	//int preMitosisNumCells = cells.size();// number of cells before mitosis
	//int postMitosisNumCells; // number of cells after mitosis
	//print Feedback Parameters
	/*
	if (CURRENTTIME == 1) {
		cout << "Mitosis Feedback Parameters:" << endl;
		cout << "p0: " << p0 << "		alpha: " << alpha << endl;
	}
	*/
	//Here is where we will later insert some criteria for mitosis
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		thisarea = (*(*cit)).GetArea();
		area_ratio = (thisarea - (*(*cit)).A0) / (*(*cit)).A0;
		if (area_ratio < 0)
		{
			p = 0;
		}
		else if (area_ratio > .01)
		{
			p = 1;
			flag++;
			mitcells.push_back(&(*(*cit)));
		}
		else
		{
			p = area_ratio;
			flag++;
			mitcells.push_back(&(*(*cit)));
		}
		// conditional to make certain specific cells divide
		//if ((*(*cit)).type > 0)
		//{
		//	mitcells.push_back(&(*(*cit)));
		//	(*(*cit)).type = 0;
		//	flag++;
		//	//TOTALMITOSES++;
		//}
	}
	if (flag == 0) {
		return false;
	}
	else {
		//Output the cell id and area for all cells undergoing Mitosis
		// To do mitosis on only one cell, change for statement to i=0; i<(mitcells.size())-(mitcells.size()-1); i++, for normal, i = 0; i<mitcells.size(); i++ 
		for (i = 0; i < mitcells.size(); i++) { // changed for statement to only include i=0, so mitosis will only effect one cell
			//NOTE: only perform mitosis on cells that are convex in order to avoid cutting more than two edges.
			convex_condition = (CheckConvex((*mitcells[i]).HelperOrderVerticesPos()) == true);
			if (convex_condition) {
				//cout << "Mitosis of cell \t" << (*mitcells[i]).id << endl;
				//cout << "Cell Area \t" << (*mitcells[i]).GetArea() << endl;
				//cout << "Current time is \t" << CURRENTTIME << endl;
				ag1 = dynamic_cast<AGPerimCell*>((mitcells[i]));
				cout << (*ag1).id << endl;
				angle = (*ag1).DivisionAngleDistribution((*ag1).DIVISIONPOTENTIAL);
				(*this).HelperMitosis(mitcells[i], angle);
				
			}
		}

		// postMitosisNumCells = cells.size(); // number of cells in packing after mitosis

		return true;
	}
}
/*
bool AGTissue::Mitosis()
{
	vector<Cell*> mitcells;
	AGCell* ag1;
	int flag = 0;
	int i;
	uniform_real_distribution<double> dist(0.0, 1.0);
	double probmitosis, randnumb;
	double carea;
	double angle;
	uniform_real_distribution<double> angledist(0.0, 2 * PI);

	//Here is where we will later insert some criteria for mitosis

	// THIS MAY NEED SOME CHANGING UP SO THAT WE MAKE SURE THE SYSTEM STAYS QUASI-STATIC

	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		carea = (*(*cit)).GetArea((*(*cit)).SavedInfoForEvolve_CWVerts);
		// probmitosis = 0.5*INITIALSTEPSIZE*exp(15*((carea/(1.25*A0))-1)); //probability of mitosis is based on your size
		probmitosis = 0.5*INITIALSTEPSIZE*exp(((carea / (5*A0))));
		randnumb = dist(generator);
		if ( (randnumb < probmitosis) )
		{
			//cout << "pmitosis  " << probmitosis  <<  "   rand numb   " << randnumb << endl;
			mitcells.push_back(&(*(*cit)));
			flag++;
		}
	}


	//if (CURRENTTIME<5)
	//	probmitosis = 0; //probability of mitosis
	//else
	//	probmitosis = 16 * INITIALSTEPSIZE;
	//randnumb = dist(generator);
	//if (randnumb < probmitosis) {
	//	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
	//		//Get the first cell that has not yet devided
	//		if ((*(*cit)).type == 0){
	//			mitcells.push_back(&(*(*cit)));
	//			flag++;
	//			break;
	//		}
	//	}
	//}

	if (flag == 0){
		return false;
	}
	else{
		//Output the cell id and area for all cells undergoing Mitosis
		for (i = 0; i<mitcells.size(); i++){
			//NOTE: only perform mitosis on cells that are convex in order to avoid cutting more than two edges.
			if (CheckConvex((*mitcells[i]).HelperOrderVerticesPos()) == true){
				ag1 = dynamic_cast<AGCell*>((mitcells[i]));
				angle = (*ag1).DivisionAngleDistribution();
				(*this).HelperMitosis(mitcells[i], angle);
			}
		}
		return true;
	}

}
*/


/*
* bool CreateCell (overload)
*
* Creates a new cell and transfers appropriate fields to the new daughter cells.
*/
Cell* AGPerimTissue::CreateCell(int typ, vector<Edge*> edgs, string potential, double tcjstrength, double thetaan, double anstrength)
{
	int largestID, newid;
	list<Cell*>::iterator cit;
	Cell* newcell{};
	//Get the id
	largestID = (*max_element(cellIDs.begin(), cellIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;								//Now make a new ID that is unique

	if (newid > 1000000000) { cout << "ERROR to many cells" << endl; }
	//create the agcell and add the agcell to the lists
	AGPerimCell* cptr = new AGPerimCell;
	(*cptr).id = newid;
	(*cptr).type = typ;
	(*cptr).edges = edgs;
	(*cptr).TCJSTRENGTH = tcjstrength;
	(*cptr).DIVISIONPOTENTIAL = potential;
	(*cptr).ANSTRENGTH = anstrength;
	(*cptr).THETAAN = thetaan;
	(*cptr).P0 = P0;
	(cells).push_back(cptr);
	cellIDs.push_back(newid);

	for (cit = cells.begin(); cit != cells.end(); cit++) {
		if ((*(*cit)).id == newid) {
			// return &(*(*cit)); - this was the way the old code worked. Linux didn't like this though, so I added newcell.
			newcell = &(*(*cit));
			break;
		}
	}
	return newcell;
}

/*
* void HelperMitosis( Cell* parent cell, double divisionangle)
*
* Called by Mitosis to do the actual pointer and object manipulation required for mitosis
*/
void AGPerimTissue::HelperMitosis(Cell* parentcell, double divisionangle) {
	int i;										//for loops
	Vect center, p1, p2, q1, q2;				//Used to find newcenteredge
	Vect rnewvtx;								//the positions of the new vertices
	vector<Edge*> empty;
	vector<Edge*> deletededges;					//vector of edge pointers to the two edges to remove
	vector<Edge*> newedges;						//vector to hold the four new edges to replace the two that are being removed
	vector<Vertex*> newvertices;				//vector of pointers to the two new vertices
	vector<Vertex*> cwverts;					//parent cell verts in clockwise order
	Edge* newedgecenter = 0;					//ptr to the edge through the center of the cell
	Vertex* v1ptr, * v2ptr, * v3ptr, * v4ptr;
	double angle;								//angle at which to cut through the cell
	vector<vector<double>> vertposforflag;		//used to define vertex positions if any flags are set
	bool newxflag = false;						//flags for which side of boundary newvtxr is on
	bool newyflag = false;						//flags for which side of boundary newvtxr is on					
	vector<Vertex*> vertsordered;				//list of pointers to vertices used for the case of flags
	double xpossum = 0;							//sum of x positions of vertices in cell
	double ypossum = 0;							//sum of y positions of vertices in cell
	double centerposx = 0;						//used for center calculation done within such that the center is consistent with the vertposforflag
	double centerposy = 0;						//used for center calculation done within such that the center is consistent with the vertposforflag
	vector<int> edgeflagstestind, edgeflagstestx, edgeflagstesty;	//this is used to test for flags on four newly created edges (not one created in center); stores i and two flags
	int vertnumb;								//lists number of vertices
	vector<int> newflagstestx, newflagstesty;	//this is used to flag the newedgecenter
	Cell* c1ptr = NULL;
	Cell* c2ptr = NULL;
	Cell* daughter1ptr = NULL;
	Cell* daughter2ptr = NULL;
	AGPerimCell* momptr = dynamic_cast<AGPerimCell*>((parentcell));

	vertsordered = (*parentcell).HelperOrderVertices();  //stores vertex information to be called after vertposforflag is checked
	vertposforflag = (*parentcell).HelperOrderVerticesPos(); //allows positions of vertices outside of domain
	vertnumb = (*parentcell).GetVertexNumber(); //stores number of vertices in cell

	//This is the geometrical center of the polygonal cell (This center can be outside of the boundaries depending on HelperVerticesPos() values).
	center = GetGeoCenter(vertposforflag[0], vertposforflag[1]);


	//get the angle at which to put the new edge
	angle = divisionangle;
	q1 = Vect(center, angle, LX);      //make two points that define a line at the given angle through the center of the cell
	q2 = Vect(center, angle, -LX);

	//create two new vertices where this angle intersects the edges of the cell with if statement that tests for flag and directs accordingly
	if ((*parentcell).GetFlag() == 0) {
		for (i = 0; i < (*parentcell).edges.size(); i++) {								//look at each edge
			p1 = (*(*(*parentcell).edges[i]).vertex1).r;								//get it's end points and put them into p1 and p2
			p2 = (*(*(*parentcell).edges[i]).vertex2).r;
			if (test_for_Intersect(q1, q2, p1, p2, &rnewvtx)) {							//if it crosses the selected angle
				deletededges.push_back((*parentcell).edges[i]);								//record the edge
				newvertices.push_back((*this).CreateVertex(rnewvtx, empty));				//make a new vertex at the interscetion and record its ptr
			}
		}

		//create a new edge between the two new vertices (must have no xflag and no yflag)
		newedgecenter = (*this).CreateEdge(newedgecenter, 0, 0, newvertices[0], newvertices[1], parentcell, parentcell);
	}

	//This is for the case where a cell crosses any boundary (WE LOOP THROUGH AND THEN PAIR THE END WITH THE BEGINNING)
	if ((*parentcell).GetFlag() == 1) {
		for (i = 0; i < vertnumb - 1; i++) { //look at each vertex
			p1.x = (vertposforflag[0])[i]; //get it's end points and put them into p1 and p2 (THESE POSITIONS CAN BE OUT OF BOUNDS)
			p1.y = (vertposforflag[1])[i];
			p2.x = (vertposforflag[0])[i + 1];
			p2.y = (vertposforflag[1])[i + 1];
			if (test_for_Intersect(q1, q2, p1, p2, &rnewvtx)) {				//if it crosses the selected angle
				deletededges.push_back((*vertsordered[i]).GetEdgeBetween(vertsordered[i + 1]));								//record the edge			
				if (rnewvtx.x > LX || rnewvtx.x < 0) { newxflag = 1; } //These six lines used later to create flags for newedgecenter (stores information about which side of boundary for each newvertex)
				else { newxflag = 0; }
				if (rnewvtx.y > LY || rnewvtx.y < 0) { newyflag = 1; }
				else { newyflag = 0; }
				newflagstestx.push_back(newxflag);
				newflagstesty.push_back(newyflag);
				edgeflagstestind.push_back(i); //This whole sections helps us later when trying to set flags for other four newly created edges (TESTS WHICH SIDE OF THE BOUNDARY ARE THE FOUR VERTICES OF DELETED EDGES)
				if (p1.x > LX || p1.x < 0) { edgeflagstestx.push_back(1); }
				else { edgeflagstestx.push_back(0); }
				if (p1.y > LY || p1.y < 0) { edgeflagstesty.push_back(1); }
				else { edgeflagstesty.push_back(0); }
				edgeflagstestind.push_back(i + 1);
				if (p2.x > LX || p2.x < 0) { edgeflagstestx.push_back(1); }
				else { edgeflagstestx.push_back(0); }
				if (p2.y > LY || p2.y < 0) { edgeflagstesty.push_back(1); }
				else { edgeflagstesty.push_back(0); } //Conclusion of section used for setting flags for the other four newly created edges
				if (rnewvtx.x > LX) { rnewvtx.x = rnewvtx.x - LX; } //This part puts the newvtx within the boundaries appropriately
				if (rnewvtx.x < 0) { rnewvtx.x = rnewvtx.x + LX; }
				if (rnewvtx.y > LY) { rnewvtx.y = rnewvtx.y - LY; }
				if (rnewvtx.y < 0) { rnewvtx.y = rnewvtx.y + LY; }
				newvertices.push_back((*this).CreateVertex(rnewvtx, empty));
			}
		}
		p1.x = (vertposforflag[0])[vertnumb - 1]; //get it's end points and put them into p1 and p2 FOR THE CASE WHERE THE BOTTOM ELEMENT IS PAIRED WITH THE TOP ELEMENT
		p1.y = (vertposforflag[1])[vertnumb - 1];
		p2.x = (vertposforflag[0])[0];
		p2.y = (vertposforflag[1])[0];
		if (test_for_Intersect(q1, q2, p1, p2, &rnewvtx)) {				//if it crosses the selected angle
			deletededges.push_back((*(vertsordered[vertnumb - 1])).GetEdgeBetween(vertsordered[0]));								//record the edge
			if (rnewvtx.x > LX || rnewvtx.x < 0) { newxflag = 1; } //These six lines used to store information about flags for newedgecenter (stores information about which side of boundary for each newvertex)
			else { newxflag = 0; }
			if (rnewvtx.y > LY || rnewvtx.y < 0) { newyflag = 1; }
			else { newyflag = 0; }
			newflagstestx.push_back(newxflag);
			newflagstesty.push_back(newyflag);
			edgeflagstestind.push_back(vertnumb - 1);
			if (p1.x > LX || p1.x < 0) { edgeflagstestx.push_back(1); } //This whole sections helps us later when trying to set flags for four other newly created edges (TESTS WHICH SIDE OF THE BOUNDARY ARE THE FOUR VERTICES OF DELETED EDGES)
			else { edgeflagstestx.push_back(0); }
			if (p1.y > LY || p1.y < 0) { edgeflagstesty.push_back(1); }
			else { edgeflagstesty.push_back(0); }
			edgeflagstestind.push_back(0);
			if (p2.x > LX || p2.x < 0) { edgeflagstestx.push_back(1); }
			else { edgeflagstestx.push_back(0); }
			if (p2.y > LY || p2.y < 0) { edgeflagstesty.push_back(1); }
			else { edgeflagstesty.push_back(0); } //Conclusion of section used for setting flags for four other newly created edges
			if (rnewvtx.x > LX) { rnewvtx.x = rnewvtx.x - LX; }//This part puts the newvtx within the boundaries appropriately
			if (rnewvtx.x < 0) { rnewvtx.x = rnewvtx.x + LX; }
			if (rnewvtx.y > LY) { rnewvtx.y = rnewvtx.y - LY; }
			if (rnewvtx.y < 0) { rnewvtx.y = rnewvtx.y + LY; }
			newvertices.push_back((*this).CreateVertex(rnewvtx, empty));
		}

		//create a new edge between the two new vertices with appropriate flags
		//test whether new vertices are on opposite side of a boundary
		if (newflagstestx[0] != newflagstestx[1] && newflagstesty[0] != newflagstesty[1])
		{
			newedgecenter = (*this).CreateEdge(newedgecenter, 1, 1, newvertices[0], newvertices[1], parentcell, parentcell);
		}
		else if (newflagstestx[0] == newflagstestx[1] && newflagstesty[0] != newflagstesty[1])
		{
			newedgecenter = (*this).CreateEdge(newedgecenter, 0, 1, newvertices[0], newvertices[1], parentcell, parentcell);
		}
		else if (newflagstestx[0] != newflagstestx[1] && newflagstesty[0] == newflagstesty[1])
		{
			newedgecenter = (*this).CreateEdge(newedgecenter, 1, 0, newvertices[0], newvertices[1], parentcell, parentcell);
		}
		else if (newflagstestx[0] == newflagstestx[1] && newflagstesty[0] == newflagstesty[1])
		{
			newedgecenter = (*this).CreateEdge(newedgecenter, 0, 0, newvertices[0], newvertices[1], parentcell, parentcell);
		}
	}

	//Assign vertex ids to newedgecenter for output purposes
	(*newedgecenter).IOVertexID[0] = (*(newvertices[0])).id;
	(*newedgecenter).IOVertexID[1] = (*(newvertices[1])).id;

	//This is just used to warn you if for some reason more than two edges of the cell are intersected
	if (deletededges.size() != 2) {
		cout << "ERROR: The Mitosis function had number of intersected edges different than two (Problem with Vertices and GetCenter or with overlap)" << endl;
	}
	(*newvertices[0]).edges.push_back(newedgecenter);
	(*newvertices[1]).edges.push_back(newedgecenter);

	//create the four new edges to replace the ones to be removed and add them to the vertices
	c1ptr = (*deletededges[0]).GetOtherCell(parentcell);
	c2ptr = (*deletededges[1]).GetOtherCell(parentcell);
	v1ptr = (*deletededges[0]).vertex1;
	v2ptr = (*deletededges[0]).vertex2;
	v3ptr = (*deletededges[1]).vertex1;
	v4ptr = (*deletededges[1]).vertex2;

	//Making sure to give all new edges the approriate IOVertexID's and that the vertices know their edges
	newedges.push_back((*this).CreateEdge(newedgecenter, 0, 0, v1ptr, newvertices[0], parentcell, c1ptr));
	(*v1ptr).edges.push_back(newedges[0]);
	(*newvertices[0]).edges.push_back(newedges[0]);
	(*c1ptr).edges.push_back(newedges[0]);
	(*parentcell).edges.push_back(newedges[0]);
	(*(newedges[0])).IOVertexID[0] = (*v1ptr).id;
	(*(newedges[0])).IOVertexID[1] = (*newvertices[0]).id;

	newedges.push_back((*this).CreateEdge(newedgecenter, 0, 0, v2ptr, newvertices[0], parentcell, c1ptr));
	(*v2ptr).edges.push_back(newedges[1]);
	(*newvertices[0]).edges.push_back(newedges[1]);
	(*c1ptr).edges.push_back(newedges[1]);
	(*parentcell).edges.push_back(newedges[1]);
	(*(newedges[1])).IOVertexID[0] = (*v2ptr).id;
	(*(newedges[1])).IOVertexID[1] = (*newvertices[0]).id;

	newedges.push_back((*this).CreateEdge(newedgecenter, 0, 0, v3ptr, newvertices[1], parentcell, c2ptr));
	(*v3ptr).edges.push_back(newedges[2]);
	(*newvertices[1]).edges.push_back(newedges[2]);
	(*c2ptr).edges.push_back(newedges[2]);
	(*parentcell).edges.push_back(newedges[2]);
	(*(newedges[2])).IOVertexID[0] = (*v3ptr).id;
	(*(newedges[2])).IOVertexID[1] = (*newvertices[1]).id;

	newedges.push_back((*this).CreateEdge(newedgecenter, 0, 0, v4ptr, newvertices[1], parentcell, c2ptr));
	(*v4ptr).edges.push_back(newedges[3]);
	(*newvertices[1]).edges.push_back(newedges[3]);
	(*c2ptr).edges.push_back(newedges[3]);
	(*parentcell).edges.push_back(newedges[3]);
	(*(newedges[3])).IOVertexID[0] = (*v4ptr).id;
	(*(newedges[3])).IOVertexID[1] = (*newvertices[1]).id;


	//This sequence of if statements checks if the edge to be deleted crosses a boundary and sees (based on if the newedgecenter crosses the boundary)
	//which of the two newly created edges on deleted edges needs a flag. This relies heavily on a testflag I set above to see if HelperVerticesPos() took a vertex out of the domain.
	if (edgeflagstestind.size() == 4) { //This just makes sure that this only runs when the indices of the four vertices of deleted edges are well known
		if ((*newedgecenter).xflag == 1) { //if the new edge crosses the x bound
			if (edgeflagstestx[0] != edgeflagstestx[1]) { //if the deleted edge vertices are on opposite sides of the boundary
				for (i = 0; i < 2; i++) { //for first deleted edge
					if (edgeflagstestx[i] != newflagstestx[0]) { //which one is opposite the new vertex
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[0])).xflag = 1;
					}
				}
			}
			if (edgeflagstestx[2] != edgeflagstestx[3]) { //if the deleted edge vertices are on opposite sides of the boundary
				for (i = 2; i < 4; i++) { //for second deleted edge
					if (edgeflagstestx[i] != newflagstestx[1]) { //which one is opposite the new vertex
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[1])).xflag = 1;
					}
				}
			}
		}
		if ((*newedgecenter).xflag == 0 && (*parentcell).GetXFlag() == 1) { //if newedgecenter doesn't cross x bound, but cell does
			if (edgeflagstestx[0] != edgeflagstestx[1]) { //if the deleted edge vertices are on opposite sides of the boundary
				for (i = 0; i < 2; i++) { //for first deleted edge
					if (edgeflagstestx[i] != newflagstestx[0]) { //which one is opposite the new vertex
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[0])).xflag = 1;
					}
				}
			}
			if (edgeflagstestx[2] != edgeflagstestx[3]) {
				for (i = 2; i < 4; i++) { //for second edge
					if (edgeflagstestx[i] != newflagstestx[1]) {
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[1])).xflag = 1;
					}
				}
			}
		}
		if ((*newedgecenter).yflag == 1) { //if newedgecenter crosses the y boundary
			if (edgeflagstesty[0] != edgeflagstesty[1]) { //if vertices are on opposite side of boundary
				for (i = 0; i < 2; i++) { //for first deleted edge
					if (edgeflagstesty[i] != newflagstesty[0]) { //which one is opposite the new vertex
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[0])).yflag = 1;
					}
				}
			}
			if (edgeflagstesty[2] != edgeflagstesty[3]) {
				for (i = 2; i < 4; i++) { //for second deleted edge
					if (edgeflagstesty[i] != newflagstesty[1]) { //which one is opposite the vertex
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[1])).yflag = 1;
					}
				}
			}
		}
		if ((*newedgecenter).yflag == 0 && (*parentcell).GetYFlag() == 1) { //if newedgecenter doesn't cross y bound, but cell does
			if (edgeflagstesty[0] != edgeflagstesty[1]) { //if vertices on opposite side of boundary
				for (i = 0; i < 2; i++) { //for first deleted edge
					if (edgeflagstesty[i] != newflagstesty[0]) { //which one is opposite new vertex
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[0])).yflag = 1;
					}
				}
			}
			if (edgeflagstesty[2] != edgeflagstesty[3]) { //if vertices on opposite side of the boundary
				for (i = 2; i < 4; i++) { //for second deleted edge
					if (edgeflagstesty[i] != newflagstesty[1]) { //which one is opposite the new vertex
						(*(*vertsordered[edgeflagstestind[i]]).GetEdgeBetween(newvertices[1])).yflag = 1;
					}
				}
			}
		}
	}
	//end of sequence of if statements. All flags of new edges should be correct.


	//Create the two new daughter cells 
	daughter1ptr = (*this).CreateCell(0, empty, (*momptr).DIVISIONPOTENTIAL, (*momptr).TCJSTRENGTH, (*momptr).THETAAN, (*momptr).ANSTRENGTH);
	daughter2ptr = (*this).CreateCell(0, empty, (*momptr).DIVISIONPOTENTIAL, (*momptr).TCJSTRENGTH, (*momptr).THETAAN, (*momptr).ANSTRENGTH);
	//PRIMARY REASON I HAD TO OVERRIDE ^^^^ TO UPDPATE CREATECELL

	//Delete the two old edges 
	(*parentcell).RemoveEdge(deletededges[0]);
	(*parentcell).RemoveEdge(deletededges[1]);
	(*c1ptr).RemoveEdge(deletededges[0]);
	(*c2ptr).RemoveEdge(deletededges[1]);
	(*v1ptr).RemoveEdge(deletededges[0]);
	(*v2ptr).RemoveEdge(deletededges[0]);
	(*v3ptr).RemoveEdge(deletededges[1]);
	(*v4ptr).RemoveEdge(deletededges[1]);
	(*this).DestroyEdge(deletededges[0]);
	(*this).DestroyEdge(deletededges[1]);

	//fill the new daughter cells edge lists and update the edges to include the new daughter cell
	(*daughter1ptr).edges.push_back(newedgecenter);
	(*daughter2ptr).edges.push_back(newedgecenter);
	(*newedgecenter).UpdateCell(parentcell, daughter1ptr);
	(*newedgecenter).UpdateCell(parentcell, daughter2ptr);

	Edge e, e0, * eptr;
	Vertex v, v0;
	e0 = *newedgecenter;
	v0 = *(newvertices[0]);
	do {							//for the first daughter				
		for (i = 0; i < (*parentcell).edges.size(); i++) {																					// look through all of the edges in the cell
			if (((*(*(*parentcell).edges[i]).vertex1) == v0 || (*(*(*parentcell).edges[i]).vertex2) == v0) && (*(*parentcell).edges[i]) != e0) {		//find the edge that has v0 as a vertex and is not the edge e0 (the one you just looked at)
				e = (*(*parentcell).edges[i]);
				eptr = (*parentcell).edges[i];
			}
		}
		(*daughter1ptr).edges.push_back(eptr);					//put the edge in the first daughter cell
		(*eptr).UpdateCell(parentcell, daughter1ptr); //and update the edge's cells
		(*daughter1ptr).IOEdgeID.push_back((*eptr).id); //build IOEdgeID for output
		v = *(v0.GetOtherVertex(eptr));							//use that edge to get the next vertex
		e0 = e;													//update e0,e,v0,v
		v0 = v;
	} while (v0 != (*newvertices[1]));
	do {
		for (i = 0; i < (*parentcell).edges.size(); i++) {																					// look through all of the edges in the cell
			if (((*(*(*parentcell).edges[i]).vertex1) == v0 || (*(*(*parentcell).edges[i]).vertex2) == v0) && (*(*parentcell).edges[i]) != e0) {		//find the edge that has v0 as a vertex and is not the edge e0 (the one you just looked at)
				e = (*(*parentcell).edges[i]);
				eptr = (*parentcell).edges[i];
			}
		}						//for the second daughter
		(*daughter2ptr).edges.push_back(eptr);						//put the edge in the second daughter cell
		(*eptr).UpdateCell(parentcell, daughter2ptr); //and update the edge's cells
		(*daughter2ptr).IOEdgeID.push_back((*eptr).id); //build IOEdgeID for output
		v = *(v0.GetOtherVertex(eptr));								//use that edge to get the next vertex
		e0 = e;														//update e0,e,v0,v
		v0 = v;
	} while (v0 != (*newvertices[0]));

	//Set A0 for cells (otherwise integrator won't know what to do)
	(*daughter1ptr).A0 = (*parentcell).A0;
	(*daughter2ptr).A0 = (*parentcell).A0;

	//Add newedgecenter to IOEDGEID for both cells
	(*daughter1ptr).IOEdgeID.push_back((*newedgecenter).id);
	(*daughter2ptr).IOEdgeID.push_back((*newedgecenter).id);

	//Define cell types (will become important when discussing compartment boundaries)
	(*daughter1ptr).type = 0;//(*parentcell).type + 1;
	(*daughter2ptr).type = 0; //(*parentcell).type + 1;

	// ^^^^^^^ THIS MIGHT NEED SOME CHANGING!

	//delete the parent cell
	(*this).DestroyCell(parentcell);
}


/*
* MatrixXf calculateHessian()
*
* calculates the Hessian of the AGPerimTissue
*
* NEEDS UPDATING IF WE ALLOW FOR MITOSIS OR APOPTOSIS
*
* CHECK HOW TO WORK OUT SYMMETRY TO MAKE THIS FASTER
*
* MOST LIKELY NEEDS TO BE DEBUGGED BY INDIVIDUALLY CALCULATING THE DERIVATIVE TERMS IN THEIR OWN FUNCTION
*/
MatrixXf AGPerimTissue::calculateHessian(ofstream* errorfile)
{
	int numofvertices = vertices.size();
	int vnum;
	MatrixXf Hessian(2 * numofvertices, 2 * numofvertices);
	Hessian.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf dAdridAdrj(2 * numofvertices, 2 * numofvertices);
	dAdridAdrj.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf dAdridAdrj_percell(2 * numofvertices, 2 * numofvertices);
	dAdridAdrj_percell.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf dPdridPdrj(2 * numofvertices, 2 * numofvertices);
	dPdridPdrj.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf dPdridPdrj_percell(2 * numofvertices, 2 * numofvertices);
	dPdridPdrj_percell.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf d2Adridrj(2 * numofvertices, 2 * numofvertices);
	d2Adridrj.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf d2Adridrj_percell(2 * numofvertices, 2 * numofvertices);
	d2Adridrj_percell.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf d2Pdridrj_percell(2 * numofvertices, 2 * numofvertices);
	d2Pdridrj_percell.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf d2Pdridrj(2 * numofvertices, 2 * numofvertices);
	d2Pdridrj.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf incidenceMatrix(2 * numofvertices, 2 * numofvertices);
	incidenceMatrix.setZero(2 * numofvertices, 2 * numofvertices);

	MatrixXf cellincidence(2 * numofvertices, 2 * numofvertices);
	cellincidence.setZero(2 * numofvertices, 2 * numofvertices);

	vector<Vertex*> cwverts;
	//vector<Edge*> cwedgesp1, cwedgesm1;
	Edge* cwedge;
	Edge* ccwedge;
	vector<Vect> dAdr, dPdr, ep1vec, em1vec; 
	vector<int> vertID, vertIDp1, vertIDm1;
	double area = 0;
	double cell_pressure = 0;
	double perimeter = 0;
	double tension_magnitude = 0;
	double a0, p0;
	Vect CWvec_pos, vec_pos, CCWvec_pos, CWedgevec, CCWedgevec, CWedgevec_UV, CCWedgevec_UV, dAdrvec, dPdrvec;
	double CWmag, CCWmag;

	// I need to reindex the vertex indices. apparently they're not just 1-N in ascending order for some reason. 
	unordered_map<int, int> vertex_index_map;
	int index = 0;
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
	{
		// The fastest way to do this is with an unordered map (the C++ equivalent to Python dictionaries)
		vertex_index_map[(*(*vit)).id] = index;
		index++;

	}
	// and now we calculate the different components, cell-wise

	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		*errorfile << "Cell id: " << (*(*cit)).id << endl;
		area = ((*agit)).GetArea();
		perimeter = ((*agit)).GetPerimeter();
		p0 = (*agit).P0;
		a0 = ((*agit)).A0;
		cwverts = (*(*cit)).SavedInfoForEvolve_CWVerts;	//Get the cw vertices in the cell from the preassigned list 
		vnum = cwverts.size();
		//cwverts.push_back(cwverts[0]); // we'll make the vertex list extra long so that it loops around.
		// We effectively have 4 terms we need in the Hessian.

		// the pressure
		cell_pressure = (*(*cit)).GetPressure(area) * a0; // I want pressure to be PRESSURE = K (A - A0)
		// tension magnitude should be simple enough
		tension_magnitude =  LAMBDA * (perimeter - p0);

		// double K = 1; // this is for debugging
		// double LAMBDA = 1; // this is for debugging

		// I want to clear the old vectors each iteration of the loop. 
		if (dAdr.size() > 0) { dAdr.clear(); }
		if (dPdr.size() > 0) { dPdr.clear(); }
		if (em1vec.size() > 0) { em1vec.clear(); }
		if (ep1vec.size() > 0) { ep1vec.clear(); }
		if (vertIDp1.size() > 0) { vertIDp1.clear(); }
		if (vertID.size() > 0) { vertID.clear(); }
		if (vertIDm1.size() > 0) { vertIDm1.clear(); }

		for (int vindex = 0; vindex < vnum ; vindex++)
		{
			// vertID, vertIDp1, vertIDm1

			// get the ID and position list for each vertex in current cell
			
			if (vindex == vnum - 1)
			{
				// CWvec_pos = (*cwverts[0]).r;
				cwedge = (*cwverts[vindex]).GetEdgeBetween(cwverts[0]);
				CWedgevec = cwedge->GetVectorDifference(cwverts[vindex],cwverts[0]);
				vertIDp1.push_back(vertex_index_map [ (*cwverts[0]).id ]);

			}
			else
			{
				// CWvec_pos = (*cwverts[vindex + 1]).r;
				cwedge = (*cwverts[vindex]).GetEdgeBetween(cwverts[vindex + 1]);
				CWedgevec = cwedge->GetVectorDifference(cwverts[vindex], cwverts[vindex + 1]);
				vertIDp1.push_back(vertex_index_map [ (*cwverts[vindex + 1]).id ] );
			}

			// vec_pos = (*cwverts[vindex]).r;
			vertID.push_back(vertex_index_map [ (*cwverts[vindex]).id ]);

			if (vindex == 0)
			{
				// CCWvec_pos = (*cwverts[vnum - 1]).r;
				ccwedge = (*cwverts[vindex]).GetEdgeBetween(cwverts[vnum - 1]);
				CCWedgevec = ccwedge->GetVectorDifference(cwverts[vnum - 1], cwverts[vindex]);
				vertIDm1.push_back(vertex_index_map [ (*cwverts[vnum - 1]).id ] );
			}
			else
			{
				// CCWvec_pos = (*cwverts[vindex - 1]).r;
				ccwedge = (*cwverts[vindex]).GetEdgeBetween(cwverts[vindex - 1]);
				CCWedgevec = ccwedge->GetVectorDifference(cwverts[vindex - 1], cwverts[vindex]);
				vertIDm1.push_back(vertex_index_map [ (*cwverts[vindex - 1]).id ]);
			}
			 dAdrvec = (CWedgevec + CCWedgevec).Rotate(PI / 2) * 0.5;


			//CWvec_pos - vec_pos;
			//vec_pos - CCWvec_pos;

			CWmag = cwedge->GetLength();
			CCWmag = ccwedge->GetLength();

			CWedgevec_UV = CWedgevec / CWmag;
			CCWedgevec_UV = CCWedgevec / CCWmag;
			dPdrvec = (CWedgevec_UV * -1 + CCWedgevec_UV);

			ep1vec.push_back(CWedgevec);
			em1vec.push_back(CCWedgevec);
			dAdr.push_back(dAdrvec);
			dPdr.push_back(dPdrvec);
		}

		int i_id, j_id, jp1_id, jm1_id;
		//vector<vector<double>> dAdidj{ { 0, 0 }, { 0, 0 } };
		//vector<vector<double>> dPdidj{ { 0, 0 }, { 0, 0 } };
		//vector<vector<double>> dHdidj_c{ { 0, 0 }, { 0, 0 } };
		// In the second loop, we'll do a nested calculation and get the 2nd order derivatives as well


		// DEBUGGING LINE TO SEE IF THE TRAVERSAL IS MAKING SENSE
		cellincidence.setZero(2 * numofvertices, 2 * numofvertices);
		d2Pdridrj_percell.setZero(2 * numofvertices, 2 * numofvertices);
		dPdridPdrj_percell.setZero(2 * numofvertices, 2 * numofvertices);
		d2Adridrj_percell.setZero(2 * numofvertices, 2 * numofvertices);
		dAdridAdrj_percell.setZero(2 * numofvertices, 2 * numofvertices);

		for (int j = 0; j < vnum; j++)
		{
			j_id = vertID[j];
			jp1_id = vertIDp1[j];
			jm1_id = vertIDm1[j];

			for (int i = 0; i < vnum; i++)
			{

				//*errorfile << "incidence = \n" << incidenceMatrix << "\n" << "\n" << endl;
				i_id = vertID[i];

				cellincidence(2 * j_id, 2 * i_id) += 1;
				cellincidence(2 * j_id, 2 * i_id + 1) += 1;
				cellincidence(2 * j_id + 1, 2 * i_id + 1) += 1;
				cellincidence(2 * j_id + 1, 2 * i_id) += 1;

				// the dAdri * dAdrj portion
				dAdridAdrj_percell(2 * j_id, 2 * i_id) += K * dAdr[j].x * dAdr[i].x;
				dAdridAdrj_percell(2 * j_id, 2 * i_id + 1) += K * dAdr[j].x * dAdr[i].y;
				dAdridAdrj_percell(2 * j_id + 1, 2 * i_id + 1) += K * dAdr[j].y * dAdr[i].y;
				dAdridAdrj_percell(2 * j_id + 1, 2 * i_id) += K * dAdr[j].y * dAdr[i].x;

				// the dPdri * dPdrj portion
				dPdridPdrj_percell(2 * j_id, 2 * i_id) += LAMBDA * dPdr[j].x * dPdr[i].x;
				dPdridPdrj_percell(2 * j_id, 2 * i_id + 1) += LAMBDA * dPdr[j].x * dPdr[i].y;
				dPdridPdrj_percell(2 * j_id + 1, 2 * i_id + 1) += LAMBDA * dPdr[j].y * dPdr[i].y;
				dPdridPdrj_percell(2 * j_id + 1, 2 * i_id) += LAMBDA * dPdr[j].y * dPdr[i].x;


			}


			// the second order derivatives section
			//*errorfile << "indices: " << jm1_id << " " << j_id << " " << jp1_id << endl;
			// first the area portion - these only act on cases where j is not equal to i and off-diagonal components
			d2Adridrj_percell(2 * jm1_id, 2 * j_id + 1) -= 0.5 * cell_pressure;
			d2Adridrj_percell(2 * jm1_id + 1, 2 * j_id) += 0.5 * cell_pressure;
			d2Adridrj_percell(2 * jp1_id, 2 * j_id + 1) += 0.5 * cell_pressure;
			d2Adridrj_percell(2 * jp1_id + 1, 2 * j_id ) -= 0.5 * cell_pressure;

			// perimeter portion - where i = j
			d2Pdridrj_percell(2 * j_id, 2 * j_id) += tension_magnitude * (1 / em1vec[j].Magnitude() + 1 / ep1vec[j].Magnitude() - pow(em1vec[j].x, 2) / pow(em1vec[j].Magnitude(), 3) - pow(ep1vec[j].x, 2) / pow(ep1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * j_id, 2 * j_id + 1) -= tension_magnitude * ((em1vec[j].x * em1vec[j].y) / pow(em1vec[j].Magnitude(), 3) + (ep1vec[j].x * ep1vec[j].y) / pow(ep1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * j_id + 1, 2 * j_id + 1) += tension_magnitude * (1 / em1vec[j].Magnitude() + 1 / ep1vec[j].Magnitude() - pow(em1vec[j].y, 2) / pow(em1vec[j].Magnitude(), 3) - pow(ep1vec[j].y, 2) / pow(ep1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * j_id + 1, 2 * j_id) -= tension_magnitude * ((em1vec[j].x * em1vec[j].y) / pow(em1vec[j].Magnitude(), 3) + (ep1vec[j].x * ep1vec[j].y) / pow(ep1vec[j].Magnitude(), 3));

			// perimeter portion - where j = i-1
			d2Pdridrj_percell(2 * jm1_id, 2 * j_id) -= tension_magnitude * (1 / em1vec[j].Magnitude() - pow(em1vec[j].x, 2) / pow(em1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * jm1_id, 2 * j_id + 1) +=  tension_magnitude * ((em1vec[j].x * em1vec[j].y) / pow(em1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * jm1_id + 1, 2 * j_id + 1) -= tension_magnitude * (1 / em1vec[j].Magnitude() - pow(em1vec[j].y, 2) / pow(em1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * jm1_id + 1, 2 * j_id) += tension_magnitude * ((em1vec[j].x * em1vec[j].y) / pow(em1vec[j].Magnitude(), 3));

			// perimeter portion - where j = i+1
			d2Pdridrj_percell(2 * jp1_id, 2 * j_id) -= tension_magnitude * (1 / ep1vec[j].Magnitude() - pow(ep1vec[j].x, 2) / pow(ep1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * jp1_id, 2 * j_id + 1) += tension_magnitude * ((ep1vec[j].x * ep1vec[j].y) / pow(ep1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * jp1_id + 1, 2 * j_id + 1) -= tension_magnitude * (1 / ep1vec[j].Magnitude() - pow(ep1vec[j].y, 2) / pow(ep1vec[j].Magnitude(), 3));
			d2Pdridrj_percell(2 * jp1_id + 1, 2 * j_id) += tension_magnitude * ((ep1vec[j].x * ep1vec[j].y) / pow(ep1vec[j].Magnitude(), 3));

		}


		//*errorfile << "H = \n" << Hessian << "\n" << "\n" << endl;
		//*errorfile << "incidence after end of cell = \n" << cellincidence << "\n" << "\n" << endl;
		//incidenceMatrix = incidenceMatrix + cellincidence;
		//*errorfile << std::setprecision(3) << "d2Adridrj = \n" << d2Adridrj_percell << "\n" << "\n" << endl;


		//*errorfile << std::setprecision(3) << "d2Pdridrj = \n" << d2Pdridrj_percell << "\n" << "\n" << endl;

		//*errorfile << std::setprecision(3) << "dPdridPdrj = \n" << dPdridPdrj_percell << "\n" << "\n" << endl;

		//*errorfile << std::setprecision(3) << "dAdridAdrj = \n" << dAdridAdrj_percell << "\n" << "\n" << endl;

		d2Adridrj = d2Adridrj + d2Adridrj_percell;
		d2Pdridrj = d2Pdridrj + d2Pdridrj_percell;
		dAdridAdrj = dAdridAdrj + dAdridAdrj_percell;
		dPdridPdrj = dPdridPdrj + dPdridPdrj_percell;
	}
	Hessian = d2Adridrj + d2Pdridrj + dAdridAdrj + dPdridPdrj;
	//*errorfile << std::setprecision(3) << "H = \n" << Hessian << "\n" << "\n" << endl;

	*errorfile << std::setprecision(6) << "incidence at the end  = \n" << incidenceMatrix << "\n" << "\n" << endl;

	*errorfile << std::setprecision(6) << "d2Adridrj = \n" << d2Adridrj << "\n" << "\n" << endl;

	*errorfile << std::setprecision(6) << "d2Pdridrj = \n" << d2Pdridrj << "\n" << "\n" << endl;

	*errorfile << std::setprecision(6) << "dPdridPdrj = \n" << dPdridPdrj << "\n" << "\n" << endl;

	*errorfile << std::setprecision(6) << "dAdridAdrj = \n" << dAdridAdrj << "\n" << "\n" << endl;

	//*errorfile << std::setprecision(6) << "dPdridPdrj - normed dAdridAdrj = \n" << dPdridPdrj - dAdridAdrj / 952.8286 << "\n" << "\n" << endl;
	//*errorfile << std::setprecision(6) << "Hessian = \n" << Hessian << "\n" << "\n" << endl;

	return Hessian;
}
/*
* double calculateShearModulusWithEigenDecomposition(errorfile)
*
* here, I will implement the shear modulus calculation with the theory developed in Merkel and Manning's 3D vertex model paper
* (A geometrically controlled rigidity transition in a model for confluent 3D tissues)
*
*
*/
double AGPerimTissue::calculateShearModulusWithEigenDecomposition(ofstream* errorfile)
{
	//double shearmodulus = 0;
	double area, perimeter, p0, a0, cell_pressure, tension_magnitude;
	int vmapid;
	vector<Vertex*> cwverts;
	unordered_map<int, int> vertex_index_map;
	int numofvertices = vertices.size();
	int numofcells = cells.size();

	MatrixXf Hessian(2 * numofvertices, 2 * numofvertices);
	Hessian = calculateHessian(errorfile);

	SelfAdjointEigenSolver<MatrixXf> hessianeigen(Hessian);
	//MatrixXf hessianeigenvalues(2 * numofvertices, 1);
	//hessianeigenvalues = eigensolver.eigenvalues();

	FullPivLU<MatrixXf> lu_decomp(Hessian);

	//MatrixXf hessianeigenvectors(2 * numofvertices, 2 * numofvertices);
	//hessianeigenvectors = eigensolver.eigenvectors();

	MatrixXf BoundaryVector(2 * numofvertices, 1);
	BoundaryVector.setZero(2 * numofvertices, 1);

	MatrixXf d2Edgamma2_mat(1, 1);
	d2Edgamma2_mat.setZero(1, 1);

	//MatrixXf shearModulus(1, 1);
	//shearModulus.setZero(1, 1);

	double shearmodulus = 0;
	double d2Edgamma2 = 0;

	Vect d2Pdrdgamma;
	Vect dAdv, dPdv;
	double d2Adrdgamma;

	double dadgcheck = 0;
	double dpdgcheck = 0;
	double dpdpdgdgcheck = 0;
	int cellid = 0;
	int relid = 0;
	int index = 0;
	int currenteigenvalue;
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
	{
		// The fastest way to do this is with an unordered map (the C++ equivalent to Python dictionaries)
		vertex_index_map[(*(*vit)).id] = index;
		index++;
	}

	vector<Cell*> xboundarycrossingcell = (*this).GetXCrossingCellVector();
	for (vector<Cell*>::iterator cit = xboundarycrossingcell.begin(); cit != xboundarycrossingcell.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		cellid = (*(*cit)).id;
		*errorfile << "Cell id: " << cellid << endl;
		area = ((*agit)).GetArea();
		perimeter = ((*agit)).GetPerimeter();
		p0 = (*agit).P0;
		a0 = ((*agit)).A0;
		cwverts = (*(*cit)).SavedInfoForEvolve_CWVerts;	//Get the cw vertices in the cell from the preassigned list 
		//cwverts.push_back(cwverts[0]); // we'll make the vertex list extra long so that it loops around.
		// We effectively have 4 terms we need in the Hessian.

		//K = 1; // this is for debugging
		//LAMBDA = 1; // this is for debugging, just as in the Hessian calculation

		// the pressure
		cell_pressure = (*(*cit)).GetPressure(area) * a0; // I want pressure to be PRESSURE = K (A - A0)
		// tension magnitude should be simple enough
		tension_magnitude = LAMBDA * (perimeter - p0);


		dadgcheck = (*agit).GetdAdgammaForSimpleShearAlongXAxis();
		dpdgcheck = (*agit).GetdPdgammaForSimpleShearAlongXAxis();
		dpdpdgdgcheck = (*agit).Getd2Pdgamma2ForSimpleShearAlongXAxis();
		// let's do d2Edgamma2 first:
		d2Edgamma2 += K * dadgcheck * dadgcheck;
		d2Edgamma2 += tension_magnitude * dpdpdgdgcheck;
		d2Edgamma2 += LAMBDA * dpdgcheck * dpdgcheck;

		for (vector<Vertex*>::iterator vit = cwverts.begin(); vit != cwverts.end(); vit++)
		{
			relid = (*(*vit)).id;
			vmapid = vertex_index_map[relid];
			d2Adrdgamma = (*agit).Getd2AdgammadyForSimpleShearAlongXAxis((*vit));
			BoundaryVector(2 * vmapid + 1, 0) += cell_pressure * d2Adrdgamma;
			d2Pdrdgamma = (*agit).Getd2PdgammadvForSimpleShearAlongXAxis((*vit));
			BoundaryVector(2 * vmapid, 0) += tension_magnitude * d2Pdrdgamma.x;
			BoundaryVector(2 * vmapid + 1, 0) += tension_magnitude * d2Pdrdgamma.y;
			dAdv = (*agit).GetdAdv((*vit));
			BoundaryVector(2 * vmapid, 0) += K * dadgcheck * dAdv.x; // double check dAdv
			BoundaryVector(2 * vmapid + 1, 0) += K * dadgcheck * dAdv.y; // double check dAdv
			dPdv = (*agit).GetdPdv((*vit));
			BoundaryVector(2 * vmapid, 0) += LAMBDA * dpdgcheck * dPdv.x; // double check dPdv
			BoundaryVector(2 * vmapid + 1, 0) += LAMBDA * dpdgcheck * dPdv.y; // double check dPdv
		}
	}
	d2Edgamma2_mat(0, 0) = d2Edgamma2;
	

	//following code is for calculating the shear modulus with Manning and Merkel's method (drops zero modes)
	* errorfile << "\n \n \n \n Boundary Vector = " << BoundaryVector << "\n \n \n \n " << endl;
	double nonaffine_shearportion = 0;
	//double dispx, dispy, rot;
	*errorfile << "Number of eigenvalues = " << hessianeigen.eigenvalues().size() << "\n \n \n" << endl;
	for (int index = 0; index < hessianeigen.eigenvalues().size(); index++)
	{
		int jndex = 0;
		//dispx = 0;
		//dispy = 0;
		//rot = 0;
		currenteigenvalue = real(hessianeigen.eigenvalues()(index, 0));
		*errorfile << hessianeigen.eigenvalues()(index, 0) << " , eigenvalue " << index << " and the vect : \n";
		// debug code for eigenvalues
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
		{
			*errorfile << "Vertex " << (*(*vit)).id << "x " << hessianeigen.eigenvectors().col(index)(jndex) << " Bvx " << BoundaryVector(jndex, 0) << endl;
			//dispx += hessianeigen.eigenvectors().col(index)(jndex) / hessianeigen.eigenvalues().size();
			//rot -= (*(*vit)).r.y * hessianeigen.eigenvectors().col(index)(jndex) / hessianeigen.eigenvalues().size();
			jndex++; 
			*errorfile << "Vertex " << (*(*vit)).id << "y " << hessianeigen.eigenvectors().col(index)(jndex) << " Bvy " << BoundaryVector(jndex, 0) << endl;
			//dispy += hessianeigen.eigenvectors().col(index)(jndex) / hessianeigen.eigenvalues().size();
			//rot += (*(*vit)).r.x * hessianeigen.eigenvectors().col(index)(jndex) / hessianeigen.eigenvalues().size();
			jndex++;
		}
		//errorfile << " \n \n COM displacement =(" << dispx << ", " << dispy << ")      ROTATION = " << rot << "\n \n" << endl;
		
		if ((index > lu_decomp.rank() - 1)) // HOW TO IDENTIFY ZERO MODES?
		{
			nonaffine_shearportion += 1 / hessianeigen.eigenvalues()(index, 0) * pow((hessianeigen.eigenvectors().col(index).transpose() * BoundaryVector)(0, 0), 2);
			//nonaffine_shearportion = hessianeigenvalues(index, 0);
		}
	}
	*errorfile << "\n \n \n The column of the Hessian is the number of vertices" << " - " << hessianeigen.eigenvalues().size() - lu_decomp.rank() << endl;
	*errorfile << "Here is a matrix whose columns form a basis of the null-space of A:\n"
		<< lu_decomp.kernel() <<"\n \n" << endl;
	//*errorfile << "Eigenvalues " << hessianeigenvalues << endl;
	shearmodulus = 1 / (K * A0 * A0 * numofcells * LX * LY) * (d2Edgamma2 - nonaffine_shearportion);//numofcells / (K * pow((LX * LY),2)) * (d2Edgamma2 - nonaffine_shearportion);
	//shearmodulus = 1 / (K * A0 * LX * LY) * (d2Edgamma2 - nonaffine_shearportion);
	return shearmodulus;

}

/*
* double CalculateShearModulusWithMatrixInversion(errorfile)
* 
*  NOT USABLE CURRENTLY (11/30/2023)
*
* here, I will implement the shear modulus calculation with the theory developed in Dapeng (Max) Bi's Multicellular Rosette's paper
* (Multicellular Rosettes Drive Fluid-Solid Transition in Epithelial Tissue)
* 
* NEEDS TO BE DOUBLE CHECKED ON A PRESTRESSED NETWORK
*/
double AGPerimTissue::CalculateShearModulusWithMatrixInversion(ofstream* errorfile)
{
	//double shearmodulus = 0;
	double area, perimeter, p0, a0, cell_pressure, tension_magnitude, K, LAMBDA;
	int vmapid;
	vector<Vertex*> cwverts;
	unordered_map<int, int> vertex_index_map;
	int numofvertices = vertices.size();
	
	MatrixXf Hessian(2 * numofvertices, 2 * numofvertices);
	Hessian = calculateHessian(errorfile);

	MatrixXf BoundaryVector(2 * numofvertices, 1);
	BoundaryVector.setZero(2 * numofvertices, 1);

	MatrixXf d2Edgamma2_mat(1, 1);
	d2Edgamma2_mat.setZero(1, 1);

	MatrixXf shearModulus(1, 1);
	shearModulus.setZero(1, 1);

	double d2Edgamma2 = 0;

	Vect dPdrdgamma;
	Vect dAdv, dPdv;
	
	int index = 0;
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
	{
		// The fastest way to do this is with an unordered map (the C++ equivalent to Python dictionaries)
		vertex_index_map[(*(*vit)).id] = index;
		index++;
	}

	vector<Cell*> xboundarycrossingcell = (*this).GetXCrossingCellVector();
	for (vector<Cell*>::iterator cit = xboundarycrossingcell.begin(); cit != xboundarycrossingcell.end(); cit++)
	{
		AGPerimCell* agit = dynamic_cast<AGPerimCell*>(*cit);
		*errorfile << "Cell id: " << (*(*cit)).id << endl;
		area = ((*agit)).GetArea();
		perimeter = ((*agit)).GetPerimeter();
		p0 = (*agit).P0;
		a0 = ((*agit)).A0;
		cwverts = (*(*cit)).SavedInfoForEvolve_CWVerts;	//Get the cw vertices in the cell from the preassigned list 
		//cwverts.push_back(cwverts[0]); // we'll make the vertex list extra long so that it loops around.
		// We effectively have 4 terms we need in the Hessian.

		K = 1; // this is for debugging
		LAMBDA = 1; // this is for debugging, just as in the Hessian calculation

		// the pressure
		cell_pressure = 1;//(*(*cit)).GetPressure(area) * a0; // I want pressure to be PRESSURE = K (A - A0)
		// tension magnitude should be simple enough
		tension_magnitude = 1; //LAMBDA * (perimeter - p0);

		

		// let's do d2Edgamma2 first:
		d2Edgamma2 += K * (*agit).GetdAdgammaForSimpleShearAlongXAxis() * (*agit).GetdAdgammaForSimpleShearAlongXAxis();
		d2Edgamma2 += tension_magnitude * (*agit).Getd2Pdgamma2ForSimpleShearAlongXAxis();
		d2Edgamma2 += LAMBDA * (*agit).GetdPdgammaForSimpleShearAlongXAxis() * (*agit).GetdPdgammaForSimpleShearAlongXAxis();

		for (vector<Vertex*>::iterator vit = cwverts.begin(); vit != cwverts.end(); vit++)
		{
			vmapid = vertex_index_map[(*(*vit)).id];
			//BoundaryVector(2 * vmapid)
			BoundaryVector(2 * vmapid + 1, 0) += cell_pressure * (*agit).Getd2AdgammadyForSimpleShearAlongXAxis((*vit));
			dPdrdgamma = (*agit).Getd2PdgammadvForSimpleShearAlongXAxis((*vit));
			BoundaryVector(2 * vmapid, 0) += tension_magnitude * dPdrdgamma.x;
			BoundaryVector(2 * vmapid + 1, 0) += tension_magnitude * dPdrdgamma.y;
			dAdv = (*agit).GetdAdv((*vit));
			BoundaryVector(2 * vmapid, 0) += K * (*agit).GetdAdgammaForSimpleShearAlongXAxis() * dAdv.x; // double check dAdv
			BoundaryVector(2 * vmapid + 1, 0) += K * (*agit).GetdAdgammaForSimpleShearAlongXAxis() * dAdv.y; // double check dAdv
			dPdv = (*agit).GetdPdv((*vit));
			BoundaryVector(2 * vmapid, 0) += LAMBDA * (*agit).GetdPdgammaForSimpleShearAlongXAxis() * dPdv.x; // double check dPdv
			BoundaryVector(2 * vmapid + 1, 0) += LAMBDA * (*agit).GetdPdgammaForSimpleShearAlongXAxis() * dPdv.y; // double check dPdv
		}
	}
	d2Edgamma2_mat(0, 0) = d2Edgamma2;
	shearModulus= 1 / numofvertices * (d2Edgamma2_mat - BoundaryVector.transpose() * Hessian.inverse() * BoundaryVector);
	*errorfile << "Boundary Vector = " << BoundaryVector << "\n \n \n \n " << endl;
	*errorfile << "Hessian = " << Hessian <<"\n \n \n \n " << endl;
	*errorfile << "Hessian det = " << Hessian.determinant() << "\n \n \n \n " << endl;
	*errorfile << "Hessian times its Inverse = " << Hessian*Hessian.inverse() << "\n \n \n \n " << endl;
	*errorfile << shearModulus << "\n \n \n \n" << endl;
	*errorfile << d2Edgamma2_mat << "\n \n \n \n" << endl;
	return shearModulus(0, 0);
}
/*
* double listAverage(inputList)
*
* function to find the average of a (c++ stl) list - used in Evolve to determine equilibrium
*
*/
double listAverage(list<double> inputList)
{
	int listSize = inputList.size();
	double average = 0;
	for (list<double>::iterator lit = inputList.begin(); lit != inputList.end(); lit++)
	{
		average += *(lit);
	}
	average = average / listSize;
	return average;
}
/*
* double getLinearRegressionSlope(inputList)
*
* function to find the slope of a list of points, stored here in C++ stl lists - used in Evolve to determine equilibrium
*
*/
double getLinearRegressionSlope(list<double> x, list<double> y)
{
	double xbar = listAverage(x);
	double ybar = listAverage(y);
	double num = 0;
	double den = 0; // slope of linear regression
	for (int index = 0; index < x.size(); index++)
	{
		num += (x.front() - xbar) * (y.front() - ybar);
		den += (x.front() - xbar) * (x.front() - xbar);
		x.pop_front();
		y.pop_front();
	}
	return num / den;
}



