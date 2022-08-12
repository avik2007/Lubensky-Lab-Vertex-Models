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


		
		(*this).currentFreeEnergy = (K*A0 + LAMBDA*P0*P0)*cells.size(); 
		(*this).pastFreeEnergy = (*this).currentFreeEnergy;

		// instantiate the free energy to be larger than the actual energy of the packing
		// just so that this doesn't mess with an adaptive stepsize changer (for finding
		// equilibrium)


		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
 
	}
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
			area_difference_sqr += -2 * a0 * log(1 + carea/a0);
			// where the energy would look like Ep = - K A0^2 Log(1 + A/A0) ([K] = pressure)
		}
		else if ((*agit).pressureType == "EXPONENTIAL")
		{
			area_difference_sqr += 2 * a0*exp(-carea / a0); 
			// where the energy would look like Ep = K A0^2 exp((A/A0-1)) ([K] = pressure)
		}
		else
		{
			area_difference_sqr += pow(a0,2)*pow(carea/a0 - 1, 2);
			// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
		}
		 
		
		perim_difference_sqr += pow( cperim - p0, 2); 
		
		// extra factor of 1/A0 is because K is written in units of pressure instead of E/length^4
	}
	freeenergy =0.5*LAMBDA*perim_difference_sqr + 0.5*K*area_difference_sqr - (*this).CalculateBoxEnergy();


	return freeenergy;
}

/*
* Edge* CreateEdge
*
* Creates a new AGPerimEdge and adds it to the list of edges and returns a pointer
*/
Edge* AGPerimTissue::CreateEdge(Edge* subclass, bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2){
	int largestID, newid;
	Edge newedge;
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
			return &(*(*eit));
			break;
		}
	}
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
	//double holding = 0;
	

	while (CURRENTTIME < endtime)
	{
		count++;
		//STEP 1: If the topology has changed get a list of the clockwise vertices in each cell
		if (reinitialize){																	//If their has been any topological change OR vertices have crossed a boundary 
			for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
				(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
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
				if (*stepsize < 0.0000001) // this could be turned into a parameter, the minumum stepsize
				//if (*stepsize <   0.000000000001)
				{
					breakVar = true;
					//*errorfile << "Stepsize got too small. Break" << endl;
					break;
					// This is a step I take so that the adaptive time stepper doesn't keep pushing the stepsize all the way down the machine precision 
				}
				//*errorfile << "Current FE \t Past FE \t Stepsize\t Energy Change \n";
				//*errorfile << std::setprecision(outputPrecision) << currentFreeEnergy << "\t" << pastFreeEnergy << "\t" << *stepsize << "\t" << currentFreeEnergy - pastFreeEnergy << endl;
				//*errorfile << "Energy Increase - change timestep and redo step\n";
				//*errorfile << "\n";
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

		*errorfile << std::setprecision(outputPrecision) << CURRENTTIME << "\t" << currentFreeEnergy << "\t" << LX << "\t" << LY << endl;

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
		

		/*
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
			else if (timeratio > 1 && timeratio > floortr && timeratio < 2 * INITIALSTEPSIZE + floortr)
				// Gives two timestep windows for mitosis
			{
				flag4 = (*this).Mitosis();
			}
			// THIS ACTIVATES MITOSIS!!!
			
			if (WRITEMITOSIS) (*this).WriteOutputDataAppend(errorfile);
			
		}*/
		if (flag1  || flag2 || flag3 || flag4) reinitialize = true;  //Determine if we need to reinitialize

		

		

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
				(*this).equilibrated = (*this).IsEquilibrated();
				*errorfile << "Equilibrium condition was met." << endl;
			}
		}
 
		if (breakVar)
		{
			(*this).equilibrated = true;
			break;
		}
	}

	cout << "TISSUE EVOLVED TO TIME " << CURRENTTIME << " in " << count << " steps" << endl;
	cout << "LX = " << LX << "  LY = " << LY << endl;
	
	// 
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
* bool CalculatePressureEnergy()
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
			area_difference_sqr += -2 * a0 * log(1 + carea / a0);
			// where the energy would look like Ep = - 1 / 2 K A0^2 Log(1 + A/A0) ([K] = pressure)
		}
		else if ((*agit).pressureType == "EXPONENTIAL")
		{
			area_difference_sqr += 2 * pow(a0,2) * exp(- carea / a0 );
			// where the energy would look like Ep = 1/2 K A0^2 exp(-A/A0) ([K] = pressure)
		}
		else
		{
			area_difference_sqr += pow(a0,2) * pow(carea / a0 - 1, 2);
			// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
		}
	}
	freeenergy = 0.5*K*area_difference_sqr;
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
			*datafile << std::setprecision(outputPrecision) << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << endl;			//write its info to the file
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
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag\ttension\t_cell1\t_cell2" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++)
		{			
			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			double L = (*(*eit)).GetLength();
			*datafile << std::setprecision(outputPrecision) << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << "\t" << (*(*eit)).GetTension_Magnitude() << "\t" << (*(*(*eit)).cell1).id << "\t" << (*(*(*eit)).cell2).id << endl;
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

