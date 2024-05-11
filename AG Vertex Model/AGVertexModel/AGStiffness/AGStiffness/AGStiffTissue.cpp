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
#include "VectClass.h"
#include "AGStiffCell.h"
#include "AGStiffTissue.h"

using namespace Eigen;
using namespace std;


/*
* bool CalculateFreeEnergy()
*
* Calculates the energy of the packing using the Hamiltonian
* in this model:
* H = (sum over cells) [ 0.5 K(x,y) (A_cells - A_0) + 0.5 LAMBDA (P_cell - P_0)^2]
*/
double AGStiffTissue::CalculateFreeEnergy()
{
	double carea, cperim;
	double a0, p0;
	double freeenergy = 0;

	double perim_difference_sqr = 0;
	// (P_c - P_0)^2
	double pressure_portion = 0;
	// (A_c - A_0)^2

	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGStiffCell* agit = dynamic_cast<AGStiffCell*>(*cit);
		carea = ((*agit)).GetArea();
		cperim = ((*agit)).GetPerimeter();
		p0 = (*agit).P0;
		a0 = ((*agit)).A0;
		
		pressure_portion += (*agit).KAM * a0 * a0 * pow(carea / a0 - 1, 2);
	
		perim_difference_sqr += pow(cperim - p0, 2);

		// extra factor of 1/A0 is because K is written in units of pressure instead of E/length^4
	}
	freeenergy = 0.5 * LAMBDA * perim_difference_sqr + 0.5 * pressure_portion - (*this).CalculateBoxEnergy();


	return freeenergy;
}

/*
* double CalculatePressureEnergy()
*
* Calculates the area energy of the packing using the Hamiltonian
* in this model:
* H = (sum over cells) [ 0.5 K(A_cells - A_0)^2 ]
*/
double AGStiffTissue::CalculatePressureEnergy()
{
	double freeenergy = 0;

	double pressure_portion = 0;
	// (A_c - A_0)^2
	double carea;
	double a0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGStiffCell* agit = dynamic_cast<AGStiffCell*>(*cit);
		carea = (*agit).GetArea();
		a0 = (*agit).A0;
		
		pressure_portion += (*agit).KAM * a0 * a0 * pow(carea / a0 - 1, 2);
			// where the energy would look like Ep = 1/2 K A0^2 (A/A0 - 1)^2 ([K] = pressure)
		
	}
	freeenergy = 0.5 * pressure_portion;
	return freeenergy;
}

/*
* void AssignAreaModuli()
*
* Sets the values of KAM for the cells in the tissue
*/
void AGStiffTissue::AssignAreaModuli()
{
	
	double xcom = 0;
	double stifffunc = 0;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++)
	{
		AGStiffCell* agit = dynamic_cast<AGStiffCell*>(*cit);
		xcom = (*agit).GetCenter().x;
		
		stifffunc = K + 0.2 * K * sin(xcom * PI / LX); // FOR NOW, THIS IS THE CHOSEN STIFFNESS PATTERN
		(*agit).SetAreaModulus(stifffunc);

	}
}

/*
* void HelperReadInputData_SetGlobalValues(ifstream* datafile)
*
* Reads the data and sets all of the global variables.
* Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
*/
void AGStiffTissue::HelperReadInputData_SetGlobalValues(ifstream* datafile) {

	string line, word, key, value;		//temporary variables to store lines and words 
	string type;						//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline;				//ssline is a string stream so we can comvert strings to stringstreams
	map <string, string> globals;



	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//

	//Read until you get to the word Global signifying the begining of the global variables
	while (line.compare("Global:") != 0) {
		getline(*datafile, line);
		cout << line << endl;
	}

	//Now read the global variables and make a map with the keywords
	while (line.compare("endGlobal") != 0) {
		getline(*datafile, line);							//get the next line
		ssline.clear();
		ssline << line;										//put it into a stream
		ssline >> key;										//the first word is the key
		ssline >> value;									//the second is the value of the key
		if (value != "") {									//if the item actually has a value
			globals[key] = value;								//stick the item and it's value in the globals map
		}
	}

	//Now assign the constants their values
	//makes sure the key is actually part of the map before giving the constant a value if the key is not defined it will have the value given at the begining of main
	// the atof is there to convert the string to a double
	if (globals.count("LX") == 1) { LX = atof(globals["LX"].c_str()); }
	else if (globals.count("LX") == 0) { cout << "WARNING: LX is not given in input file, used default value " << LX << endl; }
	else { cout << "ERROR: Global variable LX has multiple definition in input file, used default value LX = " << LX << endl; }

	if (globals.count("LY") == 1) { LY = atof(globals["LY"].c_str()); }
	else if (globals.count("LY") == 0) { cout << "WARNING: LY is not given in input file, used default value " << LY << endl; }
	else { cout << "ERROR: Global variable LX has multiple definition in input file, used default value LY = " << LY << endl; }

	if (globals.count("MINAREA") == 1) { MINAREA = atof(globals["MINAREA"].c_str()); }
	else if (globals.count("MINAREA") == 0) { cout << "WARNING: MINAREA is not given in input file, used default value " << MINAREA << endl; }
	else { cout << "ERROR: Global variable MINAREA has multiple definition in input file, used default value MINAREA = " << MINAREA << endl; }

	if (globals.count("MINLENGTH") == 1) { MINLENGTH = atof(globals["MINLENGTH"].c_str()); }
	else if (globals.count("MINLENGTH") == 0) { cout << "WARNING: LX is not given in input file, used default value " << MINLENGTH << endl; }
	else { cout << "ERROR: Global variable MINLENGTH has multiple definition in input file, used default value MINLENGTH = " << MINLENGTH << endl; }

	if (globals.count("NEWEDGELENGTH") == 1) { NEWEDGELENGTH = atof(globals["NEWEDGELENGTH"].c_str()); }
	else if (globals.count("NEWEDGELENGTH") == 0) { cout << "WARNING: NEWEDGELENGTH is not given in input file, used default value " << NEWEDGELENGTH << endl; }
	else { cout << "ERROR: Global variable NEWEDGELENGTH has multiple definition in input file, used default value NEWEDGELENGTH = " << NEWEDGELENGTH << endl; }

	if (globals.count("GAMMA") == 1) { GAMMA = atof(globals["GAMMA"].c_str()); }
	else if (globals.count("GAMMA") == 0) { cout << "WARNING: GAMMA is not given in input file, used default value " << GAMMA << endl; }
	else { cout << "ERROR: Global variable GAMMA has multiple definition in input file, used default value GAMMA = " << GAMMA << endl; }

	if (globals.count("K") == 1) { K = atof(globals["K"].c_str()); }
	else if (globals.count("K") == 0) { cout << "WARNING: K is not given in input file, used default value " << K << endl; }
	else { cout << "ERROR: Global variable K has multiple definition in input file, used default value K = " << K << endl; }

	if (globals.count("SEED") == 1) { SEED = atof(globals["SEED"].c_str()); }
	else if (globals.count("SEED") == 0) { cout << "WARNING: SEED is not given in input file, used default value " << SEED << endl; }
	else { cout << "ERROR: Global variable SEED has multiple definition in input file, used default value SEED = " << SEED << endl; }

	if (globals.count("INITIALSTEPSIZE") == 1) { INITIALSTEPSIZE = atof(globals["INITIALSTEPSIZE"].c_str()); }
	else if (globals.count("INITIALSTEPSIZE") == 0) { cout << "WARNING: INITIALSTEPSIZE is not given in input file, used default value " << INITIALSTEPSIZE << endl; }
	else { cout << "ERROR: Global variable INITIALSTEPSIZE has multiple definition in input file, used default value INITIALSTEPSIZE = " << INITIALSTEPSIZE << endl; }

	if (globals.count("TIMEINTERVAL") == 1) { TIMEINTERVAL = atof(globals["TIMEINTERVAL"].c_str()); }
	else if (globals.count("TIMEINTERVAL") == 0) { cout << "WARNING: TIMEINTERVAL is not given in input file, used default value " << TIMEINTERVAL << endl; }
	else { cout << "ERROR: Global variable TIMEINTERVAL has multiple definition in input file, used default value TIMEINTERVAL = " << TIMEINTERVAL << endl; }

	if (globals.count("TOTALTIME") == 1) { TOTALTIME = atof(globals["TOTALTIME"].c_str()); }
	else if (globals.count("TOTALTIME") == 0) { cout << "WARNING: TOTALTIME is not given in input file, used default value " << TOTALTIME << endl; }
	else { cout << "ERROR: Global variable TOTALTIME has multiple definition in input file, used default value TOTALTIME = " << TOTALTIME << endl; }

	if (globals.count("CURRENTTIME") == 1) { CURRENTTIME = atof(globals["CURRENTTIME"].c_str()); }
	else if (globals.count("CURRENTTIME") == 0) { cout << "WARNING: CURRENTTIME is not given in input file, used default value " << CURRENTTIME << endl; }
	else { cout << "ERROR: Global variable CURRENTTIME has multiple definition in input file, used default value CURRENTTIME = " << CURRENTTIME << endl; }

	if (globals.count("error_tolerance_absolute_mechanical") == 1) { error_tolerance_absolute_mechanical = atof(globals["error_tolerance_absolute_mechanical"].c_str()); }
	else if (globals.count("error_tolerance_absolute_mechanical") == 0) { cout << "WARNING: error_tolerance_absolute_mechanical is not given in input file, used default value " << error_tolerance_absolute_mechanical << endl; }
	else { cout << "ERROR: Global variable error_tolerance_absolute_mechanical has multiple definition in input file, used default value error_tolerance_absolute_mechanical = " << error_tolerance_absolute_mechanical << endl; }

	if (globals.count("ALPHAM") == 1) { alphaM = atof(globals["ALPHAM"].c_str()); }
	else if (globals.count("ALPHAM") == 0) { cout << "WARNING: ALPHAM is not given in input file, used default value " << alphaM << endl; }
	else { cout << "ERROR: Global variable ALPHAM has multiple definition in input file, used default value ALPHAM = " << alphaM << endl; }

	if (globals.count("TISSUE_ALPHA") == 1) { TISSUE_ALPHA = atof(globals["TISSUE_ALPHA"].c_str()); }
	else if (globals.count("TISSUE_ALPHA") == 0) { cout << "WARNING: TISSUE_ALPHA is not given in input file, used default value " << TISSUE_ALPHA << endl; }
	else { cout << "ERROR: Global variable TISSUE_ALPHA has multiple definition in input file, used default value TISSUE_ALPHA = " << TISSUE_ALPHA << endl; }

	if (globals.count("TISSUE_EXT_STRESS") == 1) { TISSUE_EXT_STRESS = atof(globals["TISSUE_EXT_STRESS"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS") == 0) { cout << "WARNING: TISSUE_EXT_STRESS is not given in input file, used default value " << TISSUE_EXT_STRESS << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS has multiple definition in input file, used default value TISSUE_EXT_STRESS = " << TISSUE_EXT_STRESS << endl; }

	if (globals.count("TISSUE_EXT_STRESS_X") == 1) { TISSUE_EXT_STRESS_X = atof(globals["TISSUE_EXT_STRESS_X"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_X") == 0) { cout << "WARNING: TISSUE_EXT_STRESS_X is not given in input file, used default value " << TISSUE_EXT_STRESS_X << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_X has multiple definition in input file, used default value TISSUE_EXT_STRESS_X = " << TISSUE_EXT_STRESS_X << endl; }

	if (globals.count("TISSUE_EXT_STRESS_Y") == 1) { TISSUE_EXT_STRESS_Y = atof(globals["TISSUE_EXT_STRESS_Y"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_Y") == 0) { cout << "WARNING: TISSUE_EXT_STRESS_Y is not given in input file, used default value " << TISSUE_EXT_STRESS_Y << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_Y has multiple definition in input file, used default value TISSUE_EXT_STRESS_Y = " << TISSUE_EXT_STRESS_Y << endl; }

	if (globals.count("MITOSISON") == 1) { MITOSISON = (atof(globals["MITOSISON"].c_str()) > 0); }
	else if (globals.count("MITOSISON") == 0) { cout << "WARNING: MITOSISON is not given in input file, used default value " << MITOSISON << endl; }
	else { cout << "ERROR: Global variable MITOSISON has multiple definition in input file, used default value MITOSISON = " << MITOSISON << endl; }

	if (globals.count("T1ON") == 1) { T1ON = (atof(globals["T1ON"].c_str()) > 0); }
	else if (globals.count("T1ON") == 0) { cout << "WARNING: T1ON is not given in input file, used default value " << T1ON << endl; }
	else { cout << "ERROR: Global variable T1ON has multiple definition in input file, used default value T1ON = " << T1ON << endl; }

	if (globals.count("APOPTOSISON") == 1) { APOPTOSISON = (atof(globals["APOPTOSISON"].c_str()) > 0); }
	else if (globals.count("APOPTOSISON") == 0) { cout << "WARNING: APOPTOSISON is not given in input file, used default value " << APOPTOSISON << endl; }
	else { cout << "ERROR: Global variable APOPTOSISON has multiple definition in input file, used default value APOPTOSISON = " << APOPTOSISON << endl; }

	if (globals.count("WRITEMITOSIS") == 1) { WRITEMITOSIS = atof(globals["WRITEMITOSIS"].c_str()); }
	else if (globals.count("WRITEMITOSIS") == 0) { cout << "WARNING: WRITEMITOSIS is not given in input file, used default value " << WRITEMITOSIS << endl; }
	else { cout << "ERROR: Global variable WRITEMITOSIS has multiple definition in input file, used default value WRITEMITSOSIS = " << WRITEMITOSIS << endl; }

	if (globals.count("VERBOSEON") == 1) { VERBOSEON = (atof(globals["VERBOSEON"].c_str()) > 0); }
	else if (globals.count("VERBOSEON") == 0) { cout << "WARNING: VERBOSEON is not given in input file, used default value " << VERBOSEON << endl; }
	else { cout << "ERROR: Global variable VERBOSEON has multiple definition in input file, used default value APOPTOSISON = " << VERBOSEON << endl; }

	if (globals.count("EQUILIBRIUMMITOSIS") == 1) { EQUILIBRIUMMITOSIS = (atof(globals["EQUILIBRIUMMITOSIS"].c_str()) > 0); }
	else if (globals.count("EQUILIBRIUMMITOSIS") == 0) { cout << "WARNING: EQUILIBRIUMMITOSIS is not given in input file, used default value " << EQUILIBRIUMMITOSIS << endl; }
	else { cout << "ERROR: Global variable EQUILIBRIUMMITOSIS has multiple definition in input file, used default value EQUILIBRIUMMITOSIS = " << EQUILIBRIUMMITOSIS << endl; }

	if (globals.count("P0") == 1) { P0 = (atof(globals["P0"].c_str())); }
	else if (globals.count("P0") == 0) { cout << "WARNING: P0 is not given in input file, used default value " << P0 << endl; }
	else { cout << "ERROR: Global variable P0 has multiple definition in input file, used default value P0 = " << P0 << endl; }

	if (globals.count("LAMBDA") == 1) { LAMBDA = (atof(globals["LAMBDA"].c_str())); }
	else if (globals.count("LAMBDA") == 0) { cout << "WARNING: LAMBDA is not given in input file, used default value " << LAMBDA << endl; }
	else { cout << "ERROR: Global variable LAMBDA has multiple definition in input file, used default value LAMBDA = " << LAMBDA << endl; }

	if (globals.count("A0") == 1) { A0 = (atof(globals["A0"].c_str())); }
	else if (globals.count("A0") == 0) { cout << "WARNING: A0 is not given in input file, used default value " << A0 << endl; }
	else { cout << "ERROR: Global variable A0 has multiple definition in input file, used default value A0 = " << A0 << endl; }

	if (globals.count("energyListSize") == 1) { energyListSize = (atof(globals["energyListSize"].c_str())); }
	else if (globals.count("energyListSize") == 0) { cout << "WARNING: energyListSize is not given in input file, used default value " << energyListSize << endl; }
	else { cout << "ERROR: Global variable energyListSize has multiple definition in input file, used default value energyListSize = " << energyListSize << endl; }

	if (globals.count("equilibriumtolerance") == 1) { equilibriumtolerance = (atof(globals["equilibriumtolerance"].c_str())); }
	else if (globals.count("equilibriumtolerance") == 0) { cout << "WARNING: equilibriumtolerance is not given in input file, used default value " << equilibriumtolerance << endl; }
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
void AGStiffTissue::HelperReadInputData_CreateObjects(ifstream* datafile) {

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
	while (line.compare("Vertices") != 0) {
		getline(*datafile, line);
	}
	//then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')) {
		keywords.push_back(key);
	}
	//now recorde the data and make the vertex object untill you get to the word edges signifying the end of the data
	while (true) {
		getline(*datafile, line);							//get the next line
		if ((line.compare("Edges") == 0)) { break; }
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')) {				//parse it and put the data into the objectInfo map	
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
	while (getline(ssline, key, '\t')) {
		keywords.push_back(key);
	}
	//now recorde the data and make the edge object untill you get to the word cells signifying the end of the data
	while (true) {
		getline(*datafile, line);							//get the next line
		if (line.compare("Cells") == 0) { break; }
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')) {				//parse it and put the data into the objectInfo map	
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
	while (getline(ssline, key, '\t')) {
		keywords.push_back(key);
	}
	//now recorde the data and make the cell object untill you get to the word end signifying the end of the data
	while (true) {
		getline(*datafile, line);							//get the next line
		if (line.compare("End") == 0 || line.compare("Fibers") == 0) { break; }				//break if keyworkd 'End' is found
		ssline.clear();
		ssline << line;
		objectInfo.clear();
		i = 0;
		while (getline(ssline, value, '\t')) {				//parse it and put the data into the objectInfo map	
			objectInfo[keywords[i]] = value;
			i++;
		}
		//extract the data from the map
		newid = atoi(objectInfo["id"].c_str());
		enumb = atoi(objectInfo["EdgeNumb"].c_str());
		if (objectInfo.count("Type")) {
			typenumb = atoi(objectInfo["Type"].c_str());
		}
		else { typenumb = 0; }
		if (objectInfo.count("A0")) {
			A0 = atof(objectInfo["A0"].c_str());
		}
		// you only want atoi if you have a number you want to store as a string. but we should probably have a string


		ssedges.clear();
		ssedges << objectInfo["EdgeList"];
		//now worry about the edge list
		listedges.clear();
		while (getline(ssedges, edg, ',')) {
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
		AGStiffCell* cptr = new AGStiffCell;
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
		(*this).AssignAreaModuli();
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
	(*this).currentFreeEnergy = A0 * A0 * A0 * (K * A0 + LAMBDA * P0 * P0) * cells.size(); // this just needs to start out as a bigger number than the actual tissue energy
	(*this).pastFreeEnergy = (*this).currentFreeEnergy;
}

/*
* void HelperWriteOutputData_Write2File
*
* Writes data into file specified by 'datafile'. If flag 'includeExtraInfo' is true, then additional information (such as tension, cell area, etc) is included.
*The extra information is non-essential in reconstructing the simulation state, and is just provided for convenience (e.g. for display in MATLAB).
* (should be helpful for incremental coding erro debug sessions)
*
*/

void AGStiffTissue::HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo) {
	time_t rawtime;					//to deal with writing the date out
	struct tm* timeinfo;
	*datafile << "Global:" << endl;
	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;
	//Now write the global variables you used
	*datafile << std::setprecision(outputPrecision) << "CURRENTTIME\t" << CURRENTTIME << endl;
	*datafile << "SEED\t" << SEED << endl;
	*datafile << std::setprecision(outputPrecision) << "TIMEINTERVAL\t" << TIMEINTERVAL << endl;
	*datafile << std::setprecision(outputPrecision) << "LX\t" << LX << endl;
	*datafile << std::setprecision(outputPrecision) << "LY\t" << LY << endl;
	*datafile << std::setprecision(outputPrecision) << "Abox\t" << LX * LY << endl;
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
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++) {			//loop through all of the vertices in the tissue
			*datafile << std::setprecision(outputPrecision) << (*(*vit)).id << "\t" << (*(*vit)).r.x << "\t" << (*(*vit)).r.y << endl;						//write it's info to the file
		}

		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++) {			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			*datafile
				<< (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << endl;			//write its info to the file
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;
		*datafile << "id\tA0\tEdgeNumb\tEdgeList\tType" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++) {			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++) {				//make comma-separated list of edges
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
			*datafile << std::setprecision(outputPrecision) << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << L << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << "\t" << (*(*eit)).GetTension_Magnitude() << "\t" << (*(*(*eit)).cell1).id << "\t" << (*(*(*eit)).cell2).id << endl;
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

		* datafile << std::setprecision(outputPrecision) << "id\tA0\tK\tEdgeNumb\tType\txflag\tyflag\tArea\tPerimeter\tPressure\tCenterx\tCentery\tEdgeList\tCWVertices\tCellPressureType" << endl;
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
			AGStiffCell* agit = dynamic_cast<AGStiffCell*>(*cit);
			*datafile << std::setprecision(outputPrecision) << (*(*cit)).id << "\t" << (*agit).A0 << "\t" << (*agit).KAM << "\t" << (*(*cit)).edges.size() << "\t" << (*(*cit)).type << "\t" << (*(*cit)).GetXFlag() << "\t" << (*(*cit)).GetYFlag() << "\t"
				<< (*agit).GetArea() << "\t" << (*agit).GetPerimeter() << "\t" << (*agit).GetPressure((*(*cit)).GetArea()) << "\t" << centroid.x << "\t" << centroid.y << "\t" << ssEdgeList.str() << "\t" << ssCWVertices.str() << "\t" << (*agit).pressureType << endl;
		}
	}

	*datafile << "End" << endl << endl;
}
