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

//Other files we need 
#include "Constants.h"
#include "ConstantsPerim.h"
#include "TissueClass.h"
#include "PerimTissue.h"
#include "PerimEdge.h"
#include "PerimCell.h"

// Reads the data and sets all of the global variables
void Tissue_Perim::HelperReadInputData_SetGlobalValues(ifstream* datafile){

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

	if (globals.count("P0") == 1){ P0 = atof(globals["P0"].c_str()); }
	else if (globals.count("P0") == 0){ cout << "WARNING: P0 is not given in input file, used default value " << P0 << endl; }
	else { cout << "ERROR: Global variable P0 has multiple definition in input file, used default value P0 = " << P0 << endl; }

	if (globals.count("LAMBDA") == 1){ LAMBDA = atof(globals["LAMBDA"].c_str()); }
	else if (globals.count("LAMBDA") == 0){ cout << "WARNING: LAMBDA is not given in input file, used default value " << LAMBDA << endl; }
	else { cout << "ERROR: Global variable LAMBDA has multiple definition in input file, used default value LAMBDA = " << LAMBDA << endl; }

	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//
}


// Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids.
// Starts at current position of *datafile and reads until keyword 'End'.
void Tissue_Perim::HelperReadInputData_CreateObjects(ifstream* datafile){

	string line, key, value, edg;			//temporary variables to store lines and words 
	string type;							//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline, ssedges;			//ssline is a string stream so we can comvert strings to stringstreams
	vector<string> keywords;				//vector of strings to hold the keywords
	vector<int> listedges;
	int i;
	map <string, string> objectInfo;		//map to hold the keywords and their data
	int newid, v1, v2, enumb, typenumb;					//temp varables to store the object info
	double x, y;
	bool xf, yf;
	double A0;

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
		Edge_Perim* eptr = new Edge_Perim;
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
		ssedges.clear();
		ssedges << objectInfo["EdgeList"];
		//now worry about the edge list
		listedges.clear();
		while (getline(ssedges, edg, ',')){
			listedges.push_back(atoi(edg.c_str()));
		}
		//make the cell object
		Cell_Perim* cptr = new Cell_Perim;
		(*cptr).id = newid;
		(*cptr).type = typenumb;
		(*cptr).A0 = A0;
		(*cptr).IOEdgeID = listedges;
		(*cptr).P0 = P0; 
		cells.push_back(cptr);
		cellIDs.push_back(newid);
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}
}