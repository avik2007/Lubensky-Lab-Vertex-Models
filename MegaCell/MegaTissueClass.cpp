/**	\file	MegaTissueClass.cpp
*	Implementation of MegaTissue class.
*	Contains function definitions associated with MegaTissue objects.
*	Depends on base class Tissue.
*/

#include <iostream>
#include <string>
#include <cmath>
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
#include <cassert>					//for error checking/assertion

#include "Constants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "MegaTissueClass.h"
#include "MegaCellClass.h"
#include "MegaEdgeClass.h"
#include "MegaCellConstants.h"

using namespace std;

// Creates a new edge and addes it to the list of edges, returns a pointer to that edge.
Edge* MegaTissue::CreateEdge(bool xflg, bool yflg, Vertex* v1, Vertex* v2, Cell* c1, Cell* c2){
	int largestID, newid;
	Edge newedge;
	//Get the id
	largestID = (*max_element(edgeIDs.begin(), edgeIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;								//Now make a new ID that is unique

	if (newid > 1000000000) { cout << "ERROR to many edges" << endl; }
	//create the vertex and add the vertex to the lists
	MegaEdge* eptr = new MegaEdge;
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

// Reads and constructs simulation from file. The file must be formatted like InputTemplate.txt. This function parses the input file then initializes all the objects and pointers.
void MegaTissue::ReadInputData(string InputFileName){
	ifstream datafile;
	datafile.open(InputFileName);											//open file
	if (!datafile.is_open()){											//check if open file successful
		cout << "ERROR cannot open file: " << InputFileName << endl;
	}

	(*this).HelperReadInputData_SetGlobalValues(&datafile);
	(*this).HelperReadInputData_CreateObjects(&datafile);
	(*this).HelperReadInputData_CreatePointers();
	(*this).HelperReadInputData_CheckTypes();

	datafile.close();
}

//! Helper function for checking the cell input conforms to expectations.
void MegaTissue::HelperReadInputData_CheckTypes(){
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
		if ((*(*cit)).type != 0 && (*(*cit)).type != center_cell_type && (*(*cit)).type != outer_cell_type){
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
void MegaTissue::HelperReadInputData_SetGlobalValues(ifstream* datafile){

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

	if (globals.count("ALPHAM") == 1){ alphaM = atof(globals["ALPHAM"].c_str()); }
	else if (globals.count("ALPHAM") == 0){ cout << "WARNING: ALPHAM is not given in input file, used default value " << alphaM << endl; }
	else { cout << "ERROR: Global variable ALPHAM has multiple definition in input file, used default value ALPHAM = " << alphaM << endl; }

	if (globals.count("INNERPRESSURE") == 1){ inner_pressure = atof(globals["INNERPRESSURE"].c_str()); }
	else if (globals.count("INNERPRESSURE") == 0){ cout << "WARNING: INNERPRESSURE is not given in input file, used default value " << inner_pressure << endl; }
	else { cout << "ERROR: Global variable INNERPRESSURE has multiple definition in input file, used default value INNERPRESSURE = " << inner_pressure << endl; }

	if (globals.count("OUTERPRESSURE") == 1){ outer_pressure = atof(globals["OUTERPRESSURE"].c_str()); }
	else if (globals.count("OUTERPRESSURE") == 0){ cout << "WARNING: OUTERPRESSURE is not given in input file, used default value " << outer_pressure << endl; }
	else { cout << "ERROR: Global variable OUTERPRESSURE has multiple definition in input file, used default value OUTERPRESSURE = " << outer_pressure << endl; }

	if (globals.count("FACTORINCREASEDTENSION") == 1){ factor_increased_tension_at_inner_boundary = atof(globals["FACTORINCREASEDTENSION"].c_str()); }
	else if (globals.count("FACTORINCREASEDTENSION") == 0){ cout << "WARNING: FACTORINCREASEDTENSION is not given in input file, used default value " << factor_increased_tension_at_inner_boundary << endl; }
	else { cout << "ERROR: Global variable FACTORINCREASEDTENSION has multiple definition in input file, used default value FACTORINCREASEDTENSION = " << factor_increased_tension_at_inner_boundary << endl; }

	if (globals.count("ISVORONOIPACKING") == 1){ is_voronoi_packing = atof(globals["ISVORONOIPACKING"].c_str()); }
	else if (globals.count("ISVORONOIPACKING") == 0){ cout << "WARNING: ISVORONOIPACKING is not given in input file, used default value " << is_voronoi_packing << endl; }
	else { cout << "ERROR: Global variable ISVORONOIPACKING has multiple definition in input file, used default value ISVORONOIPACKING = " << is_voronoi_packing << endl; }

	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//
}

// Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids.
// Starts at current position of *datafile and reads until keyword 'End'.
void MegaTissue::HelperReadInputData_CreateObjects(ifstream* datafile){

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
		MegaEdge* eptr = new MegaEdge;
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
		MegaCell* cptr = new MegaCell;
		(*cptr).id = newid;
		(*cptr).type = typenumb;
		(*cptr).A0 = A0;
		(*cptr).IOEdgeID = listedges;
		cells.push_back(cptr);
		cellIDs.push_back(newid);
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}
}

// Writes data into file specified by 'datafile'. If flag 'includeExtraInfo' is true, then additional information (such as tension, cell area, etc) is included.
// The extra information is non-essential in reconstructing the simulation state, and is just provided for convenience (e.g. for display in MATLAB).
void MegaTissue::HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo){
	time_t rawtime;					//to deal with writing the date out
	struct tm * timeinfo;


	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;

	//Now write the global variables you used
	*datafile << "Global:" << endl;
	*datafile << "SEED\t" << SEED << endl;
	*datafile << "TIMEINTERVAL\t" << TIMEINTERVAL << endl;
	*datafile << "LX\t" << LX << endl;
	*datafile << "LY\t" << LY << endl;
	*datafile << "MINAREA\t" << MINAREA << endl;
	*datafile << "MINLENGTH\t" << MINLENGTH << endl;
	*datafile << "NEWEDGELENGTH\t" << NEWEDGELENGTH << endl;
	*datafile << "K\t" << K << endl;
	*datafile << "GAMMA\t" << GAMMA << endl;
	*datafile << "TOTALTIME\t" << TOTALTIME << endl;
	*datafile << "CURRENTTIME\t" << CURRENTTIME << endl;
	*datafile << "INITIALSTEPSIZE\t" << INITIALSTEPSIZE << endl;
	*datafile << "error_tolerance_absolute_mechanical\t" << error_tolerance_absolute_mechanical << endl;
	*datafile << "INNERPRESSURE\t" << inner_pressure << endl;
	*datafile << "OUTERPRESSURE\t" << outer_pressure << endl;
	*datafile << "FACTORINCREASEDTENSION\t" << factor_increased_tension_at_inner_boundary << endl;
	*datafile << "ISVORONOIPACKING\t" << is_voronoi_packing << endl;
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
		*datafile << "id\tA0\tEdgeNumb\tEdgeList\tType" << endl;
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){			//loop through all of the cells in the tissue
			stringstream ssEdgeList;
			ssEdgeList.clear();
			for (i = 0; i < (*(*cit)).edges.size(); i++){				//make comma-separated list of edges
				ssEdgeList << (*(*(*cit)).edges[i]).id << ",";
			}

			*datafile << (*(*cit)).id << "\t" << (*(*cit)).A0 << "\t" << (*(*cit)).edges.size() << "\t" << ssEdgeList.str() << "\t" << (*(*cit)).type << endl;
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

		//print the edge info
		*datafile << "Edges" << endl;
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag\t_length\t_cell1\t_cell2" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			double L = (*(*eit)).GetLength();
			*datafile << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << "\t" << L << "\t" << (*(*(*eit)).cell1).id << "\t" << (*(*(*eit)).cell2).id << endl;
		}

		//print the cell info
		int i;
		*datafile << "Cells" << endl;
		*datafile << "id\tA0\tEdgeNumb\tType\t_xflag\t_yflag\t_Area\t_Pressure\t_Center_x\t_Center_y\tEdgeList\t_CWVertices" << endl;
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

			Vect centroid = (*(*cit)).GetCenter();
			*datafile << (*(*cit)).id << "\t" << (*(*cit)).A0 << "\t" << (*(*cit)).edges.size() << "\t" << (*(*cit)).type << "\t" << (*(*cit)).GetXFlag() << "\t" << (*(*cit)).GetYFlag() << "\t"
				<< (*(*cit)).GetArea() << "\t" << (*(*cit)).GetPressure((*(*cit)).GetArea()) << "\t" << centroid.x << "\t" << centroid.y << "\t" << ssEdgeList.str() << "\t" << ssCWVertices.str() << endl;
		}
	}


	*datafile << "End" << endl << endl;
}

////WARNING!!! This function assumes that cell.SavedInfoForEvolve_CWVerts has been updated correctly
//Calculates the force on every vertex and stores the result in vertex.Tissue subclasses have their own implementation
void MegaTissue::ForceFunction(ofstream* errorfile){
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double area, deltaP;
	vector<vector<int>> cwflags;
	vector<Vertex*> cwverts;
	int i, cellnumber;
	Edge* eptr;
	double Pmag, angle;
	Vect FT;
	//Declare and initialize pointers.
	Vertex* v1 = 0;
	Vertex* v2 = 0;
	Cell* c1;
	Vect zerovec(0.0, 0.0);

	//First loop through all the cells and find the pressure WARNING!! cell.SavedInfoForEvolve_CWVerts must already be correctly assigned
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
		cwverts = (*(*cit)).SavedInfoForEvolve_CWVerts;							//Get the cw vertices in the cell from the preassigned list 
		area = (*(*cit)).GetArea(cwverts);										//Get the area quickly using the pre assigned cw vertices
		(*(*cit)).SavedInfoForEvolve_Pressure = (*(*cit)).GetPressure(area);	//Get the pressure
	}

	//Calculate the force on each vertex from the pressure and tension
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
		(*(*vit)).SavedInfoForEvolve_Force = zerovec; //first set the force to zero 
		//For each edge
		for (i = 0; i < (*(*vit)).edges.size(); i++){
			eptr = &(*(*(*vit)).edges[i]);
			//Get the force from the tension on the edge
			FT = (*eptr).GetTension_ForceOnVertex((*vit));
			//Get the force from the pressure difference across the edge
			v1 = (*eptr).GetOtherVertex((*vit));
			c1 = &(*(*eptr).cell1);
			v2 = (*c1).GetNextCWVert(*vit, (*c1).SavedInfoForEvolve_CWVerts);
			//Gets the CW ordering of cells 1 and 2 across eptr
			if ((*v1).id == (*v2).id){ deltaP = (*(*eptr).cell2).SavedInfoForEvolve_Pressure - (*(*eptr).cell1).SavedInfoForEvolve_Pressure; }
			else{ deltaP = (*(*eptr).cell1).SavedInfoForEvolve_Pressure - (*(*eptr).cell2).SavedInfoForEvolve_Pressure; }
			Pmag = -0.5 * deltaP * (*eptr).GetLength();   		//magnitude of the force due to the pressure difference
			angle = (*eptr).GetAngle((*vit)) + PI / 2;	 	//get the angle zXl
			Vect FP(zerovec, angle, Pmag);
			//Add the forces on this edge to the vertex
			(*(*vit)).SavedInfoForEvolve_Force = (*(*vit)).SavedInfoForEvolve_Force + FT + FP;
		}
	}

	Vect origin_this_cell;
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
		if ((*(*cit)).type == outer_cell_type){
			origin_this_cell = (*(*cit)).GetCenter();
		}
	}

	//Calculate the tissue stress
	double GammaX, GammaY;
	GammaX = 0.0;
	GammaY = 0.0;

	//Add the force from the tissue stress to every vertex
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
		(*(*vit)).SavedInfoForEvolve_Force.x = 1 / alphaM * (*(*vit)).SavedInfoForEvolve_Force.x + (1 / TISSUE_ALPHA) * (*(*vit)).r.x * GammaX;
		(*(*vit)).SavedInfoForEvolve_Force.y = 1 / alphaM * (*(*vit)).SavedInfoForEvolve_Force.y + (1 / TISSUE_ALPHA) * (*(*vit)).r.y * GammaY;
	}

	//Update LY and LX
	ForceLX = (1 / TISSUE_ALPHA) * GammaX * LX;
	ForceLY = (1 / TISSUE_ALPHA) * GammaY * LY;
}

// Loop through list of edge_candidates to store edges that point outward (cell of type 0 radially outward from the edge).
Edge* MegaTissue::edge_candidate_outward(const vector<Edge*> &edge_candidates){
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

//Called by HelperResolveVertices_ResolveFourfoldVertex to check the stability of one configuration. Splits vertices so edges (e1,e2) and (e3,e4) are paired and cells c2 and c4 are neighbors retruns false if the vertex is unstable in that topology
double MegaTissue::HelperResolveVertices_CheckStability(Vertex* vert, Edge* edge1, Edge* edge2, Edge* edge3, Edge* edge4, Cell* cell1, Cell* cell2, Cell* cell3, Cell* cell4, Vect* Force){
	//First make new objects which are coppies of the objects in the linked lists we will minipulate these objecs and they will be destroyed when the function ends
	MegaCell* mega_c1 = dynamic_cast<MegaCell*>(cell1);
	MegaCell* mega_c2 = dynamic_cast<MegaCell*>(cell2);
	MegaCell* mega_c3 = dynamic_cast<MegaCell*>(cell3);
	MegaCell* mega_c4 = dynamic_cast<MegaCell*>(cell4);
	MegaCell c1 = (*mega_c1);
	MegaCell c2 = (*mega_c2);
	MegaCell c3 = (*mega_c3);
	MegaCell c4 = (*mega_c4);
	MegaEdge* mega_edge1 = dynamic_cast<MegaEdge*>(edge1);
	MegaEdge* mega_edge2 = dynamic_cast<MegaEdge*>(edge2);
	MegaEdge* mega_edge3 = dynamic_cast<MegaEdge*>(edge3);
	MegaEdge* mega_edge4 = dynamic_cast<MegaEdge*>(edge4);
	MegaEdge e1 = (*mega_edge1);
	MegaEdge e2 = (*mega_edge2);
	MegaEdge e3 = (*mega_edge3);
	MegaEdge e4 = (*mega_edge4);
	int stability;
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

	//Calculate the tension on the new edge
	Vect ForceNewEdge = enew.GetTension_ForceOnVertex(&v12);

	//Check the stability
	stability = (*Force).Magnitude() - ForceNewEdge.Magnitude();
	return stability;
}