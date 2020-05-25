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
#include <typeinfo>
#include <cassert>
#include<numeric>

#include "Constants.h"
#include "FiberConstants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "FiberTissueClass.h"
#include "FiberCellClass.h"
#include "FiberEdgeClass.h"
#include "NucleationSite.h"


using namespace std;

// Reads the data and sets all of the global variables
void FiberTissue::HelperReadInputData_SetGlobalValues(ifstream* datafile){

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

	if (globals.count("A0") == 1){ A0 = atof(globals["A0"].c_str()); }
	else if (globals.count("A0") == 0){ cout << "WARNING: A0 is not given in input file, used default value " << A0 << endl; }
	else { cout << "ERROR: Global variable A0 has multiple definition in input file, used default value LY = " << A0 << endl; }

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

	if (globals.count("LAMBDA") == 1){ LAMBDA = atof(globals["LAMBDA"].c_str()); }
	else if (globals.count("LAMBDA") == 0){ cout << "WARNING: LAMBDA is not given in input file, used default value " << LAMBDA << endl; }
	else { cout << "ERROR: Global variable LAMBDA has multiple definition in input file, used default value K = " << LAMBDA << endl; }

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

	if (globals.count("TISSUE_ALPHA") == 1){ TISSUE_ALPHA = atof(globals["TISSUE_ALPHA"].c_str()); }
	else if (globals.count("TISSUE_ALPHA") == 0){ cout << "WARNING: TISSUE_ALPHA is not given in input file, used default value " << TISSUE_ALPHA << endl; }
	else { cout << "ERROR: Global variable TISSUE_ALPHA has multiple definition in input file, used default value TISSUE_ALPHA = " << TISSUE_ALPHA << endl; }

	if (globals.count("ALPHAM") == 1){ alphaM = atof(globals["ALPHAM"].c_str()); }
	else if (globals.count("ALPHAM") == 0){ cout << "WARNING: ALPHAM is not given in input file, used default value " << alphaM << endl; }
	else { cout << "ERROR: Global variable ALPHAM has multiple definition in input file, used default value ALPHAM = " << alphaM << endl; }

	if (globals.count("TISSUE_EXT_STRESS_X") == 1){ TISSUE_EXT_STRESS_X = atof(globals["TISSUE_EXT_STRESS_X"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_X") == 0){ cout << "WARNING: TISSUE_EXT_STRESS_X is not given in input file, used default value " << TISSUE_EXT_STRESS_X << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_X has multiple definition in input file, used default value TISSUE_EXT_STRESS_X = " << TISSUE_EXT_STRESS_X << endl; }

	if (globals.count("TISSUE_EXT_STRESS_Y") == 1){ TISSUE_EXT_STRESS_Y = atof(globals["TISSUE_EXT_STRESS_Y"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS_Y") == 0){ cout << "WARNING: TISSUE_EXT_STRESS_Y is not given in input file, used default value " << TISSUE_EXT_STRESS_Y << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS_Y has multiple definition in input file, used default value TISSUE_EXT_STRESS_Y = " << TISSUE_EXT_STRESS_Y << endl; }

	if (globals.count("A0SKEW") == 1){ A0SKEW = atof(globals["A0SKEW"].c_str()); }
	else if (globals.count("A0SKEW") == 0){ cout << "WARNING: A0SKEW is not given in input file, used default value " << A0SKEW << endl; }
	else { cout << "ERROR: Global variable A0SKEW has multiple definition in input file, used default value A0SKEW = " << A0SKEW << endl; }

	if (globals.count("TAU") == 1){ TAU = atof(globals["TAU"].c_str()); }
	else if (globals.count("TAU") == 0){ cout << "WARNING: TAU is not given in input file, used default value " << TAU << endl; }
	else { cout << "ERROR: Global variable TAU has multiple definition in input file, used default value TAU = " << TAU << endl; }

	if (globals.count("TISSUE_EXT_STRESS") == 1){ TISSUE_EXT_STRESS = atof(globals["TISSUE_EXT_STRESS"].c_str()); }
	else if (globals.count("TISSUE_EXT_STRESS") == 0){ cout << "WARNING: TISSUE_EXT_STRESS is not given in input file, used default value " << TISSUE_EXT_STRESS << endl; }
	else { cout << "ERROR: Global variable TISSUE_EXT_STRESS has multiple definition in input file, used default value TISSUE_EXT_STRESS = " << TISSUE_EXT_STRESS << endl; }

	if (globals.count("TAU_STRESS") == 1){ TAU_STRESS = atof(globals["TAU_STRESS"].c_str()); }
	else if (globals.count("TAU_STRESS") == 0){ cout << "WARNING: TAU_STRESS is not given in input file, used default value " << TAU_STRESS << endl; }
	else { cout << "ERROR: Global variable TAU_STRESS has multiple definition in input file, used default value TAU_STRESS = " << TAU_STRESS << endl; }

	if (globals.count("SPRINGCONST") == 1){ SPRINGCONST = atof(globals["SPRINGCONST"].c_str()); }
	else if (globals.count("SPRINGCONST") == 0){ cout << "WARNING: SPRINGCONST is not given in input file, used default value " << SPRINGCONST << endl; }
	else { cout << "ERROR: Global variable SPRINGCONST has multiple definition in input file, used default value SPRINGCONST = " << SPRINGCONST << endl; }

	if (globals.count("DTENSIONCONST") == 1){ DTENSIONCONST = atof(globals["DTENSIONCONST"].c_str()); }
	else if (globals.count("DTENSIONCONST") == 0){ cout << "WARNING: DTENSIONCONST is not given in input file, used default value " << DTENSIONCONST << endl; }
	else { cout << "ERROR: Global variable DTENSIONCONST has multiple definition in input file, used default value DTENSIONCONST = " << DTENSIONCONST << endl; }

	if (globals.count("ANGLECONST") == 1){ ANGLECONST = atof(globals["ANGLECONST"].c_str()); }
	else if (globals.count("ANGLECONST") == 0){ cout << "WARNING: ANGLECONST is not given in input file, used default value " << ANGLECONST << endl; }
	else { cout << "ERROR: Global variable ANGLECONST has multiple definition in input file, used default value ANGLECONST = " << ANGLECONST << endl; }

	if (globals.count("BETAZ") == 1){ BETAZ = atof(globals["BETAZ"].c_str()); }
	else if (globals.count("BETAZ") == 0){ cout << "WARNING: BETAZ is not given in input file, used default value " << BETAZ << endl; }
	else { cout << "ERROR: Global variable BETAZ has multiple definition in input file, used default value BETAZ = " << BETAZ << endl; }

	if (globals.count("BETAGAMMA") == 1){ BETAGAMMA = atof(globals["BETAGAMMA"].c_str()); }
	else if (globals.count("BETAGAMMA") == 0){ cout << "WARNING: BETAGAMMA is not given in input file, used default value " << BETAGAMMA << endl; }
	else { cout << "ERROR: Global variable BETAGAMMA has multiple definition in input file, used default value BETAGAMMA = " << BETAGAMMA << endl; }

	if (globals.count("B") == 1){ B = atof(globals["B"].c_str()); }
	else if (globals.count("B") == 0){ cout << "WARNING: B is not given in input file, used default value " << B << endl; }
	else { cout << "ERROR: Global variable B has multiple definition in input file, used default value B = " << B << endl; }

	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//
}

//! Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids. Starts at current position of *datafile and reads until keyword 'End'.
void FiberTissue::HelperReadInputData_CreateObjects(ifstream* datafile){
	string line, key, value, edg;			//temporary variables to store lines and words 
	string type;							//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline, ssedges;			//ssline is a string stream so we can comvert strings to stringstreams
	vector<string> keywords;				//vector of strings to hold the keywords
	vector<int> listedges;
	int i;
	map <string, string> objectInfo;		//map to hold the keywords and their data
	int newid, v1, v2, enumb, typenumb;					//temp varables to store the object info
	double x, y, thisA0;
	bool xf, yf;

	cout << "I'm making a fiber tissue!" << endl;

	//VERTICES read the file untill you get to the word vertices
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
	while (1>0){
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
		(*vptr).SavedInfoForEvolve_MechForce.x = NAN;
		(*vptr).SavedInfoForEvolve_MechForce.y = NAN;
		(vertices).push_back(vptr);
		vertexIDs.push_back(newid);								//and add its id to the list of ids
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*vertices.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}

	//EDGES then get the keywords on the next line
	keywords.clear();
	getline(*datafile, line);
	ssline.clear();
	ssline << line;
	while (getline(ssline, key, '\t')){
		keywords.push_back(key);
	}
	//now recorde the data and make the edge object untill you get to the word cells signifying the end of the data
	while (1>0){
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
		Edge* eptr = new FiberEdge;
		(*eptr).id = newid;
		(*eptr).IOVertexID[0] = v1;
		(*eptr).IOVertexID[1] = v2;
		(*eptr).xflag = xf;
		(*eptr).yflag = yf;
		(edges).push_back(eptr);		// Create a vertes with an id and x,y pos and put it in the master list of all veritces
		edgeIDs.push_back(newid);								//and add its id to the list of ids
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*edges.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}
	cout << "starting to read cells" << endl;
	//CELLS then get the keywords on the next line
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
		if (line.compare("Fibers") == 0){ break; }				//break if keyworkd 'End' is found
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
			thisA0 = atof(objectInfo["A0"].c_str());
		}
		ssedges.clear();
		ssedges << objectInfo["EdgeList"];
		//now worry about the edge list
		listedges.clear();
		while (getline(ssedges, edg, ',')){
			listedges.push_back(atoi(edg.c_str()));
		}
		//make the cell object
		Cell* cptr = new FiberCell;
		(*cptr).id = newid;
		(*cptr).type = typenumb;
		(*cptr).A0 = thisA0;
		(*cptr).IOEdgeID = listedges;
		cells.push_back(cptr);
		cellIDs.push_back(newid);
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}
	cout << "I FINISHED cells" << endl;

	cout << "Starting to read fibs" << endl;
	int thiscellid;
	//FIBERS then get the keywords on the next line
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
		if (objectInfo.count("cellid")){
			thiscellid = atof(objectInfo["cellid"].c_str());
		}
		ssedges.clear();
		
		// if fiber belongs to cell with given cell id, increment the fiber counter for that given cell.
		for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){
			if ((*(*cit)).id == thiscellid){
				FiberCell* cptr = dynamic_cast<FiberCell*>((*cit)); 
				(*cptr).number_of_fibers_from_data++;
				break;
			}
		}

		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//
	}
	cout << "I finished fibs" << endl;
}

//Create a spring cell
Cell* FiberTissue::CreateCell(int typ, vector<Edge*> edgs){
	int largestID, newid;
	list<Cell*>::iterator cit;
	//Get the id
	largestID = (*max_element(cellIDs.begin(), cellIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;								//Now make a new ID that is unique

	if (newid > 1000000000) { cout << "ERROR to many cells" << endl; }
	//create the cell and add the cell to the lists
	Cell* cptr = new FiberCell;
	(*cptr).id = newid;
	(*cptr).type = typ;
	(*cptr).edges = edgs;
	(cells).push_back(cptr);
	cellIDs.push_back(newid);

	for (cit = cells.begin(); cit != cells.end(); cit++){
		if ((*(*cit)).id == newid){
			return &(*(*cit));
			break;
		}
	}
}

//Create a fiber vertex
NucleationSite* FiberTissue::CreateNucleationSite(Vertex* parentv, FiberEdge* e1, FiberEdge* e2, FiberCell* c){
	int largestID, newid;
	
	//Get the id
	largestID = (*max_element(NucleationSiteIDs.begin(), NucleationSiteIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;								//Now make a new ID that is unique

	if (newid > 1000000000) { cout << "ERROR to many fiber verticies" << endl; }

	NucleationSite* vptr = new NucleationSite(newid,parentv,e1,e2,c);
	(fibernucsites).push_back(vptr);
	NucleationSiteIDs.push_back(newid);

	return vptr;
}

//! Takes a tissue and creates fiber verticies from the parent verticies and sets their pointers 
void FiberTissue::CreateNucleationSites(){
	vector<Vertex*> cwvertices;
	vector<Cell*> cwcells;
	vector<Edge*> cwedges;
	FiberEdge *a, *b; 
	FiberCell* c; 
	NucleationSite* fv; 
	int number_of_verts; 
	int N;

	//For every parent vertex
	for (list<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++){
		Vertex* v = &(*(*it));
		N = (*v).edges.size(); 
		//For each pair of edges around the vertex 
		if (N > 2){  // if the vertex only has one edge do nothing
			
			//First get the cw edges and cells (if you walked in a CW circle around v starting at e0 you would see: e0,c0,e1,c1,e2,c2,e3,c3,e0)
			cwcells = (*v).GetCWCells();		//Get all the cells 
			cwedges = (*v).GetCWEdges();	//and the edges 
			
			for (int i = 0; i < N; i++){
				a = dynamic_cast<FiberEdge*>(cwedges[i]);								//get the pair 
				b = dynamic_cast<FiberEdge*>(cwedges[(i + 1) % N]);
				c = dynamic_cast<FiberCell*>(cwcells[i]);	//and the cell between them 
				// if the cell is NOT the boundary cell then make a NucleationSite between the two edges 
				if (!(*c).isBoundary){
					fv=CreateNucleationSite(v, a, b, c); //makes a fiber vertex and adds it to the list of fiber verticies 
					(*c).fibernucsites.push_back(fv); 
				} //end pair
			} //end vertex
		}
	}//end vertex list 

	////// CODE FROM HAYDEN: IF YOU WANT TO PUT A NUCLEATION SITE AT EACH VERTEX (INSTEAD OF JUST VERTICES THAT ARE NOT DANGLING)
	////For every cell, create nucleation sites
	//for (list<Cell *>::iterator it = cells.begin(); it != cells.end(); it++){
	//	c = dynamic_cast<FiberCell*>(&(*(*it)));
	//	if ((*c).isCompleete){
	//		cwvertices = (*c).HelperOrderVertices();
	//		number_of_verts = cwvertices.size();
	//		for (int i = 0; i < number_of_verts; i++){
	//			Vertex* this_vert = cwvertices[i];
	//			Vertex* neighbor_vert_1 = cwvertices[(i + 1) % number_of_verts];
	//			Vertex* neighbor_vert_2;
	//			if (i - 1 == -1){
	//				neighbor_vert_2 = cwvertices[number_of_verts - 1];
	//			}
	//			else{
	//				neighbor_vert_2 = cwvertices[i - 1];
	//			}
	//			Edge* edge1_ptr = (*this_vert).GetEdgeBetween(neighbor_vert_1);
	//			Edge* edge2_ptr = (*this_vert).GetEdgeBetween(neighbor_vert_2);
	//			a = dynamic_cast<FiberEdge*>(edge1_ptr);								//get the pair 
	//			b = dynamic_cast<FiberEdge*>(edge2_ptr);
	//			// if the cell is NOT the boundary cell then make a NucleationSite between the two edges 
	//			if (!(*c).isBoundary){
	//				fv = CreateNucleationSite(this_vert, a, b, c); //makes a fiber vertex and adds it to the list of fiber verticies 
	//				(*c).fibernucsites.push_back(fv);
	//			} //end pair
	//		} //end vertex
	//		if (cwvertices.size() != (*c).fibernucsites.size()){ 
	//			cout << "something is terribly wrong" << endl; 
	//		}
	//	}
	//}//end vertex list
}

//Creates a fiber and adds it to the tissue
Fiber* FiberTissue::CreateFiber(NucleationSite* v){
	int largestID, newid, i;
	FiberCell* c = (*v).cell; 
	FiberEdge* e1 = (*v).edge1; 
	FiberEdge* e2 = (*v).edge2;
	double vx, vy, direction, xloc; 
	double u = starting_distance_from_first_vert; //um to put fiber away from nucleation site

	//Get the new id for the fiber
	largestID = (*max_element(fiberIDs.begin(), fiberIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;														//Now make a new ID that is unique
	if (newid > 1000000000) { cout << "ERROR to many fibers" << endl; }

	//create the fiber and add the fiber to the lists
	Fiber* fptr = new Fiber;
	(*fptr).id = newid;
	(fibers).push_back(fptr);
	fiberIDs.push_back(newid);

	//Record the peramiters
	(*fptr).birthtime = CURRENTTIME;
	(*fptr).passedvertex = false;

	//Add the x placment and a pointer to its cell
	(*fptr).cell = &(*c);
	(*fptr).v1 = (*v).parentvertex;
	(*fptr).v2 = (*v).parentvertex;

	//Add the fiber to its cell
	(*c).fibers.push_back(fptr);

	//Add the edges to the fiber and the fiber to its edges
	(*fptr).edge1 = e1;
	(*fptr).edge2 = e2; 
	(*(*fptr).edge1).fibers.push_back(fptr);
	(*(*fptr).edge2).fibers.push_back(fptr);

	//Get the endpoint positions
	(*fptr).dedge1_x = u;
	(*fptr).dedge2_x = u;
	//Calculate distance to next vertex along each edge
	double dist_edge1_fib = (*fptr).DistanceToNextVertexonEdge1();
	double dist_edge2_fib = (*fptr).DistanceToNextVertexonEdge2();
	while (dist_edge1_fib <= distance_to_next_vert_cutoff) {
		cout << "something is horribly wrong" << endl;
	}
	while (dist_edge2_fib <= distance_to_next_vert_cutoff) { 
		cout << "something is horribly wrong" << endl;
	}
	
	return fptr;
}

////Removes a fiber from the tissue
void FiberTissue::DestroyFiber(Fiber* f){
	//remove the fiber from the tissue
	for (list<Fiber*>::iterator fit = fibers.begin(); fit != fibers.end(); fit++){			//loop through all of the edges in the tissue
		if ((*(*fit)).id == (*f).id){										//if the curent vertex is one of the vertices to be removed
			delete (*fit);
			fibers.erase(fit);												//remove it
			break;
		}
	}
}

//writes the data including fiber info
void FiberTissue::WriteOutputData(ofstream* datafile){
	time_t rawtime;					//to deal with writing the date out
	struct tm * timeinfo;
	double rx, ry, n, a, ux, uy; 

	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
		FiberCell* fc = dynamic_cast<FiberCell*>(*cit);
		if ((*fc).isCompleete){
			(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
		}
	}

	//First write the info about what you just ran and the time you ran it
	time(&rawtime);
	//	timeinfo = localtime(&rawtime);
	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;
	*datafile << "MYOSINDYNAMICS\t" << MYOSINDYNAMICS << endl;
	//	*datafile << "Date:\t" << asctime(timeinfo) << endl;

	//Now write the global variables you used
	*datafile << "Global:" << endl;
	*datafile << "CURRENTTIME\t" << CURRENTTIME << endl;
	*datafile << "SEED\t" << SEED << endl;
	*datafile << "TIMEINTERVAL\t" << TIMEINTERVAL << endl;
	*datafile << "LX\t" << LX << endl;
	*datafile << "LY\t" << LY << endl;
	*datafile << "TOTALTIME\t" << TOTALTIME << endl;
	*datafile << "INITIALSTEPSIZE\t" << INITIALSTEPSIZE << endl;
	*datafile << "endGlobal" << endl;

	//print the vertex info
	*datafile << "Vertices" << endl;
	*datafile << "id\tx\ty\tn\ta" << endl;
	for (list<NucleationSite*>::iterator vit = fibernucsites.begin(); vit != fibernucsites.end(); vit++){			//loop through all of the vertices in the tissue
		NucleationSite* v = &(*(*vit)); 
		rx = (*v).GetLocation().x; 
		ry = (*v).GetLocation().y;
		n = (*v).GetOutwardAngle(); 
		a = (*v).GetAngle();
		*datafile << (*v).id << "\t" << rx << "\t" << ry << "\t" << n << "\t" << a << endl;						//write it's info to the file
	}

	//print the edge info
	*datafile << "Edges" << endl;
	*datafile << "id\tv1x\tv1y\tv2x\tv2y" << endl;
	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
		FiberEdge* e = dynamic_cast<FiberEdge*>((*eit));
		Vertex* v1 = (*e).vertex1;
		Vertex* v2 = (*e).vertex2;
		rx = (*v1).r.x; 
		ry = (*v1).r.y;
		ux = (*v2).r.x;
		uy = (*v2).r.y;
		*datafile << (*(*eit)).id << "\t" << rx << "\t" << ry << "\t" << ux << "\t" << uy << endl;			//write its info to the file
	}

	//print the cell info
	int i,j;
	vector<vector<double>> xypos; 
	vector<double> xpos;
	vector<double> ypos;

	*datafile << "Cells" << endl;
	*datafile << "id\tsize\tarea\txextent\tyextent\txdisttodeathvert\tNneuc\tNdeath\tNDatafibs\tNfib\txpos\typos\tcos2nvect\topeninganglevect\txdisttodeathvertvect" << endl;
	for (list<Cell*>::iterator it = cells.begin(); it != cells.end(); it++){			//loop through all of the edges in the tissue
		FiberCell* c = dynamic_cast<FiberCell*>((*it));
		if ((*c).isCompleete){
			xypos = (*c).HelperOrderVerticesPos();
			xpos = xypos[0];
			ypos = xypos[1];
			// vector for storing whether a specified vertex within a given cell is a nucleating site (with proper orientation)
			vector<double> cos2n_vect;
			vector<double> opening_angle_vect;
			vector<double> x_dist_to_death_vert_vect;
			cos2n_vect.resize((*c).fibernucsites.size());
			opening_angle_vect.resize((*c).fibernucsites.size());
			x_dist_to_death_vert_vect.resize((*c).fibernucsites.size());

			//only print out cells that are not dangling (check that each vertex within each cell has been checked for proper orientation).
			if (xpos.size() == (*c).fibernucsites.size()){
				for (int i = 0; i < cos2n_vect.size(); i++){

					// Store the pointer to nucleation site (vertex within cell).
					NucleationSite* this_nuc_site_ptr = (*c).fibernucsites.at(i);

					// Store the parent vertex.
					Vertex* this_parent_vertex_ptr = (*this_nuc_site_ptr).parentvertex;

					//Store position of this parent vertex.
					double parent_vertex_xpos = (*this_parent_vertex_ptr).r.x;
					double parent_vertex_ypos = (*this_parent_vertex_ptr).r.y;

					// Initialize int for storing which vertex (in xpos and ypos) corresponds to this nucleating site
					int corresponding_vertex_index = -1;

					// Check if has proper orientation.
					int is_nuc_site_test = ((*this_nuc_site_ptr).ProbFiberBirth() > 0.0);
					// set cos2n
					double cos2n_calc = (*this_nuc_site_ptr).Calculate_Cos2OutwardAngle();

					double opening_angle_calc = (*this_nuc_site_ptr).GetAngle();

					double x_dist_to_death_vert = (*this_nuc_site_ptr).Calc_XDist_toDeathVertex();

					// Loop through xpos, ypos to check which vertex this nuc site corresponds to.
					for (int j = 0; j < cos2n_vect.size(); j++){
						if (parent_vertex_xpos == xpos[j] && parent_vertex_ypos == ypos[j]){ corresponding_vertex_index = j; }
					}
					
					// Store vertex orientation.
					cos2n_vect[corresponding_vertex_index] = cos2n_calc;
					opening_angle_vect[corresponding_vertex_index] = opening_angle_calc;
					x_dist_to_death_vert_vect[corresponding_vertex_index] = x_dist_to_death_vert;
				}

				*datafile << (*c).id << "\t" << (*c).sizeclass << "\t" << (*c).GetArea() << "\t" << (*c).GetMaxHorizontalLength() << "\t" << (*c).GetMaxVerticalLength() << "\t" << (*c).CalcAverageXDist_toDeathVertex() << "\t" << (*c).CountNucleatingSites() << "\t" << (*c).CountFiberDeathSites() << "\t" << (*c).number_of_fibers_from_data << "\t" << (*c).fibers.size() << "\t";
				for (i = 0; i < xpos.size(); i++){
					*datafile << xpos[i] << ",";
				}
				*datafile << "\t";
				for (i = 0; i < ypos.size(); i++){
					*datafile << ypos[i] << ",";
				}
				*datafile << "\t";
				for (i = 0; i < (*c).fibernucsites.size(); i++){
					*datafile << cos2n_vect[i] << ",";
				}
				*datafile << "\t";
				for (i = 0; i < (*c).fibernucsites.size(); i++){
					*datafile << opening_angle_vect[i] << ",";
				}
				*datafile << "\t";
				for (i = 0; i < (*c).fibernucsites.size(); i++){
					*datafile << x_dist_to_death_vert_vect[i] << ",";
				}
				*datafile << endl;
			}
		}
	}
	
	//print fiber info
	*datafile << "Fibers" << endl;
	*datafile << "id\tr1x\tr1y\tr2x\tr2y" << endl;
	for (list<Fiber*>::iterator fit = fibers.begin(); fit != fibers.end(); fit++){			//loop through all of the edges in the tissue
		*datafile << (*(*fit)).id << "\t" << (*(*fit)).r1().x << "\t" << (*(*fit)).r1().y << "\t" << (*(*fit)).r2().x << "\t" << (*(*fit)).r2().y << endl;			//write its info to the file
	}

	*datafile << "End" << endl << endl;
}

//writes the data including fiber info
void FiberTissue::WriteOutputStatistics(ofstream* datafile){
	time_t rawtime;					//to deal with writing the date out
	struct tm * timeinfo;
	double rx, ry, n, a, ux, uy;
	fiberstats fs; 

	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
		FiberCell* fc = dynamic_cast<FiberCell*>(*cit);
		if ((*fc).isCompleete){
			(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
		}
	}

	//First write the info about what you just ran and the time you ran it
	time(&rawtime);
	//	timeinfo = localtime(&rawtime);
	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;
	*datafile << "MYOSINDYNAMICS\t" << MYOSINDYNAMICS << endl;
	//	*datafile << "Date:\t" << asctime(timeinfo) << endl;

	//Now write the global variables you used
	*datafile << "Global:" << endl;
	*datafile << "CURRENTTIME\t" << CURRENTTIME << endl;
	*datafile << "SEED\t" << SEED << endl;
	*datafile << "TIMEINTERVAL\t" << TIMEINTERVAL << endl;
	*datafile << "LX\t" << LX << endl;
	*datafile << "LY\t" << LY << endl;
	*datafile << "TOTALTIME\t" << TOTALTIME << endl;
	*datafile << "INITIALSTEPSIZE\t" << INITIALSTEPSIZE << endl;
	*datafile << "INITIALVEL\t" << INITIALVEL << endl;           
    *datafile << "SLOWVEL" << SLOWVEL<< endl;            
	*datafile << "NUCLEATIONRATE" << NUCLEATIONRATE << endl; 		//Nucleation rate per neucleating site in min-1
	*datafile << "INITIALDEATHRATE" << INITIALDEATHRATE << endl;   // deathrate for fibers that have not passed a vertex in min-1
	*datafile << "PASSEDDEATHRATE" << PASSEDDEATHRATE << endl;     // death rate for fibers that have passed a vertex in min-1
	*datafile << "endGlobal" << endl;

	*datafile << "FiberData" << endl;

	*datafile << "id\t" << "cellid\t" << "cellA\t" << "cellS\t" << "lifetime\t" << "passedvertex" << endl;
	for (int i=0; i < fiberstatistics.size(); i++){
		fs = fiberstatistics[i]; 
		*datafile << fs.id << " \t" << fs.cellid << "\t" << fs.cellarea << "\t" << fs.cellsize << "\t" << fs.lifetime << "\t" << fs.passedvertex << endl;
	}

	*datafile << "End" << endl << endl;
}

//Initializes the tissue after it is read in (removes dangeling edges and cells ...)
void FiberTissue::InitializeFromData(){
	vector<Vertex*> delverts;
	vector<Cell*> delcells;
	vector<Edge*> deledges;
	vector<Edge*> deledges2;
	vector<Edge*> nulledges;

	cout << "Starting to initialize" << endl;
	//Set the initial forces on vertices to 0
	Vect o(0.0, 0.0);
	for (list<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++){
		Vertex* v = &(*(*it));
		(*v).SavedInfoForEvolve_MechForce = o;
	}
	cout << "Initial forces set" << endl;
	//Remove any edges that cross the boundary
	for (list<Edge*>::iterator it = edges.begin(); it != edges.end(); it++){
		Edge* e = &(*(*it));
		if ((*e).yflag || (*e).xflag){
			//cout << "eid " << (*e).id << endl;
			deledges2.push_back(e);
			Vertex* v1 = (*e).vertex1;
			Vertex* v2 = (*e).vertex2;
			Cell* c1 = (*e).cell1;
			Cell *c2 = (*e).cell2;
			(*v1).RemoveEdge(e);
			(*v2).RemoveEdge(e);
			(*c1).RemoveEdge(e);
			(*c2).RemoveEdge(e);
		}
	}
	for (int i = 0; i < deledges2.size(); i++){
		//cout << "destroying edge " << (*deledges2[i]).id << endl;
		(*this).DestroyEdge(deledges2[i]);
	}
	cout << "edges removed" << endl;

	//Set the isCompleete flag
	for (list<Cell*>::iterator it = cells.begin(); it != cells.end(); it++){
		Cell* c = &(*(*it));
		FiberCell* sc = dynamic_cast<FiberCell*>(c);
		(*sc).SetIsCompleeteFlag();
	}
	cout << "IsCompleeteFlag set" << endl;
	//Get all the verticies in CW order 
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
		Cell* c = &(*(*cit));
		FiberCell* sc = dynamic_cast<FiberCell*>(c);
		if ((*sc).isCompleete){															//if the cells are compleete
			(*(*cit)).SavedInfoForEvolve_CWVerts = (*(*cit)).GetClockWiseVertices();	//and get the vertices in the correct clockwise order 
			(*sc).type = 2;
		}
	}

	//Remove noncompleete cells
	FiberCell* outercell; //Create a spring cell
	outercell = dynamic_cast<FiberCell*>(CreateCell(1, nulledges));
	(*outercell).isCompleete = false;
	(*outercell).isBoundary = true;

	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++){		//Loop through the cells 
		Cell* c = &(*(*cit));
		FiberCell* sc = dynamic_cast<FiberCell*>(c);
		if (!(*sc).isCompleete && (*sc).id != (*outercell).id){															//if the cells are  not compleete
			delcells.push_back(sc); //mark them for removal 
			//for every edge in the cell 
			for (int i = 0; i < (*sc).edges.size(); i++){
				(*(*sc).edges[i]).UpdateCell(sc, outercell);    //replace the pointer to theis cell with one to the outer cell
				(*outercell).edges.push_back(&(*(*sc).edges[i])); //add the edge to the outercell 
			}
		}
	}

	for (int i = 0; i < delcells.size(); i++){
		DestroyCell(delcells[i]);
	}
	cout << "destroyed non-compleete cells" << endl;

	//Remove dangeling edges and vertices
	Vertex* keepv = 0;
	for (list<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++){		//for all vertices
		if ((*(*it)).edges.size() == 1){												//if the vetres has only one edge
			delverts.push_back(&(*(*it))); //put the vertex in the delete pile
			deledges.push_back(&(*(*(*it)).edges[0])); //and its single edge
			//get the other vertex wich is atached to the main body of cells
			keepv = (*(*(*it)).edges[0]).GetOtherVertex((*it));
			(*keepv).RemoveEdge(&(*(*(*it)).edges[0])); //remove the edge from it
		}
	}
	for (int i = 0; i < deledges.size(); i++){
		DestroyEdge(deledges[i]);
	}
	for (int i = 0; i < delverts.size(); i++){
		DestroyVertex(delverts[i]);
	}

	//Remove dangeling edges and vertices again
	keepv = 0;
	for (list<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++){		//for all vertices
		if ((*(*it)).edges.size() == 1){												//if the vetres has only one edge
			delverts.push_back(&(*(*it))); //put the vertex in the delete pile
			deledges.push_back(&(*(*(*it)).edges[0])); //and its single edge
			//get the other vertex wich is atached to the main body of cells
			keepv = (*(*(*it)).edges[0]).GetOtherVertex((*it));
			(*keepv).RemoveEdge(&(*(*(*it)).edges[0])); //remove the edge from it
		}
	}
	for (int i = 0; i < deledges.size(); i++){
		DestroyEdge(deledges[i]);
	}
	for (int i = 0; i < delverts.size(); i++){
		DestroyVertex(delverts[i]);
	}
	cout << "removed dangeling edges" << endl << endl;

	//Set length scale (initaly 1 unit length = 0.1 mu) 
	for (list<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++){		//for all vertices
		(*(*it)).r.x = (*(*it)).r.x / 10; 
		(*(*it)).r.y = (*(*it)).r.y / 10;
	}
	cout << "Set lengthscale: 1 unit length = 1  micrometre" << endl; 
}

void FiberTissue::RecordCellSize(){
	vector<double> cellareas; 
	double ncells = 0; 
	double meansize; 

	for (list<Cell*>::iterator it = cells.begin(); it != cells.end(); it++){
		FiberCell* c = dynamic_cast<FiberCell*>(&(*(*it)));
		if ((*c).isCompleete){
			cellareas.push_back( (*c).GetArea() );
			ncells++; 
		}
	}

	size_t n = cellareas.size() / 2;
	nth_element(cellareas.begin(), cellareas.begin() + n, cellareas.end());
	meansize=cellareas[n];

	for (list<Cell*>::iterator it = cells.begin(); it != cells.end(); it++){
		FiberCell* c = dynamic_cast<FiberCell*>(&(*(*it)));
		if ((*c).isCompleete){
			if ((*c).GetArea() < meansize)
				(*c).sizeclass = 1;
			else
				(*c).sizeclass = 2;
		}
		else
			(*c).sizeclass = 0; 
	}


}

//Evolves the fibers by giving them some velocity 
void FiberTissue::EvolveFibers(double dt, double tstep){
	double prob, p; 
	double tend; 
	double vel; 
	double dx; 
	FiberEdge *e1, *e2;
	FiberEdge *en = 0;
	FiberEdge *etemp = 0;
	Edge* edge=0; 
	FiberCell* c; 
	Vertex *v, *vn;
	vector<Fiber*> brokenfibers; 

	//Step 0: initialaze 
	tend = CURRENTTIME + tstep; 

	while (CURRENTTIME < tend){
		//Step 1: Birth new fibers 
		for (list<NucleationSite*>::iterator vit = fibernucsites.begin(); vit != fibernucsites.end(); vit++){
			NucleationSite* v = &(*(*vit));
			prob = (*v).ProbFiberBirth()*dt;
			p = distribution(generator);
			if (p < prob) {
				CreateFiber(v);
			}
		}

		//Step 2: Move all of the fibers 
		for (list<Fiber*>::iterator fit = fibers.begin(); fit != fibers.end(); fit++){
			Fiber* f = &(*(*fit));
			vel = (*f).GetVelocity(); 
			dx = dt*vel; 
			(*f).dedge1_x += dx; 
			(*f).dedge2_x += dx;
		}

		//Step 3: Check for any edge changes and make them 
		for (list<Fiber*>::iterator fit = fibers.begin(); fit != fibers.end(); fit++){
			Fiber* f = &(*(*fit));
			e1 = (*f).edge1;
			e2 = (*f).edge2;
			c = (*f).cell; 
			//Calculate distance to next vertex along each edge
			double dist_edge1_fib = (*f).DistanceToNextVertexonEdge1();
			double dist_edge2_fib = (*f).DistanceToNextVertexonEdge2();

			//if the fiber has walked allong edge 1 
			if (dist_edge1_fib <= distance_to_next_vert_cutoff){
				//if this is not the last edge 
				if (!(*f).IsLastEdge()){
					//Get the next edge 
					v = (*f).v1;
					vn = (*e1).GetOtherVertex(v);
					edge = (*vn).GetOtherEdgeInCell(c, e1);
					en = dynamic_cast<FiberEdge*>(edge);
					//swap the edges and the verts 
					(*f).v1 = vn;
					(*f).edge1 = en;
					(*e1).RemoveFiber(f);
					(*en).fibers.push_back(f);
					//reset the distance along the edge 
					(*f).dedge1_x = starting_distance_from_first_vert;
					(*f).passedvertex = true;
				}
				//Otherwise if it is the last edge, push it back.
				else{
					vel = (*f).GetVelocity();
					dx = dt*vel;
					(*f).dedge1_x -= dx; 
				}
			}
			//if the fiber has walked allong edge 2 
			if (dist_edge2_fib <= distance_to_next_vert_cutoff){
				//if this is not the last edge 
				if (!(*f).IsLastEdge()){
					//Get the next edge 
					v = (*f).v2;
					vn = (*e2).GetOtherVertex(v);
					edge = (*vn).GetOtherEdgeInCell(c, e2);
					en = dynamic_cast<FiberEdge*>(edge);
					//swap the edges and the verts 
					(*f).v2 = vn;
					(*f).edge2 = en;
					(*e2).RemoveFiber(f);
					(*en).fibers.push_back(f);
					(*f).dedge2_x = starting_distance_from_first_vert;
					(*f).passedvertex = true;
				}
				//Otherwise if this is the last edge, push it back. 
				else{
					vel = (*f).GetVelocity();
					dx = dt*vel;
					(*f).dedge2_x -= dx;
				}
			}
		}


		//Step 4: Destroy appropreate fibers 
		brokenfibers.clear(); 
		for (list<Fiber*>::iterator fit = fibers.begin(); fit != fibers.end(); fit++){
			//make a list of fibers to destroy 
			Fiber* f = &(*(*fit));
			prob = (*f).GetBreakProb()*dt;
			p = distribution(generator);
			if (p < prob) { brokenfibers.push_back(f); }
		}
		for (int i = 0; i < brokenfibers.size(); i++){
			//remove the fibers listed in broken fibers 
			Fiber* f = &(*brokenfibers[i]);

			//record its statistics 
			fiberstats fstats;
			fstats.id = (*f).id;
			fstats.cellid = (*(*f).cell).id;
			fstats.cellarea = (*(*f).cell).GetArea(); 
			fstats.cellsize = (*(*f).cell).sizeclass;
			fstats.lifetime = CURRENTTIME - (*f).birthtime; 
			fstats.passedvertex = (*f).passedvertex;
			fiberstatistics.push_back(fstats); 
			
			//then remove it 
			e1 = (*f).edge1;
			e2 = (*f).edge2;
			c = (*f).cell;
			(*e1).RemoveFiber(f);
			(*e2).RemoveFiber(f);
			(*c).RemoveFiber(f);
			DestroyFiber(f);
		}
		


		//Step 5: update time
		CURRENTTIME += dt; 
	}
}

