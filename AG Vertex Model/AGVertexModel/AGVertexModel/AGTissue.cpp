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


#include "Constants.h"
#include "AGConstants.h"
#include "TissueClass.h"
#include "AGTissue.h"
#include "AGCell.h"
#include "CellClass.h"
#include "VertexClass.h"

/*
* void HelperReadInputData_SetGlobalValues(ifstream* datafile)
*
* Reads the data and sets all of the global variables.
* Starts at current position of *datafile and look for keyword 'Global', reads until keyword 'endGlobal'.
*/
void AGTissue::HelperReadInputData_SetGlobalValues(ifstream* datafile){

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

	if (globals.count("SIMTIMERES") == 1){ SIMTIMERES = atof(globals["SIMTIMERES"].c_str()); }
	else if (globals.count("SIMTIMERES") == 0){ cout << "WARNING: SIMTIMERES is not given in input file, used default value " << SIMTIMERES << endl; }
	else { cout << "ERROR: Global variable SIMTIMERES has multiple definition in input file, used default value SIMTIMERES = " << SIMTIMERES << endl; }


	//******ADD MORE GLOBAL VARIABLES HERE ALSO ADD THEM TO THE CONSTANTS.H FILE AND GIVE THEM A DEFAULT VALUE IN MAIN*********************//
}

 
 /*
 * void HelperReadInputData_CreateObjects
 *
 * Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's 
 * list of ids.
 * Starts at current position of *datafile and reads until keyword 'End'.
 */
void AGTissue::HelperReadInputData_CreateObjects(ifstream* datafile){

	string line, key, value, edg;			//temporary variables to store lines and words 
	string type;							//what type of file is it (in case latter versions need to do slightly different things based on type)
	stringstream ssline, ssedges;			//ssline is a string stream so we can comvert strings to stringstreams
	vector<string> keywords;				//vector of strings to hold the keywords
	vector<int> listedges;
	int i,j;
	map <string, string> objectInfo;		//map to hold the keywords and their data
	int newid, v1, v2, enumb, typenumb;					//temp varables to store the object info
	double x, y;
	bool xf, yf;
	double A0;
	string potential;
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
		Edge* eptr = new Edge;
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

		//make the cell object
		AGCell* cptr = new AGCell;
		(*cptr).id = newid;
		(*cptr).type = typenumb;
		(*cptr).A0 = A0;
		(*cptr).IOEdgeID = listedges;
		(*cptr).DIVISIONPOTENTIAL = potential;

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
		//******************HOW TO ADD A NEW DATA MEMBER ***********************//
		//NEWDATA = atof( objectInfo["NEWDATA"].c_str() );
		//(*cells.back()).NEWDATANAME = NEWDATA
		//*************************************************************************//

		// this is a bit out of place, but I'll also set the initial box dimensions here
		LX0 = LX;
		LY0 = LY;
	}
}

/*
* void HelperWriteOutputData_Write2File
*
* Writes data into file specified by 'datafile'. If flag 'includeExtraInfo' is true, then additional information (such as tension, cell area, etc) is included.
*The extra information is non-essential in reconstructing the simulation state, and is just provided for convenience (e.g. for display in MATLAB).
* (should be helpful for incremental coding erro debug sessions)
*
*/

void AGTissue::HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo){
	time_t rawtime;					//to deal with writing the date out
	struct tm * timeinfo;

	*datafile << "INTEGRATORNAME\t" << INTEGRATORNAME << endl;

	//Now write the global variables you used
	*datafile << "CURRENTTIME\t" << CURRENTTIME << endl;
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
	*datafile << "INITIALSTEPSIZE\t" << INITIALSTEPSIZE << endl;
	*datafile << "error_tolerance_absolute_mechanical\t" << error_tolerance_absolute_mechanical << endl;
	*datafile << "TISSUE_ALPHA\t" << TISSUE_ALPHA << endl;
	*datafile << "TISSUE_EXT_STRESS_X\t" << TISSUE_EXT_STRESS_X << endl;
	*datafile << "TISSUE_EXT_STRESS_Y\t" << TISSUE_EXT_STRESS_Y << endl;
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
		*datafile << "id\tvertex1\tvertex2\txflag\tyflag\ttension\t_cell1\t_cell2" << endl;
		for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++){			//loop through all of the edges in the tissue
			Vertex* v1 = (*(*eit)).vertex1;
			Vertex* v2 = (*(*eit)).vertex2;
			double L = (*(*eit)).GetLength();
			*datafile << (*(*eit)).id << "\t" << (*v1).id << "\t" << (*v2).id << "\t" << (*(*eit)).xflag << "\t" << (*(*eit)).yflag << "\t" << (*(*eit)).GetTension_Magnitude() << "\t" << (*(*(*eit)).cell1).id << "\t" << (*(*(*eit)).cell2).id << endl;
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
		
		*datafile << "End" << endl << endl;
	}


	//*********************************************************************************************************************
	//******************* FUNCTIONS FOR PREFORMING TOPOLOGICAL CHANGES (Mostly Overloads) *********************************
	//*********************************************************************************************************************
	
/*
* bool Mitosis (overload)
*
* Do mitosis on cells. Picks cells for mitosis based on certain criteria and calls HelperMitosis(Cell * parentcell) to perform mitosis
* on those cells. This overloaded cell uses anisotropy potentials to determine the orientation of the cleavage plane.
*/
bool AGTissue::Mitosis(){
	vector<Cell*> mitcells;
	AGCell* ag1;
	int flag = 0;
	int i;
	double pressure;
	double area;
	double p0 = .01;
	
	double f = 0; //feedback strength - between 0 and 1
	// double alpha = p0 / (K*(A0 - f*A0)); //this defines feedback so that if a cell is smaller than f times the preferred area of cells, it can't divide
	double alpha = 0;

	double p;
	double angle; 

	//print Feedback Parameters
	/*
	if (CURRENTTIME == 1) {
		cout << "Mitosis Feedback Parameters:" << endl;
		cout << "p0: " << p0 << "		alpha: " << alpha << endl;
	}
	*/
	//Here is where we will later insert some criteria for mitosis
	for (list<Cell*>::iterator cit = cells.begin(); cit != cells.end(); cit++) {
		//division with feedback like normal
		area = (*(*cit)).GetArea();
		pressure = (*(*cit)).GetPressure(area);
		p = (p0 - alpha*pressure);
		if (p < 0)
		{ //if negative probability, then transfers over to similar procedure for apoptosis and p is set to 0 so mitosis will not occur
			//ApoptosisSelection(p, *cit); //comment this out to stop apoptosis from happening like this
			p = 0;
		}

		bernoulli_distribution distribution(p);        //distribution for weighting chance of success
		if (distribution(generator)) {
			mitcells.push_back(&(*(*cit)));
			flag++;
			//TOTALMITOSES++;
		}
	}
	if (flag == 0){
		return false;
	}
	else{
		//Output the cell id and area for all cells undergoing Mitosis
		// To do mitosis on only one cell, change for statement to i=0; i<(mitcells.size())-(mitcells.size()-1); i++, for normal, i = 0; i<mitcells.size(); i++ 
		for (i = 0; i<mitcells.size(); i++){ // changed for statement to only include i=0, so mitosis will only effect one cell
			//NOTE: only perform mitosis on cells that are convex in order to avoid cutting more than two edges.
			if (CheckConvex((*mitcells[i]).HelperOrderVerticesPos()) == true){
				//cout << "Mitosis of cell \t" << (*mitcells[i]).id << endl;
				//cout << "Cell Area \t" << (*mitcells[i]).GetArea() << endl;
				//cout << "Current time is \t" << CURRENTTIME << endl;
				ag1 = dynamic_cast<AGCell*>((mitcells[i]));
				angle = (*ag1).DivisionAngleDistribution((*ag1).DIVISIONPOTENTIAL);
				(*this).HelperMitosis(mitcells[i], angle);
			}
		}
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
Cell* AGTissue::CreateCell(int typ, vector<Edge*> edgs, string potential, double tcjstrength, double thetaan, double anstrength)
{
	int largestID, newid;
	list<Cell*>::iterator cit;
	//Get the id
	largestID = (*max_element(cellIDs.begin(), cellIDs.end()));			//first look through all of the current vertex ids and get the highest
	newid = largestID + 1;								//Now make a new ID that is unique

	if (newid > 1000000000) { cout << "ERROR to many cells" << endl; }
	//create the agcell and add the agcell to the lists
	AGCell* cptr = new AGCell;
	(*cptr).id = newid;
	(*cptr).type = typ;
	(*cptr).edges = edgs;
	(*cptr).TCJSTRENGTH = tcjstrength;
	(*cptr).DIVISIONPOTENTIAL = potential;
	(*cptr).ANSTRENGTH = anstrength;
	(*cptr).THETAAN = thetaan;
	(cells).push_back(cptr);
	cellIDs.push_back(newid);

	for (cit = cells.begin(); cit != cells.end(); cit++){
		if ((*(*cit)).id == newid){
			return &(*(*cit));
			break;
		}
	}
}

/*
* void HelperMitosis( Cell* parent cell, double divisionangle)
*
* Called by Mitosis to do the actual pointer and object manipulation required for mitosis 
*/
void AGTissue::HelperMitosis(Cell * parentcell, double divisionangle) {
	int i;										//for loops
	Vect center, p1, p2, q1, q2;				//Used to find newcenteredge
	Vect rnewvtx;								//the positions of the new vertices
	vector<Edge*> empty;
	vector<Edge*> deletededges;					//vector of edge pointers to the two edges to remove
	vector<Edge*> newedges;						//vector to hold the four new edges to replace the two that are being removed
	vector<Vertex*> newvertices;				//vector of pointers to the two new vertices
	vector<Vertex*> cwverts;					//parent cell verts in clockwise order
	Edge* newedgecenter = 0;					//ptr to the edge through the center of the cell
	Vertex *v1ptr, *v2ptr, *v3ptr, *v4ptr;
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
	AGCell* momptr = dynamic_cast<AGCell*>((parentcell));

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
		for (i = 0; i < (*parentcell).edges.size(); i++){								//look at each edge
			p1 = (*(*(*parentcell).edges[i]).vertex1).r;								//get it's end points and put them into p1 and p2
			p2 = (*(*(*parentcell).edges[i]).vertex2).r;
			if (test_for_Intersect(q1, q2, p1, p2, &rnewvtx)){							//if it crosses the selected angle
				deletededges.push_back((*parentcell).edges[i]);								//record the edge
				newvertices.push_back((*this).CreateVertex(rnewvtx, empty));				//make a new vertex at the interscetion and record its ptr
			}
		}

		//create a new edge between the two new vertices (must have no xflag and no yflag)
		newedgecenter = (*this).CreateEdge(newedgecenter, 0, 0, newvertices[0], newvertices[1], parentcell, parentcell);
	}

	//This is for the case where a cell crosses any boundary (WE LOOP THROUGH AND THEN PAIR THE END WITH THE BEGINNING)
	if ((*parentcell).GetFlag() == 1) {
		for (i = 0; i < vertnumb - 1; i++){ //look at each vertex
			p1.x = (vertposforflag[0])[i]; //get it's end points and put them into p1 and p2 (THESE POSITIONS CAN BE OUT OF BOUNDS)
			p1.y = (vertposforflag[1])[i];
			p2.x = (vertposforflag[0])[i + 1];
			p2.y = (vertposforflag[1])[i + 1];
			if (test_for_Intersect(q1, q2, p1, p2, &rnewvtx)){				//if it crosses the selected angle
				deletededges.push_back((*vertsordered[i]).GetEdgeBetween(vertsordered[i + 1]));								//record the edge			
				if (rnewvtx.x>LX || rnewvtx.x<0){ newxflag = 1; } //These six lines used later to create flags for newedgecenter (stores information about which side of boundary for each newvertex)
				else { newxflag = 0; }
				if (rnewvtx.y>LY || rnewvtx.y<0){ newyflag = 1; }
				else { newyflag = 0; }
				newflagstestx.push_back(newxflag);
				newflagstesty.push_back(newyflag);
				edgeflagstestind.push_back(i); //This whole sections helps us later when trying to set flags for other four newly created edges (TESTS WHICH SIDE OF THE BOUNDARY ARE THE FOUR VERTICES OF DELETED EDGES)
				if (p1.x>LX || p1.x < 0) { edgeflagstestx.push_back(1); }
				else { edgeflagstestx.push_back(0); }
				if (p1.y>LY || p1.y<0)  { edgeflagstesty.push_back(1); }
				else { edgeflagstesty.push_back(0); }
				edgeflagstestind.push_back(i + 1);
				if (p2.x>LX || p2.x < 0) { edgeflagstestx.push_back(1); }
				else { edgeflagstestx.push_back(0); }
				if (p2.y>LY || p2.y < 0)  { edgeflagstesty.push_back(1); }
				else { edgeflagstesty.push_back(0); } //Conclusion of section used for setting flags for the other four newly created edges
				if (rnewvtx.x>LX){ rnewvtx.x = rnewvtx.x - LX; } //This part puts the newvtx within the boundaries appropriately
				if (rnewvtx.x<0){ rnewvtx.x = rnewvtx.x + LX; }
				if (rnewvtx.y>LY){ rnewvtx.y = rnewvtx.y - LY; }
				if (rnewvtx.y<0){ rnewvtx.y = rnewvtx.y + LY; }
				newvertices.push_back((*this).CreateVertex(rnewvtx, empty));
			}
		}
		p1.x = (vertposforflag[0])[vertnumb - 1]; //get it's end points and put them into p1 and p2 FOR THE CASE WHERE THE BOTTOM ELEMENT IS PAIRED WITH THE TOP ELEMENT
		p1.y = (vertposforflag[1])[vertnumb - 1];
		p2.x = (vertposforflag[0])[0];
		p2.y = (vertposforflag[1])[0];
		if (test_for_Intersect(q1, q2, p1, p2, &rnewvtx)){				//if it crosses the selected angle
			deletededges.push_back((*(vertsordered[vertnumb - 1])).GetEdgeBetween(vertsordered[0]));								//record the edge
			if (rnewvtx.x>LX || rnewvtx.x<0){ newxflag = 1; } //These six lines used to store information about flags for newedgecenter (stores information about which side of boundary for each newvertex)
			else { newxflag = 0; }
			if (rnewvtx.y>LY || rnewvtx.y<0){ newyflag = 1; }
			else { newyflag = 0; }
			newflagstestx.push_back(newxflag);
			newflagstesty.push_back(newyflag);
			edgeflagstestind.push_back(vertnumb - 1);
			if (p1.x > LX || p1.x < 0) { edgeflagstestx.push_back(1); } //This whole sections helps us later when trying to set flags for four other newly created edges (TESTS WHICH SIDE OF THE BOUNDARY ARE THE FOUR VERTICES OF DELETED EDGES)
			else { edgeflagstestx.push_back(0); }
			if (p1.y>LY || p1.y<0)  { edgeflagstesty.push_back(1); }
			else { edgeflagstesty.push_back(0); }
			edgeflagstestind.push_back(0);
			if (p2.x>LX || p2.x < 0) { edgeflagstestx.push_back(1); }
			else { edgeflagstestx.push_back(0); }
			if (p2.y>LY || p2.y < 0)  { edgeflagstesty.push_back(1); }
			else { edgeflagstesty.push_back(0); } //Conclusion of section used for setting flags for four other newly created edges
			if (rnewvtx.x>LX){ rnewvtx.x = rnewvtx.x - LX; }//This part puts the newvtx within the boundaries appropriately
			if (rnewvtx.x<0){ rnewvtx.x = rnewvtx.x + LX; }
			if (rnewvtx.y>LY){ rnewvtx.y = rnewvtx.y - LY; }
			if (rnewvtx.y<0){ rnewvtx.y = rnewvtx.y + LY; }
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
	if (edgeflagstestind.size() == 4){ //This just makes sure that this only runs when the indices of the four vertices of deleted edges are well known
		if ((*newedgecenter).xflag == 1){ //if the new edge crosses the x bound
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
		if ((*newedgecenter).yflag == 1){ //if newedgecenter crosses the y boundary
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

	Edge e, e0, *eptr;
	Vertex v, v0;
	e0 = *newedgecenter;
	v0 = *(newvertices[0]);
	do{							//for the first daughter				
		for (i = 0; i < (*parentcell).edges.size(); i++){																					// look through all of the edges in the cell
			if (((*(*(*parentcell).edges[i]).vertex1) == v0 || (*(*(*parentcell).edges[i]).vertex2) == v0) && (*(*parentcell).edges[i]) != e0){		//find the edge that has v0 as a vertex and is not the edge e0 (the one you just looked at)
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
	do{
		for (i = 0; i < (*parentcell).edges.size(); i++){																					// look through all of the edges in the cell
			if (((*(*(*parentcell).edges[i]).vertex1) == v0 || (*(*(*parentcell).edges[i]).vertex2) == v0) && (*(*parentcell).edges[i]) != e0){		//find the edge that has v0 as a vertex and is not the edge e0 (the one you just looked at)
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
	(*daughter1ptr).type = (*parentcell).type + 1;
	(*daughter2ptr).type = (*parentcell).type + 1;

	// ^^^^^^^ THIS MIGHT NEED SOME CHANGING!

	//delete the parent cell
	(*this).DestroyCell(parentcell);
}

/*
* void Evolve( double interval, double *stepsize, ofstream* errorfile)
*
* Evolves the tissue one time interval taking steps no larger than dt and updates the vertex positions
*/
void AGTissue::Evolve(double interval, double *stepsize, ofstream* errorfile)
{
	//find the end point of this evolve
	double endtime = CURRENTTIME + interval;
	bool reinitialize = true;
	bool flag0, flag1, flag2, flag3, flag4;
	int count = 0;

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
		for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++) { (*(*vit)).rold = (*(*vit)).r; }
		LXold = LX;				//Record the old box size
		LYold = LY;

		//******CHOOSE YOUR INTEGRATOR******
		if (INTEGRATORNAME == "RK4_Adaptive"){
			//4th order Runge-Kutta wi th adaptive stepsize
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

		//Step 3: Update
		LX = LXnew;				//Update the current box size 
		LY = LYnew;
		(*this).HelperEvolve_UpdatePosition();  //Update the vertex positions (this automaticaly updates edge flags if a vertex has crossed the boundary)

		//STEP 4: look for topological changes and set flags if any have occured 

		if (APOPTOSISON) flag1 = (*this).Apoptosis(errorfile);

		if (T1ON)  flag2 = (*this).T1(errorfile);

		flag3 = (*this).ResolveVertices(errorfile);

		double timeratio = CURRENTTIME / TAU_STRESS;
		double floortr = floor(CURRENTTIME / TAU_STRESS);
		// This needs a proper calculate equilibrium function
		if (MITOSISON)
		{
			if (CURRENTTIME > TAU_STRESS && timeratio > floortr && timeratio < floortr + 2 * INITIALSTEPSIZE)
			{
				// Gives a half second window for mitosis
				flag4 = (*this).Mitosis();
				// THIS ACTIVATES MITOSIS!!!

				if (SIMTIMERES <= 0) (*this).WriteOutputDataAppend(errorfile);
				// if there is a non-zero "SIMULATION TIME RESOLUTION" inputted, don't write at every mitosis
			}
			else
			{
				flag4 = false;
			}
		}
		if (flag1 || flag2 || flag3 || flag4) reinitialize = true;  //Determine if we need to reinitialize
	}

	cout << "TISSUE EVOLVED TO TIME " << CURRENTTIME << " in " << count << " steps" << endl;
	cout << "LX = " << LX << "  LY = " << LY << endl;
}

/*
* void ForceFunction(ofstream* errorfile)
*
* Calculates the force on every vertex and stores the result in vertex.
*
* WARNING!!! This function assumes that cell.SavedInfoForEvolve_CWVerts has been updated correctly
*/
void AGTissue::ForceFunction(ofstream* errorfile){
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
		(*(*vit)).SavedInfoForEvolve_MechForce = zerovec; //first set the force to zero 
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
			if ((*v1).id == (*v2).id)
			{
				deltaP = (*(*eptr).cell2).SavedInfoForEvolve_Pressure - (*(*eptr).cell1).SavedInfoForEvolve_Pressure;
			}
			else
			{
				deltaP = (*(*eptr).cell1).SavedInfoForEvolve_Pressure - (*(*eptr).cell2).SavedInfoForEvolve_Pressure;
			}
			Pmag = -0.5 * deltaP * (*eptr).GetLength();   		//magnitude of the force due to the pressure difference
			angle = (*eptr).GetAngle((*vit)) + PI / 2;	 	//get the angle zXl
			Vect FP(zerovec, angle, Pmag);
			//Add the forces on this edge to the vertex
			(*(*vit)).SavedInfoForEvolve_MechForce = (*(*vit)).SavedInfoForEvolve_MechForce + FT + FP;
		}
	}

//Calculate the tissue stress
// This needs to be updated. Maybe get rid of GammaX/Y and figure out how to replace them
double GammaX = 0;
double GammaY = 0;
// calculateexternalstress_linear(TAU_STRESS);
/*
double cutlineX = distribution(generator) * LX; //Continue to generate random numbers until you cut no concave cells
GammaX = TISSUE_EXT_STRESS_X + (*this).CalculateTissueStressX(cutlineX, errorfile);

double cutlineY = distribution(generator) * LY; //Continue to generate random numbers until you cut no concave cells
GammaY = TISSUE_EXT_STRESS_Y + (*this).CalculateTissueStressY(cutlineY, errorfile);
*/
// new code, 5/20/2020
// p_ext is external pressure
double p_ext = (TISSUE_EXT_STRESS_X + TISSUE_EXT_STRESS_Y) / 2;

GammaX = -1 * (LY*LX0) * ( p_ext + GetInternalTissueStressX(errorfile) ); // + internalstress (placeholder)
GammaY = -1 * (LY0*LX) * ( p_ext + GetInternalTissueStressY(errorfile) );

//Update LY and LX
ForceLX = LX0 * GammaX;
ForceLY = LY0 * GammaY;

//Add the force from the tissue stress to every vertex
//Added the length squared scaling in this step
for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
{
	(*(*vit)).SavedInfoForEvolve_TotalForce.x = LX * LX * (*(*vit)).SavedInfoForEvolve_MechForce.x + (alphaM / TISSUE_ALPHA) * ((*(*vit)).r.x / LX) * ForceLX;
	(*(*vit)).SavedInfoForEvolve_TotalForce.y = LY * LY * (*(*vit)).SavedInfoForEvolve_MechForce.y + (alphaM / TISSUE_ALPHA) * ((*(*vit)).r.y / LX) * ForceLY;
}
}

/*
* void CalculateTissueStressX(ofstream* errorfile)
*
* Calculates the x-component of the external stress in the tissue.
* 1. Adds internal stresses that are independent of the boundary
* 2. Finds crossing edges and crossing cells
* 3. Finds contributions to stress from crossing edges (happens as crossing edges are identified)
* 4. Calculate area contribution to stress
*	4.1 Loop over crossing cells
*   4.2 Identify crossing edges specific to each cell
*   4.3 Calculate dA/dLx using said crossing edges
* 5. stress = boundary_independent_stress + areacontribution + edgecontribution + A0 *Ly * Lx 
* See Avik's notes for derivation of all this
*/
double AGTissue::GetInternalTissueStressX(ofstream* errorfile)
{
	double stressx = 0;
	// total stress on tissue
	double edgestress = 0;
	// contribution only from edges
	double areastress = 0;
	// contribution from cell areas
	double dAdLx = 0;
	// derivative of the cell area with respect to Lx
	bool xflag, c1, c2;
	double directioncos, angle;
	Vertex *v1;
	vector<Edge*> crossingEdges;									//vector of all crossing edges
	vector<Edge*> cellEdges;										// vector of edges for a specific cell
	vector<Cell*> crossingCells;									//vector of all crossing cells
	int i;
	Cell *cptr;
	double length;
	vector<double> edgecenter;

	// add in the contributions from within the tissue to the total stress
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
	{
		stressx += (*(*vit)).SavedInfoForEvolve_MechForce.x * (*(*vit)).r.x;
	}

	// find crossing edges and crossing cells
	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++)
	{
		xflag = (*(*eit)).xflag;
		if (xflag)
		{
			crossingEdges.push_back(&(*(*eit)));
			//if the cells are not in the list of crossing cells already add them 
			long int c1id = (*(*(*eit)).cell1).id;
			long int c2id = (*(*(*eit)).cell2).id;
			c1 = false;
			c2 = false;
			for (i = 0; i < crossingCells.size(); i++)
			{
				if ((*crossingCells[i]).id == c1id) c1 = true;
				if ((*crossingCells[i]).id == c2id) c2 = true;
			}
			if (!c1) crossingCells.push_back(&(*(*(*eit)).cell1));
			if (!c2) crossingCells.push_back(&(*(*(*eit)).cell2));

			// here, I'll add the edge contribution to the stress
			v1 = (*(*eit)).vertex1;
			angle = (*(*eit)).GetAngle(v1);
			directioncos = abs(cos(angle));
			edgestress += GAMMA * directioncos ;
		}
	}

	for (i = 0; i < crossingCells.size(); i++)
	{
		cellEdges = (*crossingCells[i]).edges;
		 
		for (vector<Edge*>::iterator eit = cellEdges.begin(); eit != cellEdges.end(); eit++)
		{
			if ((*(*eit)).xflag)
			{ //crossing edge if statement
				edgecenter.push_back( (*(*eit)).GetCenter().y );
			} // crossing edge if statement
		}
		dAdLx = abs(edgecenter[1] - edgecenter[2]);
		areastress += (*crossingCells[i]).GetArea() * K * dAdLx;

		edgecenter.pop_back();
		edgecenter.pop_back();
		// empty the averageedgecenter list instead of making a new one
		if (!edgecenter.empty())
		{
			cout << "ERROR: Cell with more than two crossing edges \n";
			*errorfile << "Found cell with more than two crossing edges: " << (*crossingCells[i]).id << "\n";
			return NULL;
		}
	} //end for each cell in crossingCells
	stressx += ( edgestress + areastress - K * LX0 * LY0 * LY ) * LX;
	return stressx;
}
// PBC handling code for crossing edges - possibly made redundant by Meryl's GetCenter in EdgeClass
/*********************************************************************************************************************
v3 = (*(*eit)).vertex1;
v4 = (*(*eit)).vertex2;
if ((*v3).r.x > (*v4).r.x) // look at edges that cross x = Lx, not x = 0
{
y3 = (*v3).r.y;
y4 = (*v4).r.y;
}
else
{
y4 = (*v3).r.y;
y3 = (*v4).r.y;
}
if ((*(*eit)).yflag) // "undo" pbc
{
if (y3 > y4)
{
y4 = y4 + LY;
}
else if (y3 < y4)
{
y4 = y4 - LY;
}
}
edgecenter.push_back( (y4 + y3) / 2.0 );
*/ 
// ******************************************************************************************************************

/*
* void CalculateTissueStressY(ofstream* errorfile)
*
* Calculates the y-component of the external stress in the tissue.
*/
double AGTissue::GetInternalTissueStressY(ofstream* errorfile)
{
	double stressy = 0;
	// total stress on tissue
	double edgestress = 0;
	// contribution only from edges
	double areastress = 0;
	// contribution from cell areas
	double dAdLy = 0;
	// derivative of the cell area with respect to Ly
	bool yflag, c1, c2;
	double directionsin, angle;
	Vertex *v1;
	vector<Edge*> crossingEdges;									//vector of all crossing edges
	vector<Edge*> cellEdges;										// vector of edges for a specific cell
	vector<Cell*> crossingCells;									//vector of all crossing cells
	int i;
	Cell *cptr;
	double length;
	vector<double> edgecenter;

	// add in the contributions from within the tissue to the total stress
	for (list<Vertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
	{
		stressy += (*(*vit)).SavedInfoForEvolve_MechForce.y * (*(*vit)).r.y;
	}

	// find crossing edges and crossing cells
	for (list<Edge*>::iterator eit = edges.begin(); eit != edges.end(); eit++)
	{
		yflag = (*(*eit)).yflag;
		if (yflag)
		{
			crossingEdges.push_back(&(*(*eit)));
			//if the cells are not in the list of crossing cells already add them 
			long int c1id = (*(*(*eit)).cell1).id;
			long int c2id = (*(*(*eit)).cell2).id;
			c1 = false;
			c2 = false;
			for (i = 0; i < crossingCells.size(); i++)
			{
				if ((*crossingCells[i]).id == c1id) c1 = true;
				if ((*crossingCells[i]).id == c2id) c2 = true;
			}
			if (!c1) crossingCells.push_back(&(*(*(*eit)).cell1));
			if (!c2) crossingCells.push_back(&(*(*(*eit)).cell2));

			// here, I'll add the edge contribution to the stress
			v1 = (*(*eit)).vertex1;
			angle = (*(*eit)).GetAngle(v1);
			directionsin = abs(sin(angle));
			edgestress += GAMMA * directionsin;
		}
	}

	for (i = 0; i < crossingCells.size(); i++)
	{
		cellEdges = (*crossingCells[i]).edges;

		for (vector<Edge*>::iterator eit = cellEdges.begin(); eit != cellEdges.end(); eit++)
		{
			if ((*(*eit)).yflag)
			{ //crossing edge if statement
				edgecenter.push_back((*(*eit)).GetCenter().x);
			} // crossing edge if statement
		}
		dAdLy = abs(edgecenter[1] - edgecenter[2]);
		areastress += (*crossingCells[i]).GetArea() * K * dAdLy;

		edgecenter.pop_back();
		edgecenter.pop_back();
		// empty the averageedgecenter list instead of making a new one
		if (!edgecenter.empty())
		{
			cout << "ERROR: Cell with more than two crossing edges \n";
			*errorfile << "Found cell with more than two crossing edges: " << (*crossingCells[i]).id << "\n";
			return NULL;
		}
	} //end for each cell in crossingCells
	stressy += (edgestress + areastress - K * LX0 * LY0 * LX) * LY;
	return stressy;
}