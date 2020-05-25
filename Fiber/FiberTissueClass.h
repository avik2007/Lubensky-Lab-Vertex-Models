// include guard
#ifndef __FIBERTISSUECLASS_H_INCLUDED__
#define __FIBERTISSUECLASS_H_INCLUDED__

//Forward declared dependences
class Fiber;
class FiberCell;
class FiberEdge;
class Tissue; 
class Cell;
class Edge;
class Vertex;
class Vect;

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>


#include "VectClass.h"
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h" 
#include "FiberClass.h"

using namespace std;

struct fiberstats {
	int id;
	double lifetime;
	int cellid;
	int cellsize; 
	double cellarea;
	bool passedvertex;
};

//! Class for tissue object.
class FiberTissue : public Tissue{
public:
	//ADDITIONALDATA
	list<NucleationSite*> fibernucsites; 
	list<long int> NucleationSiteIDs;						//lists of all of the ID lables ever used for the vertices edges and cells so new lables can be made unique
	list<Fiber*> fibers;								//!< List of all fiber objects (not pointers!) in the simulation. This is where the real objects are stored.			
	list<long int> fiberIDs;						//lists of all of the ID lables ever used for the vertices edges and cells so new lables can be made unique
	vector<fiberstats> fiberstatistics;

	//FUNCTIONS
	//CONSTRUCTOR
	FiberTissue(){
		fiberIDs.push_back(0);  //start the list if ids so new ones can be generated for future fibers
		NucleationSiteIDs.push_back(0);  //start the list if ids so new ones can be generated for future fibers
	};							
	virtual ~FiberTissue(){};

	// Reads the data and sets all of the global variables
	virtual void HelperReadInputData_SetGlobalValues(ifstream* datafile); 

	//! Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids. Starts at current position of *datafile and reads until keyword 'End'.
	virtual void HelperReadInputData_CreateObjects(ifstream* datafile);

	//Creates a fiber cell and adds it to the tissue
	Cell* CreateCell(int typ, vector<Edge*> edgs);

	//Create a fiber vertex
	NucleationSite* CreateNucleationSite(Vertex* parentv, FiberEdge* e1, FiberEdge* e2, FiberCell* c);

	//! Takes a tissue and creates nucleation sites from the parent verticies and sets their pointers 
	void CreateNucleationSites();
	
	//Creates a fiber and adds it to the tissue
	Fiber* CreateFiber(NucleationSite* v);

	//Removes a fiber from the tissue
	void DestroyFiber(Fiber* f);

	//writes the output data for the movie
	void WriteOutputData(ofstream* datafile);

	//writes aditional infor to file like fibernumber v cell area
	void WriteOutputStatistics(ofstream* datafile); 

	//updates the fiber f and returns true if the fiber has changed topology
	bool UpdateFiber(Fiber* f); 

	//updates all fibers and returns true if any fiber has changed topology
	bool UpdateFibers();

	//Writes the current state of the system as an adajacency list to run graph algorythems on 
	void WriteAsGraph(ofstream* digraphfile); 

	//Initializes the tissue after it is read in (removes dangeling edges and cells ...)
	void InitializeFromData();

	//Reads in the point that jesus clicked and sprouts a fiber where they are
	bool ReadFiberData(string filename);

	//Reads a file that includes fiber info at the end
	void ReadInputFibers(string filename); 

	//Evolves the fibers by giving them some velocity 
	void EvolveFibers(double dt, double tstep); 

	//Puts the cells into small or large catagories 
	void RecordCellSize();

};

#endif