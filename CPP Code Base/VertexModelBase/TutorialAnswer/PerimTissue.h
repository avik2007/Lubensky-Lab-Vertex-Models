// guard
#ifndef __PERIMTISSUE_H_INCLUDED__
#define __PERIMTISSUE_H_INCLUDED__

// Forward declared dependencies
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
#include <stddef.h>

#include "VectClass.h"
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h"
#include "PerimEdge.h"
#include "PerimCell.h"

using namespace std;

class Tissue_Perim : public Tissue{
public:
	//ADDITIONAL DATA(no additional data to add here)

	//CONSTRUCTOR
	Tissue_Perim(){};

	//DESTRUCTOR
	virtual ~Tissue_Perim(){};

	//ADDITIONAL FUNCTIONS

	//! Reads and constructs simulation from file. The file must be formatted like InputTemplate.txt. This function parses the input file then initializes all the objects and pointers.
	bool ReadInputData(string InputFileName);

	//Updates read input function to read in new variables P0 and LAMBDA
	void HelperReadInputData_SetGlobalValues(ifstream* datafile);

	//! Make all of the objects with the correct starting positions and flags and recordds their ids in the tissue's list of ids. Starts at current position of *datafile and reads until keyword 'End'.
	void HelperReadInputData_CreateObjects(ifstream* datafile);
};

#endif