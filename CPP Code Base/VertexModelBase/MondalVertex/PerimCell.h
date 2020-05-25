// include guard
#ifndef __PERIMCELL_H_INCLUDED__
#define __PERIMCELL_H_INCLUDED__

//Forward declared dependences
class Cell;


#include <stddef.h>
#include <unordered_map>

#include "CellClass.h"

using namespace std;

//! Class for cell objects.
class Cell_Perim : public Cell{
public:
	//ADDITIONAL DATA
	double P0;

	//CONSTRUCTORS AND DESTRUCTORS
	Cell_Perim(){};
	virtual ~Cell_Perim(){};
	Cell_Perim(int idnumb, vector<Edge*> edglist, int typeinput, double A0input, double P0input){
		id = idnumb;
		edges = edglist;
		type = typeinput;
		A0 = A0input;
		P0 = P0input;
	}

};
#endif