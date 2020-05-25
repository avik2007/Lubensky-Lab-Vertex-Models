
// include guard
#ifndef __FIBERCELLCLASS_H_INCLUDED__
#define __FIBERCELLCLASS_H_INCLUDED__

//Forward declared dependences
class Fiber;
class FiberCell;
class FiberEdge;
class NucleationSite; 
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
#include <random>
#include <unordered_map>
#include "Constants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "TissueClass.h" 

using namespace std;


//! Class for cell objects.
class FiberCell : public Cell{
public:
	//ADDITIONAL DATA
	vector<Fiber*> fibers; 
	vector<NucleationSite*> fibernucsites; 
	bool isCompleete; 
	bool isBoundary; 
	int sizeclass; // 1 is small 2 is large (0 is incompleete/boundary cell) 
	int number_of_fibers_from_data;

	//ADDITIONAL FUNCTIONS

	//CLONE
	FiberCell* Clone(){
		return new FiberCell(*this);
	}

	//CONSTRUCTORS AND DESTRUCTORS
	FiberCell(){};
	virtual ~FiberCell(){};
	FiberCell(int idnumb, vector<Edge*> edglist, int typeinput, double A0input){ id = idnumb; edges = edglist; type = typeinput; A0 = A0input; number_of_fibers_from_data = 0; }

	//Functions

	//! Removes fiber from the cell.
	void RemoveFiber(Fiber* fptr);

	//returns a 2 element vector where the first number is the left side of the cell and the second is the right side obeying pbc
	vector<double> GetXExtent();

	//returns a 2 element vector where the first number is the bottom of the cell and the second is the top of the cell obeying pbc
	vector<double> GetYExtent(); 

	//Sets the isCompleete flag
	bool SetIsCompleeteFlag(); 

	//takes and edge and a vertex on that edge and returns true if that vertex has another edge in this cell
	bool GetNextEdge(Edge* e, Vertex* v, int* result);

	//calculate the number of nucleating sites with appropriate opening angle.
	double CountNucleatingSites(); 

	double CountFiberDeathSites();

	//returns the max horizontal extent of the cell
	double GetMaxHorizontalLength();

	double GetMaxVerticalLength();

	double CalcAverageXDist_toDeathVertex();
};

double inverseGaussian(double mu, double lambda, mt19937_64* generator);

#endif