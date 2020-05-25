// guard
#ifndef __AGTISSUE_H_INCLUDED__
#define __AGTISSUE_H_INCLUDED__

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
#include "AGConstants.h"
#include "EdgeClass.h"
#include "CellClass.h"
#include "AGCell.h"
#include "VertexClass.h"
#include "TissueClass.h"

using namespace std;

class AGTissue : public Tissue{
public:
	//ADDITIONAL DATA(no additional data to add here)									// deformation gradient tensor component (22)
	double LX0;									// initial box length in x dimension
	double LY0;									// initial box length in y dimension

	//NOTE: the list cells that is instantiated in the TissueClass will contain AGCells when AGVertex Model code
	// is being used 
	//CONSTRUCTOR
	AGTissue(){};

	//DESTRUCTOR
	virtual ~AGTissue(){};

	//ADDITIONAL FUNCTIONS

	// Updates read input function to accept new global variables 
	virtual void HelperReadInputData_SetGlobalValues(ifstream* datafile);

	//Updates read input function to create new objects (AGCells)
	virtual void HelperReadInputData_CreateObjects(ifstream* datafile);

	//! Writes data into file specified by 'datafile'
	virtual void HelperWriteOutputData_Write2File(ofstream* datafile, bool includeExtraInfo);

	//*********************************************************************************************************************
	//******************* FUNCTIONS FOR PREFORMING TOPOLOGICAL CHANGES (Mostly Overloads) *********************************
	//*********************************************************************************************************************
	// Do mitosis on cells. Picks cells for mitosis based on certain criteria and calls HelperMitosis(Cell * parentcell) to perform mitosis on those cells.
	virtual bool Mitosis();

	//virtual void DestroyCell(Cell* cptr);
	virtual Cell* CreateCell(int typ, vector<Edge*> edgs, string potential, double tcjstrength, double thetaan, double anstrenght);

	//! This is the helper function that actually do mitosis on cells.
	virtual void HelperMitosis(Cell * parentcell, double divisionangle);

	//*********************************************************************************************************************
	//***************** FUNCTIONS FOR PERFORMING NUMERICAL INTEGRATION TO EVOLVE THE TISSUE ******* ************************
	//*********************************************************************************************************************

	//! Evolves the tissue one time interval taking steps no larger than dt and updates the vertex positions
	virtual void Evolve(double interval, double *stepsize, ofstream* errorfile);

	//! Calculates the force on every vertex and stores the result in vertex.
	virtual void ForceFunction(ofstream* errorfile);

	//! Calculates the x component of the internal stress of the tissue (so sigma_xx)
	virtual double GetInternalTissueStressX(ofstream* errorfile);

	//! Calculates the y component of the internal stress of the tissue (so sigma_yy)
	virtual double GetInternalTissueStressY(ofstream* errorfile);

};
#endif