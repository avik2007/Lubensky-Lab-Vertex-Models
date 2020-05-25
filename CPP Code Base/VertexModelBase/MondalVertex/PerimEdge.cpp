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
#include <numeric>
#include <stddef.h>

#include "Constants.h"
#include "ConstantsPerim.h"
#include "PerimEdge.h"
#include "CellClass.h"
#include "PerimCell.h"


//Returns the magnitude of the tension on the edge 
//Should return gamma+Lambda(cell1.P0-cell1.perimeter) + Lambda(cell2.P0- cell2.perimeter)

double Edge_Perim::GetTension_Magnitude(){
	
	double perimdiff1; //difference between the perimeter and P0 for cell1
	double perimdiff2; //difference between the perimeter and P0 for cell2
	double tension; 

	//Get pointers to the neighboring cells (note that we need to cast the cell as the appropriate cell subclass) 
	Cell_Perim* c1 = dynamic_cast<Cell_Perim*>(cell1);
	Cell_Perim* c2 = dynamic_cast<Cell_Perim*>(cell2);
	/*
	if (LAMBDA != 0) {
		
		//Calculate P_i - P_0
		perimdiff1 = (*c1).GetPerimeter() - (*c1).P0;
		perimdiff2 = (*c2).GetPerimeter() - (*c2).P0;

		// changing up this formula based on Manning theory in paper ("A density-independent rigidity transition in biological tissues");
		// THIS WILL NEED A VERY DETAILED EXPLANATION WRITTEN AT SOME POINT TO SUPPLEMENT MERYL'S THESIS 
		tension = LAMBDA*perimdiff1 + LAMBDA*perimdiff2;
		
	}
	else {
		tension = GAMMA;
		perimdiff1 = 0;
		perimdiff2 = 0;
	}
	*/
	
	// Meryl's original code 
	//Calculate P_i - P_0
	perimdiff1 = (*c1).GetPerimeter() - (*c1).P0;
	perimdiff2 = (*c2).GetPerimeter() - (*c2).P0;
	
	//Calculate the tension
	tension = GAMMA + LAMBDA*perimdiff1 + LAMBDA*perimdiff2;
	
	return tension; 
}
