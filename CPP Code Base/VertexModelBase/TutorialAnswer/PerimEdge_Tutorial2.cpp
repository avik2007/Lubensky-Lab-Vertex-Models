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
#include "PerimEdge.h"
#include "CellClass.h"
#include "PerimCell.h"


//Returns the magnitude of the tension on the edge 
//Should return gamma+Lambda(cell1.P0-cell1.perimeter) + Lambda(cell2.P0- cell2.perimeter
double Edge_Perim::GetTension_Magnitude(){
	double perimdiff1; //difference between the perimeter and P0 for cell1
	double perimdiff2; //difference between the perimeter and P0 for cell2
	double tension; 

	//Get pointers to the neighboring cells 
	Cell_Perim* c1 = dynamic_cast<Cell_Perim*>(cell1); //pointer to cell1 (note that we need to cast the cell as the appropreate cell subclass) 
	Cell_Perim* c2 = dynamic_cast<Cell_Perim*>(cell2);

	//Calculate (cell.P0- cell.perim)
	perimdiff1 = (*c1).P0 - (*c1).GetPerimeter(); 
	perimdiff2 = (*c2).P0 - (*c2).GetPerimeter(); 

	//Calculate the tension 
	tension = GAMMA + LAMBDA*perimdiff1 + LAMBDA*perimdiff2; 
	return tension; 
}
