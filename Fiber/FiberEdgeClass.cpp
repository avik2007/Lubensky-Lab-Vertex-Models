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
#include "Constants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "FiberCellClass.h"
#include "FiberConstants.h"
#include "FiberClass.h"
#include "FiberEdgeClass.h"
#include "FiberTissueClass.h"

using namespace std;

//**************************************************************************************************************************************
//********************************************FIBER EDGE*************************************************************************
//**************************************************************************************************************************************

//Returns true if the edge has a fiber
bool FiberEdge::HasFiber(){
	if (fibers.size() > 0){ return true; }
	else{ return false; }
}

//Removes the fiber from the edge
void FiberEdge::RemoveFiber(Fiber* f){
	int i;
	for (i = 0; i<fibers.size(); i++){
		if ((*fibers[i]).id == (*f).id){
			fibers.erase(fibers.begin() + i);
		}
	}
}
