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
#include "AGPerimTissue.h"
#include "AGPerimConstants.h"
#include "AGPerimCell.h"
#include "AGPerimEdge.h"
#include "VectClass.h"
#include "AGStiffCell.h"
#include "AGStiffTissue.h"
#include "AGStiffConstants.h"


/*
* void SetAreaModulus(double K)
*
* designed to be used by the tissue to set the area modulus of the AGStiffCells
*/
void AGStiffCell::SetAreaModulus(double K)
{
	(*this).KAM = K;
}

/*
* void SetPreferredArea(double PA)
*
* designed to be used by the tissue to set the preferred area of the AGStiffCells
*/
void AGStiffCell::SetPreferredArea(double PA)
{
	(*this).A0 = PA;
}

/*
* void SetPreferredArea(double PM)
*
* designed to be used by the tissue to set the perimeter modulus of the AGStiffCells
*/
void AGStiffCell::SetPM(double PM)
{
	(*this).PM = PM;
}

/*
* double GetPressure(double area)
*
* returns the pressure of AGStiffCell
*/
double AGStiffCell::GetPressure(double area)
{
	double pressure = (*this).KAM * ((*this).A0 - area);
	return pressure;
}

/*
* double GetEnergy()
*
* returns the energy of an individual AGStiffCell
*/
double AGStiffCell::GetEnergy()
{
	double energy = 0.5 * ((*this).KAM * pow((*this).GetArea() - (*this).A0, 2)  + LAMBDA * pow((*this).GetPerimeter() - P0, 2));
	return energy;
}
/*
* double GetPressureNormalStress();
*
* Calculates the normal stress due to pressure in the model
*/
double AGStiffCell::GetPressureNormalStress()

{
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	// clockwiseVertex.push_back(clockwiseVertex[0]);
	Vertex* v1, * v2;
	double F1, F2;
	double x1, x2;
	double stressxx = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		// v2 = clockwiseVertex[i + 1];
		x1 = (*v1).r.x;
		F1 = (*v1).SavedInfoForEvolve_TotalForce.x;
		// F2 = (*v1).SavedInfoForEvolve_TotalForce.x;
		stressxx += (F1 * x1);
	}
	return stressxx / (*this).GetArea();
}

/*
* double GetTensionNormalStress();
*
* Calculates the normal stress due to tension in the model
*/
double AGStiffCell::GetTensionNormalStress()
{
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	Vertex* v1, * v2;
	double F1, F2;
	double y1, y2;
	double stressyy = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		y1 = (*v1).r.y;
		F1 = (*v1).SavedInfoForEvolve_TotalForce.y;
		stressyy += (F1 * y1);
	}
	return stressyy / (*this).GetArea();
}