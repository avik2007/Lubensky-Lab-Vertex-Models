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

/*
* virtual double GetEnergy()
*
* returns the energy of the cell
*
*/
double AGPerimCell::GetEnergy()
{
	double energy = 0.5 * ( K * pow((*this).GetArea() - (*this).A0 ,2) / (*this).A0  + LAMBDA * pow( (*this).GetPerimeter() - P0,2) );
	return energy;
	// before using this again you have to update to account for possible log or Gaussian area potentials
}
/*
* double GetdPdLx()
*
* Gets the derivative of a cell's perimeter with respect to the position of the right box wall	
*/
double AGPerimCell::GetdPdLx()
{
	int xcrossingEdges = 0;
	vector<Vect> edgeVects;
	double dPdLx = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double angle = 0;
	double directioncos = 0;
	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVects.push_back((*e1).GetVectorDifference(v1, v2));
		
		if ( ((*e1).xflag) )
		{
			angle = (*e1).GetAngle(v1);
			directioncos = abs(cos(angle));
			dPdLx += directioncos;
			if (xcrossingEdges < 1)
			{
				xcrossingEdges++;
			}
			else
			{
				i = vertexNum;
			}

		}
	}
	return dPdLx
}

double AGPerimCell::GetdPdLy()
{
	int ycrossingEdges = 0;
	vector<Vect> edgeVects;
	double dPdLy = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double angle = 0;
	double directionsin = 0;
	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVects.push_back((*e1).GetVectorDifference(v1, v2));

		if (((*e1).xflag))
		{
			angle = (*e1).GetAngle(v1);
			directionsin = abs(sin(angle));
			dPdLy += directionsin;
			if (ycrossingEdges < 1)
			{
				ycrossingEdges++;
			}
			else
			{
				i = vertexNum;
			}

		}
	}
	return dPdLy;
}

double AGPerimCell::Getd2PdLxdy(Vertex* V)
{
	unsigned long int targetvertexid = (*V).id;
	bool id_match = false;
	int ycrossingEdges = 0;
	vector<Vect> edgeVects;
	double d2PdLxdy = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double angle = 0;
	double directionsin = 0;

	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVects.push_back((*e1).GetVectorDifference(v1, v2));
		id_match = ((v1->id) == targetvertexid) || ((v2->id) == targetvertexid);

		if (((*e1).xflag))
		{
			angle = (*e1).GetAngle(v1);
			directionsin = abs(sin(angle));
			d2PdLxdy += directionsin;
			if (ycrossingEdges < 1)
			{
				ycrossingEdges++;
			}
			else
			{
				i = vertexNum;
			}

		}
	}
	return d2PdLxdy;
}

double AGPerimCell::GetdPdgammaForSimpleShearAlongXAxis()
{

}

double AGPerimCell::Getd2PdgammadyForSimpleShearAlongXAxis(Vertex* v)
{

}