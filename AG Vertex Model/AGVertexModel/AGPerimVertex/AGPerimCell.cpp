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
	return dPdLx;
}
/*
* double Getd2PdLx2()
*
* (7/9/24): DOUBLE CHECK THIS PLEASE!!!!
* Gets the derivative of a cell's perimeter with respect to the position of the right box wall
*/
double AGPerimCell::Getd2PdLx2()
{
	int xcrossingEdges = 0;
	vector<Vect> edgeVects;
	double d2PdLx2 = 0;
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
			d2PdLx2 += 1 / (*e1).GetLength() * pow(directionsin,2); 
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
	return d2PdLx2;
}
/*
* double GetdPdLy()
*
* * (7/9/24): DOUBLE CHECK THIS PLEASE!!!!
* Gets the derivative of a cell's perimeter with respect to the position of the upper box wall
*/
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

if (((*e1).yflag))
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
/*
* double Getd2PdLy2()
*
* * (7/9/24): DOUBLE CHECK THIS PLEASE!!!!
* Gets the derivative of a cell's perimeter with respect to the position of the right box wall
*/
double AGPerimCell::Getd2PdLy2()
{
	int ycrossingEdges = 0;
	vector<Vect> edgeVects;
	double d2PdLy2 = 0;
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

		if (((*e1).yflag))
		{
			angle = (*e1).GetAngle(v1);
			directioncos = abs(cos(angle));
			d2PdLy2 += 1 / (*e1).GetLength() * pow(directioncos, 2);
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
	return d2PdLy2;
}
/*
* Vect Getd2PdLxdy(Vertex* V)
*
* derivative of perimeter with respect to Lx and to the y of a vertex V
* *
* *
* *
*  MAKE SURE YOU RECHECK THIS ONE * (7/9/24): DOUBLE CHECK THIS PLEASE!!!!
*  THIS DOESN'T ACCOUNT FOR THE POTENTIAL FOR SHEARING
* *
* * AS OF 9/29/24 THIS IS ABSOLUTELY WRONG!
*/
Vect AGPerimCell::Getd2PdLxdv(Vertex* V)
{
	unsigned long int targetvertexid = (*V).id;
	bool id_match = false;
	int xcrossingEdges = 0;
	double d2PdLxdy = 0;
	double d2PdLxdx = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	Vect v;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double directionsin = 0;
	double directioncos = 0;

	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		if ((v1->id) == targetvertexid || (v2->id) == targetvertexid)
		{
			e1 = (*v1).GetEdgeBetween(v2);
			if ((v1->id) == targetvertexid)
			{
				v = (*e1).MakeUnitVector(v1);
			} 
			else
			{
				v = (*e1).MakeUnitVector(v2);
			}
			
			clockwiseEdges.push_back(e1);

			if (((*e1).xflag))
			{
				// directioncos = abs(v.x); // these absolute values are part of a convention needed to get this derivative right I believe
				// directionsin = abs(v.y);
				// (07/09/2024 update: I think the abs() needs to be removed, we commit to perimeter being calculated by traversing edges CCW (add a minus sign to everything))
				directioncos = -1 * v.x; 
				directionsin = -1 * v.y;
				d2PdLxdy += 1 / (*e1).GetLength() * (-1* directioncos * directionsin );
				d2PdLxdx += 1 / (*e1).GetLength() * (directioncos * directioncos); //(1 - directionsin * directionsin ) ;
				if (xcrossingEdges < 1)
				{
					xcrossingEdges++;
				}
				else
				{
					i = vertexNum-1;
				}

			}
		}
	}
	Vect d2PdLxdv(d2PdLxdx, d2PdLxdy);
	return d2PdLxdv;
}

/*
* double GetdPdgamma()
*
* derivative of perimeter of a cell with respect to gamma, where gamma is the shear parameter for Lee Edwards Boundary conditions (see Merkel Manning 2018). We shear the unit cells in the x-direction
* * (7/9/24): DOUBLE CHECK THIS PLEASE!!!!
*/
double AGPerimCell::GetdPdgammaForSimpleShearAlongXAxis()
{
	int ycrossingEdges = 0;
	Vect edgeVect;
	double dPdgamma = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double angle = 0;
	double directioncos = 1;
	clockwiseVertex.push_back(clockwiseVertex[0]);
	int updownfactor = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVect = (*e1).GetVectorDifference(v2, v1);

		if (((*e1).yflag))
		{
			// this calculation requires the designation of the top and bottom vertex. We are considering edges that
			// cross the top boundary. The "top" vertex is in the unit cell above the primary unit cell (the one we 
			// visualize). The "bottom" vertex is in the primary unit cell. 

			// if the x coordinate of the top vertex is less than that of the bottom vertex, 
			directioncos = edgeVect.x / edgeVect.Magnitude();
			if (edgeVect.y < 0)
			{ 
				updownfactor = -1;
			}
			else
			{
				updownfactor = 1;
			}
		
			//directioncos = edgeVect.x / edgeVect.Magnitude();
			dPdgamma += LY * updownfactor * directioncos;
			if (ycrossingEdges < 1)
			{
				ycrossingEdges++;
			}
			else
			{
				i = vertexNum-1;
			}
		}
	}
	return dPdgamma;
}
/*
* double GetdPdgamma()
*
* second derivative of perimeter of a cell with respect to gamma, where gamma is the shear parameter for Lee Edwards Boundary conditions (see Merkel Manning 2018). We shear the unit cells in the x-direction
* * (7/9/24): DOUBLE CHECK THIS PLEASE!!!!
*/
double AGPerimCell::Getd2Pdgamma2ForSimpleShearAlongXAxis()
{
	int ycrossingEdges = 0;
	vector<Vect> edgeVects;
	double d2Pdgamma2 = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double angle = 0;
	double directionsin = 0;
	clockwiseVertex.push_back(clockwiseVertex[0]);
	double edgelength = 0;
	Vect evec;
	for (int i = 0; i < vertexNum ; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		evec = (*e1).GetVectorDifference(v2, v1);
		if (((*e1).yflag))
		{
			//angle = (*e1).GetAngle(v1);
			edgelength = (*e1).GetLength();
			// directionsin = abs(sin(angle));
			// (07/09/2024 update: I think the abs() needs to be removed, we commit to perimeter being calculated by traversing edges CCW (add a minus sign to everything))
			directionsin = abs(evec.y) / edgelength; // the x-magnitude of the edge scales this derivative, must be positive
			d2Pdgamma2 += LY * LY * pow(directionsin,2) / edgelength;
			ycrossingEdges++;
			if (ycrossingEdges > 1)
			{
				i = vertexNum - 1;
			}
		}
	}
	
	return d2Pdgamma2;
}
Vect AGPerimCell::Getd2PdgammadvForSimpleShearAlongXAxis(Vertex* v)
{
	unsigned long int targetvertexid = (*v).id;
	bool id_match = false;
	int ycrossingEdges = 0;
	double d2PdLydy = 0;
	double d2PdLydx = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	Vect unitvec;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double directionsin = 0;
	double directioncos = 0;
	int updownfactor = -1;
	int CCWVertFactor = -1;
	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		if ((v1->id) == targetvertexid || (v2->id) == targetvertexid)
		{
			e1 = (*v1).GetEdgeBetween(v2);
			unitvec = (*e1).MakeUnitVector(v2);

			if ((v1->id) == targetvertexid)
			{
				CCWVertFactor = 1;
			}
			else
			{
				CCWVertFactor = -1;
			}
			clockwiseEdges.push_back(e1);
			
			if (((*e1).yflag))
			{
				//directioncos = abs(vec.x); // these absolute values are part of a convention needed to get this derivative right I believe
				//directionsin = abs(vec.y);
				// (07/09/2024 update: I think the abs() needs to be removed, we commit to perimeter being calculated by traversing edges CCW (add a minus sign to everything))
				directioncos = unitvec.x;//-1*vec.x;
				directionsin = unitvec.y;//-1*vec.y;
				if (unitvec.y < 0)
				{
					updownfactor = -1;
				}
				else
				{
					updownfactor = 1;
				}
				// updownfactor is positive if the vector that the edge is on points up (in the plus y hat direction) when we traverse the polygon
				// 
				d2PdLydx +=  LY / (*e1).GetLength() * CCWVertFactor * updownfactor * (1 - directioncos * directioncos); //(1 - directionsin * directionsin ) ;
				d2PdLydy +=  LY / (*e1).GetLength() * CCWVertFactor * updownfactor * (-1 * directioncos * directionsin);
				
				if (ycrossingEdges < 1)
				{
					ycrossingEdges++;
				}
				else
				{
					i = vertexNum-1;
				}

			}
		}
	}
	Vect d2PdLydv(d2PdLydx, d2PdLydy);
	return d2PdLydv;
}
