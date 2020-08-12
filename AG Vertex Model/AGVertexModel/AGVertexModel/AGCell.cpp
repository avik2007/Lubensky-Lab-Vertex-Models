/**	\file	AGCell.cpp
*	Implementation of AGCell (anisotropically growing tissue cell)	class.
*	Contains function definitions associated with AGCell objects.
*/

#include <iostream>
#include <cassert>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <algorithm>
#include <stddef.h>

#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "VectClass.h"
#include "Constants.h"
#include "AGConstants.h"
#include "AGCell.h"


using namespace std;

//***************************************************************************************************************
//*********** FUNCTIONS FOR CALCULATING CELL PROPERTIES **********************************************************************
//***************************************************************************************************************

/*
 * double GetThetaTCJ
 *
 * Calculates ThetaTCJ (angle that determines division axes)
 * see Bosveld, F.et al Y. (2016). Epithelial tricellular junctions act as interphase 
 * cell shape sensors to orient mitosis. Nature, 530(7591), 495–498. https://doi.org/10.1038/nature16970
 */
double AGCell::GetThetaTCJ()
{
	double ThetaTCJ;
	double vlongx; // x component of scaled eigenvector of bipolarity tensor associated with the cell's long axis 
	// (scaled such that the y component is 1)
	double eplus, eminus; //eigenvalues of bipolarity tensor
	// eplus is supposed to be the bigger of the two at the end of the function
	vector<vector<double>> bipolarityTensor = (*this).GetBipolarityTensor();
	double s11 = bipolarityTensor[0][0];
	double s12 = bipolarityTensor[0][1];
	double s21 = bipolarityTensor[1][0];
	double s22 = bipolarityTensor[1][1];
	//tensor elements 

	double trace = s11 + s22;
	//trace of bipolarity tensor
	double det = s11*s22 - s12*s21;
	//determinant of bipolarity tensor

	eplus = (trace / 2.0) + sqrt(pow(trace / 2.0, 2) - det);
	eminus = (trace / 2.0) - sqrt(pow(trace / 2.0, 2) - det);

	if (eplus > eminus)
	{
		vlongx = eplus - s22;
	}
	else 
	{
		vlongx = eminus - s22;
	}
	
	ThetaTCJ = atan(1.0 / vlongx);
	return ThetaTCJ;
}
/*
* double GetTCJAnisotropy
*
* Calculates cellular anisotropy, as determined by the TCJ distribution (angle that determines division axes)
* Helper function for GetThetaTCJ()
* see Bosveld, F.et al Y. (2016). Epithelial tricellular junctions act as interphase
* cell shape sensors to orient mitosis. Nature, 530(7591), 495–498. https://doi.org/10.1038/nature16970
*/
double AGCell::GetTCJAnisotropy()
{
	double TCJAnisotropy;
	double eplus, eminus; //eigenvalues of bipolarity tensor
	// eplus is supposed to be the bigger of the two at the end of the function
	vector<vector<double>> bipolarityTensor = (*this).GetBipolarityTensor();
	double s11 = bipolarityTensor[0][0];
	double s12 = bipolarityTensor[0][1];
	double s21 = bipolarityTensor[1][0]; 
	double s22 = bipolarityTensor[1][1];
	//tensor elements 

	double trace = s11 + s22;
	//trace of bipolarity tensor
	double det = s11*s22 - s12*s21;
	//determinant of bipolarity tensor

	eplus = (trace / 2.0) + sqrt(pow(trace / 2.0, 2) - det);
	eminus = (trace / 2.0) - sqrt(pow(trace / 2.0, 2) - det);

	// if eplus is the smaller of the two eigenvalues, swap
	if (eplus < eminus)
	{
		eplus = eplus + eminus;
		eminus = eplus - eminus;
		eplus = eplus - eminus;
	}

	TCJAnisotropy = 1 - eminus / eplus;

	return TCJAnisotropy;
}
/*
* list<Vect> GetBarycentriUnitVectors
*
* Calculates unit vectors that point from cell center to TCJ (vertices)
* This is a helper function for GetThetaTCJ
* see Bosveld, F.et al Y. (2016). Epithelial tricellular junctions act as interphase
* cell shape sensors to orient mitosis. Nature, 530(7591), 495–498. https://doi.org/10.1038/nature16970
*/
list<Vect> AGCell::GetBarycentricUnitVectors()
{
	int i;
	Vect GeoCent, tempvec;
	vector<vector<double>> vertposforflag;
	list<Vect> unitvectors; // CHECK IF THIS CAN STAY A LIST, IT MAY NEED TO BEOCOME A VECTOR!
	double NTCJ, tempveclength;

	// Calculate center of cell
	vertposforflag = (*this).HelperOrderVerticesPos(); //allows positions of vertices outside of domain
	// Undoes PBC, and gives coordinates of TCJ's (vertices)
	GeoCent = GetGeoCenter(vertposforflag[0], vertposforflag[1]);
	// Calculates the center of the vertices 

	NTCJ = (*this).GetVertexNumber();
	

	// moves points from arrays of arrays of doubles to an array of vects
	for (i = 0; i < NTCJ; i++)
	{
		tempvec.x = vertposforflag[0][i] - GeoCent.x;
		tempvec.y = vertposforflag[1][i] - GeoCent.y;
		tempveclength = tempvec.Magnitude();
		tempvec.x = tempvec.x / tempveclength; // normalize, to length 1
		tempvec.y = tempvec.y / tempveclength;

		unitvectors.push_back(tempvec);
	}
	return unitvectors;
}


/*
* vector<vector<double>> GetBipolarityTensor()
*
* Calculates bipolarity tensor, whose eigenvectors orient the cell's TCJ distribution
* This is a helper function for GetThetaTCJ
*
* I'm using an array of arrays 
* see Bosveld, F.et al Y. (2016). Epithelial tricellular junctions act as interphase
* cell shape sensors to orient mitosis. Nature, 530(7591), 495–498. https://doi.org/10.1038/nature16970
*/
vector<vector<double>> AGCell::GetBipolarityTensor()
{
	int xd, yd;
	list<Vect> barycentricUnitVectors = (*this).GetBarycentricUnitVectors();
	vector<vector<double>> bipolarityTensor{ { 0, 0 }, { 0, 0 } };
	vector<vector<double>> tempVector{ { 0, 0 }, { 0, 0 } };
	int Ndim = 2; //This function really can only work for the 2 dimensional case

	for (list<Vect>::iterator bit = barycentricUnitVectors.begin(); bit != barycentricUnitVectors.end(); bit++)
	{
		//BREAKPOINT ERROR HERE????
		// cout << (*bit).x << " " << (*bit).y << "\n"; // debug line
		tempVector[0][0] = (*bit).x * (*bit).x;
		tempVector[0][1] = (*bit).x * (*bit).y;
		tempVector[1][0] = (*bit).y * (*bit).x;
		tempVector[1][1] = (*bit).y * (*bit).y;
		for (xd = 0; xd < Ndim; xd++)
		{
			for (yd = 0; yd < Ndim; yd++)
			{
				bipolarityTensor[xd][yd] += tempVector[xd][yd];
			}
		}
	}

	return bipolarityTensor;
}

//************************************************************************************************************************
//********************** Division Angle for Mitosis Functions ***********************************************************
//***********************************************************************************************************************

/*
* double TCJDivPotential(double strength, double phi)
*
* Returns the probability of a certain division angle using Boltzman distribution with the potential
* developed in the Bosveld et al. paper, I'm using the first harmonic ( cos(2(phi))) where 
* phi = theta - thetatcj
*
* Possible change: remove thetatcj as a parameter and use (*this).GetThetaTCJ
*
* 
* see Bosveld, F.et al Y. (2016). Epithelial tricellular junctions act as interphase
* cell shape sensors to orient mitosis. Nature, 530(7591), 495–498. https://doi.org/10.1038/nature16970
*/
double AGCell::TCJDivProbability(double phi)
{
	return exp(-1 * (*this).TCJSTRENGTH * cos(2 * phi));
	// double cosphi = 100 * cos(phi) * cos(phi);
	// return exp(-1 * cosphi) / (exp(-100 / 2.0) * PI * b.bessi(0, 100 / 2.0));
	// for debugging, we add a "strength" of 10
}
/*
* double AnisotropicDivPotential(double strength, double phi)
*
* Returns the probability of a certain division angle using a probability distribution inspired by 
* magnetic anisotropic materials (U = k Cos^2(phi), with phi = theta - theta_an. I need a function
* that effectively parameterizes the directions of anisotropy as well as the strength of the anisotropy
*/

double AGCell::AnisotropicDivProbability( double phi)
{
	double cosphi = ANSTRENGTH * cos(phi) * cos(phi);
	return exp( -1 * cosphi );
}

/*
* double CombTCJAnProbability( double phi)
*
* Returns the probability of a certain division angle using a probability distribution that combines both
* the MudTCJ potential and the one used for the AnisotropicDivPotential. The ultimate form of the
* potential is a bit messy because U = Cos(2 ( theta - thetaTCJ)) + k Cos^2(theta - thetaan). Normalizing 
* this probability distribution is non-trivial.
*/
double AGCell::CombTCJANProbability( double phi)
{
	double ThetaTCJ = (*this).GetThetaTCJ();
	double cosConst = (TCJSTRENGTH) * cos(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * cos(2 * THETAAN);
	double sinConst = (TCJSTRENGTH) * sin(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * sin(2 * THETAAN);
	double U = cosConst * cos(2 * phi) + sinConst * sin(2 * phi);
	return exp(-1 * U) ;

}
/*
* double CombTCJAnProbability( double phi, double thetatcj)
*
* Returns the probability of a certain division angle using a probability distribution that combines both
* the MudTCJ potential and the one used for the AnisotropicDivPotential. The ultimate form of the
* potential is a bit messy because U = Cos(2 ( theta - thetaTCJ)) + k Cos^2(theta - thetaan). Normalizing
* this probability distribution is non-trivial. Adds thetatcj as a parameter to speed up DivisionAngleDistribution
*/
double AGCell::CombTCJANProbability(double phi, double thetatcj)
{
	// double ThetaTCJ = (*this).GetThetaTCJ();
	double ThetaTCJ = thetatcj;
	double cosConst = cos(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * cos(2 * THETAAN);
	double sinConst = sin(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * sin(2 * THETAAN);
	double U = cosConst * cos(2 * phi) + sinConst * sin(2 * phi);
	return exp(-1 * U) ;

}

/*
* double CombTCJANProbMaximum()
*
* Returns the maximum of the combined (tcj + anisotropic) division potential
*/
double AGCell::CombTCJANProbMaximum()
{
	double ThetaTCJ = (*this).GetThetaTCJ();
	double cosConst = (TCJSTRENGTH)* cos(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * cos(2 * THETAAN);
	double sinConst = (TCJSTRENGTH)* sin(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * sin(2 * THETAAN);
	double Umax = sqrt(pow(cosConst, 2) + pow(sinConst, 2));
	// maximum of potential 
	return exp(-1 * Umax) ;
}
/*
* double CombTCJANProbMaximum()
*
* Returns the maximum of the combined (tcj + anisotropic) division potential
*/
double AGCell::CombTCJANProbMaximum(double thetatcj)
{
	// double ThetaTCJ = (*this).GetThetaTCJ();
	double ThetaTCJ = thetatcj;
	double cosConst = cos(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * cos(2 * THETAAN);
	double sinConst = sin(2 * ThetaTCJ) + (ANSTRENGTH / 2.0) * sin(2 * THETAAN);
	double Umax = sqrt(pow(cosConst,2) + pow(sinConst,2));
	return exp(-1 * Umax);
}
/*
* double DivisionAngleDistribution(string potential)
*
* Generates a random angle that will be used by Mitosis function to help determine the cleavage plane of 
* dividing cells. Uses rejection sampling algorithm 
*/
double AGCell::DivisionAngleDistribution(string potential)
{
	// in-code labels for the different potentials available currently
	// to add a new one, add a new string declaration 
	string tcjpotential("tcj");
	string anisotropic("ani");
	string combine("combined");


	// random value generators for the rejection sampler
	uniform_real_distribution<double> ydist(0.0, 1.0);
	uniform_real_distribution<double> phidist(0.0, PI);
	// proposal values for rejection sampling
	double propy, propphi, propPofPhi; 
	// boolean value that will end the rejection sampling loop
	bool flag = true;
	// final angle that will be returned
	double angle; 


    double ThetaTCJ = (*this).GetThetaTCJ();
	// double ThetaTCJ = 0;

	// rejection sampler code
	do
	{
		propphi = phidist(generator);
		propy = ydist(generator);
		//selecting for the different probability distributions

		// first, tcj potential
		if (potential.compare(tcjpotential) == 0 )
		{
			propPofPhi = TCJDivProbability(propphi);
			// to properly do rejection sampling, we have to sample such that propy ranges from 0 to the max of PofPhi
			// PofPhi's max changes depending on the properties of the distribution, so we need to scale propy accordingly
			propy = TCJDivProbability(PI / 2.0) * propy;
		}
		// anisotropic division potential
		else if ( potential.compare(anisotropic) == 0 )
		{
			propPofPhi = AnisotropicDivProbability(propphi);
			// to properly do rejection sampling, we have to sample such that propy ranges from 0 to the max of PofPhi
			// PofPhi's max changes depending on the properties of the distribution, so we need to scale propy accordingly
			propy = AnisotropicDivProbability(PI / 2.0) * propy;
		}
		// combine anisotropic / tcj potential
		else if (potential.compare(combine) == 0)
		{
			propPofPhi = CombTCJANProbability(propphi, ThetaTCJ);
			angle = propphi;
			propy = propy * CombTCJANProbMaximum(ThetaTCJ);
			// remember, for rejection sampling, you have to have your "y-proposal" variable be a uniform random variable from 
			// 0 to the maximum of the desired potential. For this combined potential, the location of the max can change as a 
			// function of the ThetaTCJ and ThetaAn, so it was worth making a new function

		}
		else
		{
			// if no valid potential is given by user, default to isotropic divisions
			// might be worth changing default to TCJ
			propPofPhi = propphi;
			propy = propphi;
			// by setting these equal, we gurantee that the rejection sample loop will be exited
		}

	
		if ( propy <= propPofPhi ) 
		{
			flag = false;
			// if we get an acceptable value, end the loop

			if (potential.compare(tcjpotential) == 0)
			{
				angle = propphi + ThetaTCJ;
			}
			else if (potential.compare(anisotropic) == 0)
			{
				angle = propphi + THETAAN;
			}
			else
			{
				angle = propphi;
			}
			// make the appropriate angle declaration
		}
	} while (flag);

	return angle;
}
/*
* double DivisionAngleDistribution()
*
* Generates a random angle that will be used by Mitosis function to help determine the cleavage plane of
* dividing cells. Same as DivisionAngleDistribution(potential), but will use the field "DIVISIONPOTENTIAL"
*/
double AGCell::DivisionAngleDistribution()
{
	double theta = (*this).DivisionAngleDistribution(DIVISIONPOTENTIAL);
	return theta;
}
/*
* bool IsCrossingX()
*
* Checks whether a cell crosses the X-boundary
*/
bool AGCell::IsCrossingX()
{
	bool result = false; // true is the edge crosses X-boundary, false otherwise
	for (vector<Edge*>::iterator eit = (*this).edges.begin(); eit != (*this).edges.end(); eit++)
	{
		result = ( (*(*eit)).xflag || result );
	}
	return result;
}
/*
* bool IsCrossingY()
*
* Checks whether a cell crosses the Y-boundary
*/
bool AGCell::IsCrossingY()
{
	bool result = false; // true is the edge crosses X-boundary, false otherwise
	for (vector<Edge*>::iterator eit = (*this).edges.begin(); eit != (*this).edges.end(); eit++)
	{
		result = ( (*(*eit)).yflag || result ); 
	}
	return result; 
}
/*
* double GetdAdLx();
*
* Calculates dA/dLx for this cell
*/
double AGCell::GetdAdLx()
{
	int xcrossingEdges = 0;
	vector<Vect> edgeVects;
	double dAdLx = 0;
	Vertex *v1, *v2;
	Edge *e1, *e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;

	clockwiseVertex.push_back(clockwiseVertex[0]);
	
	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVects.push_back((*e1).GetVectorDifference(v1, v2));
	}

	for (int j = 0; j < vertexNum; j++)
	{
		e2 = clockwiseEdges[j];
		vectDiff = (edgeVects[j]).y; 
		if (xcrossingEdges == 0)
		{
			if ((*e2).xflag)
			{
				dAdLx += 0.5 * vectDiff;
				xcrossingEdges++;
			}
		}
		else if (xcrossingEdges == 1)
		{
			if ((*e2).xflag)
			{
				dAdLx += 0.5 * vectDiff;
				xcrossingEdges++;
			}
			else
			{
				dAdLx += vectDiff;
			}
		}
		else if (xcrossingEdges == 2)
		{
			j = vertexNum;
		}
	}
	return abs(dAdLx);
}
/*
* double GetdAdLy();
*
* Calculates dA/dLy for this cell
*/
double AGCell::GetdAdLy()
{
	int ycrossingEdges = 0;
	vector<Vect> edgeVects;
	double dAdLy = 0;
	Vertex *v1, *v2;
	Edge *e1, *e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;

	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVects.push_back((*e1).GetVectorDifference(v1, v2));
	}

	for (int j = 0; j < vertexNum; j++)
	{
		e2 = clockwiseEdges[j];
		vectDiff = (edgeVects[j]).x;
		if (ycrossingEdges == 0)
		{
			if ((*e2).yflag)
			{
				dAdLy += 0.5 * vectDiff;
				ycrossingEdges++;
			}
		}
		else if (ycrossingEdges == 1)
		{
			if ((*e2).yflag)
			{
				dAdLy += 0.5 * vectDiff;
				ycrossingEdges++;
			}
			else
			{
				dAdLy += vectDiff;
			}
		}
		else if (ycrossingEdges == 2)
		{
			j = vertexNum;
		}
	}
	return abs(dAdLy);
}

