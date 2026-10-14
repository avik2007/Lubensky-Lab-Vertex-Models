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
#include "AGTissue.h"


using namespace std;

//***************************************************************************************************************
//*********** FUNCTIONS FOR CALCULATING CELL PROPERTIES **********************************************************************
//***************************************************************************************************************
/*
* double GetPressure
*
*
* Calculates pressure of the cell based on the cell potential type (here, we're adding an option for the pressure 
* taking a log or a exponential potential 
*/
double AGCell::GetPressure(double Area)
{
	if ((*this).pressureType == "LOG")
	{
		return K * (1 / (1 + Area / (*this).A0));
		// where the energy would look like Ep = - K A0 Log(1 + A/A0) ([K] = pressure)
	}
	else if ((*this).pressureType == "EXPONENTIAL")
	{
		return K * exp(-Area / (*this).A0);
		// where the energy would look like Ep = K A0 exp(-(A/A0-1)) ([K] = pressure)
		// and the exp(1) in the denominator is put in so that the local concavity of the function
	}
	else if ((*this).pressureType == "QUARTIC")
	{
		return K * ((1 - Area / (*this).A0) + 2 * C4 * pow(1 - Area / (*this).A0, 3));
	}
	else if ((*this).pressureType == "QUARTIC2.0")
		return K * ((1 - Area / (*this).A0) * (1 - 6 * C4) + 2 * C4 * pow(1 - Area / (*this).A0, 3));
	// this is an attempt to get QUARTIC behaving properly. I'm guessing the concavity of the pressure energy 1/2 gets the 
	// vertex model behaving as expected. Which we did not find in the first quartic potential.
	else
	{
		return K * (1 - Area / (*this).A0);
		// where the energy would look like Ep = 1/2 K A0 (A/A0 - 1)^2 ([K] = pressure/L^2)
	}
}

/*
* double GetNegativeDPressureDArea
*
*
* Calculates pressure of the cell based on the cell potential type (here, we're adding an option for the pressure
* taking a log or a exponential potential
*/
double AGCell::GetNegativeDPressureDArea(double Area)
{
	if ((*this).pressureType == "LOG")
	{
		return ( K / (*this).A0 ) * (1 / pow((1 + Area / (*this).A0), 2));
		// where the energy would look like Ep = - K A0 Log(1 + A/A0) and the pressure looks like
		// P = K * (1 / (1 + Area / (*this).A0)) ([K] = pressure)
	}
	else if ((*this).pressureType == "EXPONENTIAL")
	{
		return (K / (*this).A0) * (exp(-Area / (*this).A0));
		// where the energy would look like Ep = K A0 exp(-(A/A0-1)) ([K] = pressure)
		// and the exp(1) in the denominator is put in so that the local concavity of the function
		// matches the others 
		// and pressure looks like P = K * exp(-Area / (*this).A0);
	}
	else if ((*this).pressureType == "QUARTIC")
	{
		return (K / (*this).A0) * ( 1 + 6 * C4 * pow(1 - Area / (*this).A0, 2));
	}
	else if ((*this).pressureType == "QUARTIC2.0")
		return (K / (*this).A0) * ( (1 - 6 * C4) + 6 * C4 * pow(1 - Area / (*this).A0, 2));
	// this is an attempt to get QUARTIC behaving properly. I'm guessing the concavity of the pressure energy 1/2 gets the 
	// vertex model behaving as expected. Which we did not find in the first quartic potential.
	else
	{
		return (K / (*this).A0);
		// where the energy would look like Ep = 1/2 K A0 (A/A0 - 1)^2 ([K] = pressure/L^2)
	}
}
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
* Vect GetdAdv(Vertex* v);
* 
*  Calculates the derivative of a cell's area with respect to one of its vertices
* "Plus" and "Minus" in the variable names refers to counterclockwise and clockwise
*/
Vect AGCell::GetdAdv(Vertex* v)
{
	vector<Vertex*> cwverts;
	int vnum,specm1, specp1;
	specm1 = 0;
	specp1 = 0;
	bool isvertexincell = false;
	Vertex *vplus, *vminus;
	Vect dAdv(0, 0);
	Edge * eplus, *eminus;
	Vect vTovplus, vminusTov;
	cwverts = (*this).SavedInfoForEvolve_CWVerts;
	vnum = cwverts.size();
	
	for (int vindex = 0; vindex < vnum; vindex++)
	{
		if ((*cwverts[vindex]).id == (*v).id)
		{
			if (vindex == 0)
			{
				specm1 = vnum - 1;
				specp1 = 1;
			} 
			else if (vindex == vnum - 1)
			{
				specp1 = 0;
				specm1 = vnum - 2;
			}
			else
			{
				specp1 = vindex + 1;
				specm1 = vindex - 1;
			}
			isvertexincell = ((*cwverts[vindex]).id == (*v).id);
		}
	}
	vplus = cwverts[specp1];
	vminus = cwverts[specm1];
	if (isvertexincell)
	{
		eplus = (*v).GetEdgeBetween(vplus);
		eminus = (*v).GetEdgeBetween(vminus);
		vTovplus = (*eplus).GetVectorDifference(v, vplus);
		vminusTov = (*eminus).GetVectorDifference(vminus, v);
		dAdv =  (vTovplus + vminusTov).Rotate(PI / 2) * (0.5); //NEGATIVE SIGN IS ATTEMPT TO MAINTAIN CCW convention
	}

	return dAdv;
}
/*
* Vect GetdPdv(Vertex* v);
*
*  Calculates the derivative of a cell's perimeter with respect to one of its vertices
*/
Vect AGCell::GetdPdv(Vertex* v)
{
	vector<Vertex*> cwverts;
	int vnum, specm1, specp1;
	specm1 = 0;
	specp1 = 0;
	bool isvertexincell = false;
	Vect dPdv(0, 0);
	Vertex* vplus, * vminus;
	Vect dAdv(0, 0);
	Edge* eplus, * eminus;
	Vect vTovplus, vminusTovplus;
	cwverts = (*this).SavedInfoForEvolve_CWVerts;
	vnum = cwverts.size();

	//for (int vindex = 0; vindex < vnum; vindex++)
	int vindex = 0;

	// Note: this do while loop is one way to get this done. But dAdv() uses a smart for loop which could be used here.
	// 
	do
	{
		isvertexincell = ((*cwverts[vindex]).id == (*v).id);
		if (isvertexincell)
		{
			if (vindex == 0)
			{
				specm1 = vnum - 1;
				specp1 = 1;
			}
			else if (vindex == vnum - 1)
			{
				specp1 = 0;
				specm1 = vnum - 2;
			}
			else
			{
				specp1 = vindex + 1;
				specm1 = vindex - 1;
			}
			
		}
		vindex++;
	} while ((!isvertexincell) && (vindex < vnum));

	vplus = cwverts[specp1];
	vminus = cwverts[specm1];
	if (isvertexincell)
	{
		eplus = (*v).GetEdgeBetween(vplus);
		eminus = (*v).GetEdgeBetween(vminus);
		vTovplus = (*eplus).GetVectorDifference(v, vplus);
		vminusTovplus = (*eminus).GetVectorDifference(vminus, v);
		dPdv = vTovplus / vTovplus.Magnitude() * (-1) + vminusTovplus / vminusTovplus.Magnitude() ;
		// SIGN CONVENTION OF dPdv MEANT TO MAINTAIN COUNTERCLOCKWISE CONVENTION, DESPITE THE FACT THAT I'M USING CW CONVETION IN CODE
	}

	return dPdv;
}
/*
* double GetdAdLx();
*
* Calculates dA/dLx for this cell
* 
* 
*  11/2/2023 - there is a simple way to improve this. in fact, I'm not sure why I wrote it this in way in the first place. But in the first for loop, check the xflag of the edge
*  and then take the average y-coordinate of that edge if the xflag is 1. 
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
	double y1 = 0;
	double y2 = 0;
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
			if (xcrossingEdges < 1)
			{ 
				y1 = 0.5 * (v1->r.y + v2->r.y);
				xcrossingEdges++;
			}
			else
			{
				y2 = 0.5 * (v1->r.y + v2->r.y);
				i = vertexNum;
			}
		
		}
		

	}
	
	return abs(y2 - y1);
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
	double x1 = 0;
	double x2 = 0;
	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVects.push_back((*e1).GetVectorDifference(v1, v2));

		if ( ((*e1).yflag) )
		{
			if (ycrossingEdges < 1)
			{
				x1 = 0.5 * (v1->r.x + v2->r.x);
				ycrossingEdges++;
			}
			else
			{
				x2 = 0.5 * (v1->r.x + v2->r.x);
				i = vertexNum;
			}

		}// if dAdLxdy works, then this should work
		
	}
	/*
	// everything below can be deleted if the commented snippet above works
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
	*/
	return abs(x2 - x1);
}

/*
* double Getd2AdLxdy(Vertex v)
*
* Only returns y component because x component is always zero for a rectangular box
*/
double AGCell::Getd2AdLxdy(Vertex* V)
{
	unsigned long int targetvertexid = (*V).id;
	bool id_match;
	int xcrossingEdges = 0;
	vector<Vect> edgeVects;
	double d2AdLxdy = 0;
	Vertex* v1, * v2;
	Edge* e1, * e2;
	vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	int firstEdgeContribution = 0;
	int secondEdgeContribution = 0;
	int total_contribution = 0;
	double firstEdgey = 0;
	double secondEdgey = 0;
	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		id_match = ((v1->id) == targetvertexid) || ((v2->id) == targetvertexid);
		e1 = (*v1).GetEdgeBetween(v2);
		clockwiseEdges.push_back(e1);
		edgeVects.push_back((*e1).GetVectorDifference(v1, v2));

		if ((*e1).xflag)
		{
			if (xcrossingEdges < 1)
			{
				firstEdgey = (*e1).GetCenter().y;
			}	
			else
			{
				secondEdgey = (*e1).GetCenter().y;
			}
		}
		if (((*e1).xflag) && id_match)
		{
			if (xcrossingEdges < 1)
			{
				firstEdgeContribution = 1;
				xcrossingEdges++;
			}
			else
			{
				secondEdgeContribution = 1;
				i = vertexNum;
			}
		}
	}
	if (firstEdgey < secondEdgey)
	{
		firstEdgeContribution *= -1;
	} 
	else if (firstEdgey > secondEdgey)// take into the account the equal sign
	{
		secondEdgeContribution *= -1;
	}
	return 0.5 * (firstEdgeContribution + secondEdgeContribution);
	
}
/*
* double GetdAdgammaForSimpleShearAlongXaxis()
*
*  Calculates dAdgamma for a simple shear where the deformation occurs along the x-axis (so the unit cell above the main unit cell is shifted to the right)
* 
*  THIS IS INCORRECT - IT DOESN'T HANDLE CORNER EDGES VERY WELL. AND SOMETIMES NOT EVEN EDGE CELLS.
*/
double AGCell::GetdAdgammaForSimpleShearAlongXAxis()
{
	bool id_match;
	int ycrossingEdges = 0;
	double d2AdLxdy = 0;
	Vertex* v1, * v2;
	Edge* e1;
	//vector<Edge*> clockwiseEdges;
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	double vectDiff;
	double firstEdgeContribution = 0;
	double secondEdgeContribution = 0;
	double totalContribution = 0;
	double x1 = 0;
	double x2 = 0;
	Vect ev;
	Vect crosscellvec(0, 0);
	bool evecadder = false;
	double edgeslope;
	double xcw, ycw; // coordinates of vertex 1
	double xint1 = 0; // the x-coordinates of the point at which edges cross the top of the box
	double xint2 = 0; // the x-coordinates of the point at which edges cross the top of the box
	int xcrossings = 0;
	clockwiseVertex.push_back(clockwiseVertex[0]);

	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		e1 = (*v1).GetEdgeBetween(v2);
		ev = (*e1).GetVectorDifference(v1, v2);
		xcw = v1->r.x;
		ycw = v1->r.y;


		// ADJUST BOUNDARY CONDITIONS FOR CROSSING LEFT OR RIGHT BOUNDARY TO GET INTERSECTION RIGHT 
		if (((*e1)).xflag)
		{
			if (ev.x > 0)
			{
				xcw = xcw - LX;
			}
			else
			{
				xcw = LX + xcw;
			}
		}
		// ADJUST BOUNDARY CONDITIONS FOR CROSSING TOP BOUNDARY
		if ((*e1).yflag)
		{
			if (ev.y < 0)
			{
				ycw = LY + ycw;
			}
		}
		// THIS IS THE SECTION WHERE INTERSECTIONS ARE CALCULATED
		if (((*e1)).yflag)
		{
			edgeslope = ev.y / ev.x; // SLOPE OF CROSSING EDGE, TO SOLVE FOR WHERE IT CROSSES Y-BOUNDARY
			if (ycrossingEdges < 1)
			{
				firstEdgeContribution = ycw + 0.5 * ev.y;
				xint1 = xcw + (LY - ycw) / edgeslope; // CALCULATES THE X-COORDINATE OF THE INTERCEPTION OF THE FIRST CROSSING EDGE
			}
			else
			{
				secondEdgeContribution = ycw + 0.5 * ev.y;
				xint2 = xcw + (LY - ycw) / edgeslope; // CALCULATES THE X-COORDINATE OF THE INTERCEPTION OF THE SECOND CROSSING EDGE
			}
			ycrossingEdges++;
		}
		
		if (((*e1)).xflag)
		{
			if (ycrossingEdges == 1)
			{
				if (!((*e1).yflag))
				{
					xcrossings++; // MIGHT BE UNNECESSARY
				}
				else
				{
					if (xint2 > 0)
					{
						xcrossings++; //MIGHT BE UNECESSARY
					}
				}
			}
		}
		
	}
	// xcrossing count might be unnecessary because you got the boundary conditions right before
	if (ycrossingEdges > 0)
	{
		if (xint2 < xint1 && xcrossings == 0)
		{
			totalContribution = secondEdgeContribution - firstEdgeContribution;
		}
		else
		{
			totalContribution = firstEdgeContribution - secondEdgeContribution;
		}
	}
	 return LY * totalContribution; 
}

/*
* static void HelperLEBCUnwrap(cell, cwverts, level, yunwrapped, areasign)
*
* Shared helper for the two Lees-Edwards (simple shear along x) area derivatives below.
* (Added 2026-10-14 to fix sign errors found by a finite-difference check; see
* vm_claude_actions/vm_claude_todo.md, bugs 2 and 3.)
*
* WHAT IT DOES
* Walks the cell's clockwise vertex list v_0, v_1, ..., v_{n-1} and "unwraps" the polygon
* across the periodic boundaries using Edge::GetVectorDifference (which already knows each
* edge's xflag / yflag). For every vertex it records
*   level[k]      = the vertical periodic image the vertex sits in, relative to v_0, i.e. the
*                   running sum of n_y = (number of box heights the edge v_{k-1} -> v_k jumps).
*                   Under Lees-Edwards simple shear gamma, a vertex in image level l is
*                   displaced by l * gamma * LY in x.
*   yunwrapped[k] = the y coordinate of v_k in the unwrapped (continuous) polygon.
*   areasign      = +1 if the unwrapped shoelace sum is positive, -1 if negative (it is
*                   negative for the clockwise list), so that areasign * shoelace = A > 0.
*
* WHY
* With these, dA/dgamma and d2A/(dgamma dy_j) follow exactly from the shoelace formula,
* with no geometric case analysis (corner cells, cells crossing more than once, edge
* direction heuristics), which is where the previous implementations went wrong.
*/
static void HelperLEBCUnwrap(vector<Vertex*>& cwverts, vector<int>& level, vector<double>& yunwrapped, double& areasign)
{
	int n = cwverts.size();
	vector<double> xu(n), yu(n);
	level.assign(n, 0);
	xu[0] = cwverts[0]->r.x;
	yu[0] = cwverts[0]->r.y;
	for (int k = 0; k < n - 1; k++)
	{
		Vertex* va = cwverts[k];
		Vertex* vb = cwverts[k + 1];
		Edge* e = va->GetEdgeBetween(vb);
		Vect d = e->GetVectorDifference(va, vb);                 // r_b - r_a across the boundaries
		int ny = (int)lround((d.y - (vb->r.y - va->r.y)) / LY);  // box heights jumped by this edge
		level[k + 1] = level[k] + ny;
		xu[k + 1] = xu[k] + d.x;
		yu[k + 1] = yu[k] + d.y;
	}
	double shoelace = 0;
	for (int k = 0; k < n; k++)
	{
		int kp = (k + 1) % n;
		shoelace += 0.5 * (xu[k] * yu[kp] - xu[kp] * yu[k]);
	}
	areasign = (shoelace < 0) ? -1.0 : 1.0;
	yunwrapped = yu;
}

/*
* double GetdAdgammaForSimpleShearAlongXAxisALTERNATE()
*
* dA/dgamma for Lees-Edwards simple shear along x: the periodic image above the main box is
* shifted to the right by gamma * LY.
*
* REWRITTEN 2026-10-14 (bug 3 of the FD check). The previous version took |dA/dgamma| from the
* y-midpoints of the two crossing edges, which was right, but chose the sign from the x
* component of a partial edge sum ("crosscellvec"). That heuristic gave the wrong sign for
* about 2.7% of boundary-crossing cells (41 of 1546 in the T2ss tissues), including cells
* that do not touch a corner.
*
* HOW IT WORKS NOW
* In the unwrapped polygon, vertex k is displaced by level[k] * gamma * LY in x (see
* HelperLEBCUnwrap). By the chain rule,
*     dA/dgamma = LY * sum_k level[k] * dA/dx_k,   with   dA/dx_k = areasign * 1/2 (y_{k+1} - y_{k-1})
* (shoelace formula; y in unwrapped coordinates, indices cyclic). Cells that do not cross the
* top/bottom boundary have all levels equal, so the result is exactly 0 for them. The
* function name is kept so that existing callers are unchanged.
*/
double AGCell::GetdAdgammaForSimpleShearAlongXAxisALTERNATE()
{
	vector<Vertex*> cwverts = (*this).GetClockWiseVertices();
	int n = cwverts.size();
	vector<int> level;
	vector<double> yu;
	double areasign;
	HelperLEBCUnwrap(cwverts, level, yu, areasign);
	double dAdgamma = 0;
	for (int k = 0; k < n; k++)
	{
		double dAdxk = areasign * 0.5 * (yu[(k + 1) % n] - yu[(k + n - 1) % n]);
		dAdgamma += level[k] * dAdxk;
	}
	return LY * dAdgamma;
}
/*
* double Getd2AdgammadyForSimpleShearAlongXAxis(Vertex* V)
*
* d2A/(dgamma dy_V): derivative of this cell's dA/dgamma (above) with respect to the y
* coordinate of vertex V. (The x derivative is exactly zero, because dA/dx_k depends only on
* y coordinates, so the shear-modulus code only needs this y component.)
*
* REWRITTEN 2026-10-14 (bug 2 of the FD check). The previous version located the crossing
* points on the top boundary and picked the sign from their x ordering; it had the wrong sign
* for corner cells (cells that also cross the left/right boundary): 122 of 122 failures.
*
* HOW IT WORKS NOW
* Differentiating dA/dgamma = LY * sum_k level[k] * areasign * 1/2 (y_{k+1} - y_{k-1}) with
* respect to y_j (j = position of V in the clockwise list) leaves only the k = j-1 and k = j+1
* terms:
*     d2A/(dgamma dy_j) = areasign * LY/2 * (level[j-1] - level[j+1])
* Returns 0 if V is not a vertex of this cell.
*/
double AGCell::Getd2AdgammadyForSimpleShearAlongXAxis(Vertex* V)
{
	vector<Vertex*> cwverts = (*this).GetClockWiseVertices();
	int n = cwverts.size();
	int j = -1;
	for (int k = 0; k < n; k++)
	{
		if (cwverts[k]->id == V->id) { j = k; }
	}
	if (j < 0) { return 0.0; }
	vector<int> level;
	vector<double> yu;
	double areasign;
	HelperLEBCUnwrap(cwverts, level, yu, areasign);
	return areasign * 0.5 * LY * (level[(j + n - 1) % n] - level[(j + 1) % n]);
}
/*
* double GetCellStressXX();
*
* Calculates the xx component of the Cauchy stress tensor for a cell
*/
double AGCell::GetCellStressXX()
{
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	// clockwiseVertex.push_back(clockwiseVertex[0]);
	Vertex *v1, *v2;
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
		stressxx += ( F1 * x1 ) ;
	}
	return stressxx / (*this).GetArea();
}

/*
* double GetCellStressYY();
*
* Calculates the yy component of the Cauchy stress tensor for a cell
*/
double AGCell::GetCellStressYY()
{
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	Vertex *v1, *v2;
	double F1, F2;
	double y1, y2;
	double stressyy = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		y1 = (*v1).r.y;
		F1 = (*v1).SavedInfoForEvolve_TotalForce.y;
		stressyy += (F1 * y1) ;
	}
	return stressyy / (*this).GetArea();
}
/*
* double GetCellStressInterpolateXX();
*
* Calculates the xx component of the Cauchy stress tensor for a cell using the 
* interpolating equation defined in Lee and Morishita 2017 
*
* THIS IS INHERENTLY DIFFERENT THAT LEE AND MORISHITA'S STRESS TENSOR BECAUSE
* YOU OVERCOUNT FORCES - LEE AND MORISHITA THINK ONLY ABOUT FORCES FROM AN INDIVIDUAL
* CELL, MEANWHILE, YOU COUNTED ALL OF THEM
*/
double AGCell::GetCellStressInterpolateXX()
{
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	clockwiseVertex.push_back(clockwiseVertex[0]);
	Vertex *v1, *v2;
	double F1, F2;
	double x1, x2;
	double stressxx = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		x1 = (*v1).r.x;
		x2 = (*v2).r.x;
		F1 = (*v1).SavedInfoForEvolve_TotalForce.x;
		F2 = (*v1).SavedInfoForEvolve_TotalForce.x;
		stressxx += ( F1 * x1 + F2 * x2 ) / ( 3.0 ) + (F1 * x2 + F2 * x1 )  / ( 6.0 ) ;
	}
	return stressxx / (*this).GetArea();
}

/*
* double GetCellStressInterpolateYY();
*
* Calculates the xx component of the Cauchy stress tensor for a cell using the
* interpolating equation defined in Lee and Morishita 2017
*
* THIS IS INHERENTLY DIFFERENT THAT LEE AND MORISHITA'S STRESS TENSOR BECAUSE
* YOU OVERCOUNT FORCES - LEE AND MORISHITA THINK ONLY ABOUT FORCES FROM AN INDIVIDUAL
* CELL, MEANWHILE, YOU COUNTED ALL OF THEM
*/
double AGCell::GetCellStressInterpolateYY()
{
	vector<Vertex*> clockwiseVertex = (*this).GetClockWiseVertices();
	int vertexNum = clockwiseVertex.size();
	clockwiseVertex.push_back(clockwiseVertex[0]);
	Vertex *v1, *v2;
	double F1, F2;
	double y1, y2;
	double stressyy = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		v1 = clockwiseVertex[i];
		v2 = clockwiseVertex[i + 1];
		y1 = (*v1).r.x;
		y2 = (*v2).r.x;
		F1 = (*v1).SavedInfoForEvolve_TotalForce.y;
		F2 = (*v1).SavedInfoForEvolve_TotalForce.y;
		stressyy += (F1 * y1 + F2 * y2) / (3.0) + (F1 * y2 + F2 * y1) / (6.0);
	}
	return stressyy / (*this).GetArea();
}

double AGCell::GetEnergy()
{
	double energy = 0.5 * (K * pow((*this).GetArea() - (*this).A0, 2) / (*this).A0) + 0.5*GAMMA * (*this).GetPerimeter();
	return energy;

}