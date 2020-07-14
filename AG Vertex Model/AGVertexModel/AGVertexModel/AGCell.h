// define guard
#ifndef __AGCELL_H_INCLUDED__
#define __AGCELL_H_INCLUDED__

// forward class dependencies

class Cell; 
class Vertex;
class Vect;
class Edge;

//includes 
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <random>

#include "CellClass.h"


using namespace std;
// AGCell = developed for anisotropically growing models.
class AGCell : public Cell{
public:

//***************************************************************************************************************
//********** Additional Data ************************************************************************************
//***************************************************************************************************************
	//double ThetaTCJ; 
	//double TCJAnisotropy;
	double TCJSTRENGTH;
	double THETAAN;
	double ANSTRENGTH;
	string DIVISIONPOTENTIAL;
	bool xflag;
	bool yflag;
//***************************************************************************************************************
//*********** CONSTRUCTORS AND DESTRUCTORS **********************************************************************
//***************************************************************************************************************
	AGCell(){};
	virtual ~AGCell(){};
	//destructor: does not need to be called by you unless you created the object using the word "new" otherwise the compiler will take care of this for you
	AGCell(int idnumb, vector<Edge*> edglist, int typeinput, double A0input){
		id = idnumb;
		edges = edglist;
		type = typeinput;
		A0 = A0input;
	}
	
	AGCell(int idnumb, vector<Edge*> edglist, int typeinput, double A0input, double tcjstrength, double thetaan, double anstrength, string divp){
		id = idnumb;
		edges = edglist;
		type = typeinput;
		A0 = A0input;
		TCJSTRENGTH = tcjstrength;
		THETAAN = thetaan;
		ANSTRENGTH = anstrength;
		DIVISIONPOTENTIAL = divp;
	}
	

//***************************************************************************************************************
//*********** FUNCTIONS FOR CALCULATING CELL PROPERTIES **********************************************************************
//***************************************************************************************************************

	double GetThetaTCJ(); // helps determine orientation of cell division using TCJ distribution

	double GetTCJAnisotropy(); // measure of anisotropy in cells 

	list<Vect> GetBarycentricUnitVectors(); //calculates unit vectors between TCJ's and cell center
	// this is a helper for getThetaTC() and GetBiopolarityTensor()

	vector<vector<double>> GetBipolarityTensor(); // uses barycentric unit vectors to calculate bipolarity tensor
	// this is a helper for GetThetaTCJ()

	bool IsCrossingX(); // checks whether a cell is crossing the X-boundary

	bool IsCrossingY(); // checks whether a cell is crossing the Y-boundary

	double GetdAdLx(); // gets the derivative of the cell's area w.r.t LX

	double GetdAdLy(); // gets the derivative of the cell's area w.r.t. LY

//************************************************************************************************************************
//********************** Division Angle for Mitosis Functions ***********************************************************
//***********************************************************************************************************************
	double TCJDivProbability(double phi);
	// outputs the probability of phi being selected as a cleavage angle, using a Boltzmann distribution with
	// the MudTCJ potential (see Bosveld et al, 2016)

	double AnisotropicDivProbability(double phi);
	// outputs the probability of a certain angle phi being selected as cleavage angle, taking "strength" as a 
	// parameter for a potential inspired by that of magnetic materials (  U = STRENGTH * cos^2(phi)  )

	double CombTCJANProbability(double phi);
	// outputs the probability of a certain angle phi being selected as cleavage angle, taking a strength
	// parameter to indicate strength of anisotropy. The potential for this is the combination of the 
	// MudTCJ potential and the anisotropic material potential used in AnisotropicDivProbability

	double CombTCJANProbability(double phi, double thetatcj);
	// outputs the probability of a certain angle phi being selected as cleavage angle, taking a strength
	// parameter to indicate strength of anisotropy. The potential for this is the combination of the 
	// MudTCJ potential and the anisotropic material potential used in AnisotropicDivProbability
	// adding thetatcj as an extra parameter for convenience

	double CombTCJANProbMaximum();
	// outputs the maximum of the probability distribution of the combined TCJ/anisotropic division potentials
	// this is necessary for the rejection sampling algorithm (in DivisionAngleDistribution) so that we can
	// properly scale proposal random variable 

	double CombTCJANProbMaximum(double thetatcj);
	// outputs the maximum of the probability distribution of the combined TCJ/anisotropic division potentials
	// this is necessary for the rejection sampling algorithm (in DivisionAngleDistribution) so that we can
	// properly scale proposal random variable 
	// we add thetatcj as a parameter for making DivisionAngleDistribution()

	virtual double DivisionAngleDistribution(string potential);
	// uses rejection sampling to produce a random number distributed according to a Boltzmann distribution,
	// using a specific potential as an energy. 

	double DivisionAngleDistribution();
	// uses rejection sampling to produce a random number distributed according to a Boltzmann distribution,
	// using a specific potential as an energy. Same as DivisionAngleDistribution(string potential), just uses
	// AGCell's field "potential"


};
#endif
