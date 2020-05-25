// include guard
#ifndef __FIBERCONSTANTS_H_INCLUDED__
#define __FIBERCONSTANTS_H_INCLUDED__

#include <string>
#include <random>

using namespace std;

extern double LAMBDA;						//!< LAMBDA=SUM(tansion*length) of each edge and fiber in the cell
extern double TAU;					//controls the rate of myosin accumulation tau*gammadot=f(x)
extern double A0SKEW;
extern double SPRINGCONST;			//for use with fiberedge_spring gives k in F=-kx
extern double DTENSIONCONST;
extern double ANGLECONST;	//tenson goes like GAMMA(1-ANGLECONST*cos(1*theta))
extern bool MAKEFIBERS;				//If true then the tissue can make fibers, otherwise their formation is inhibited. 
extern string MYOSINDYNAMICS; 
extern double FIB_BIRTHRATE;
extern double FIB_DEATHRATE;
extern double FIB_INITPROB;
extern double STRAINDIST; 
extern double LX0; 
extern double LY0; 
extern bool StrainY; 
extern double RELAXINTERVAL; 
extern double BETAZ;
extern double BETAGAMMA;
extern double B;
extern double T1PERCENT;			//probability of t1 per edge per t0 time unit
extern double L0; 
extern int MITOSISANGLE;   //for oriented devisions 0=random angle 1=horizontal 2=vertical
extern double PULLINGFORCETOP; 
extern double PULLINGFORCEBOTTOM;  
extern bool USESPRINGDYNAMICS; //if true use spring dynamics if false use gamma dot=escess force 
extern int TOTALT1; 
extern bool CFO; 
extern double MAXGAMMA; 
extern double MINPERPFORCE;
extern double distance_to_next_vert_cutoff; //added by Hayden
extern double starting_distance_from_first_vert; //added by Hayden
extern double fiber_nuc_cutoff_for_death; //added by Hayden
extern double fiber_nuc_cutoff_for_nucleation; //added by Hayden

#endif