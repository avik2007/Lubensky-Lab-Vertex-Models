//! \file AGPerimConstants.h
//! This file contains declaration of global variables/parameters for the anisotropic growth vertex model
//!  extension, AGPerimVertex.
//! The following variables are global and thus are accessible from any function in the simulation.
//! Note that this file only contains declaration of the global variables. Their values should be set in the Main.cpp file or through file/user input.

// include guard
#ifndef __AGPERIMCONSTANTS_H_INCLUDED__
#define __AGPERIMCONSTANTS_H_INCLUDED__

#include <string>
#include <random>

using namespace std;

extern double LAMBDA;
extern double P0;
extern double DEFORM_NUM;
extern int energyListSize;
extern double equilibriumtolerance;
#endif

