//! \file AGStiffConstants.h
//! This file contains declaration of global variables/parameters for the anisotropic growth vertex model
//!  extension, AGStiffness
//! The following variables are global and thus are accessible from any function in the simulation.
//! Note that this file only contains declaration of the global variables. Their values should be set in the Main.cpp file or through file/user input.

// include guard
#ifndef __AGSTIFFCONSTANTS_H_INCLUDED__
#define __AGSTIFFCONSTANTS_H_INCLUDED__

#include <string>
#include <random>

using namespace std;

extern double Kstdratio;
extern double A0stdratio;
extern double LAMBDAstdratio;
#endif