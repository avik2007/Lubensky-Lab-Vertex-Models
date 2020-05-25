//! \file MegaCellConstants.h
//! This file contains declaration of global variables/parameters for the simulation.
//! The following variables are global and thus are accessible from any function in the simulation.
//! Note that this file only contains declaration of the global variables. Their values should be set in the Main.cpp file or through file/user input.
//! This only includes the variables needed for MegaCell dynamics (consult Constants.h for all other constants)

// include guard
#ifndef __MEGACELLCONSTANTS_H_INCLUDED__
#define __MEGACELLCONSTANTS_H_INCLUDED__

#include "VectClass.h"

using namespace std;

extern double factor_increased_tension_at_inner_boundary;	//!< Allow for increased tension at inner boundary (where PCP will be initialized)
extern double inner_pressure;								//!< Allow for a pressure outward from the inner cell.
extern double outer_pressure;								//!< Allow for a pressure inward from the outer cell.
extern Vect origin_mega_cell;								//!< Store origin of mega cell for calculating angle-dependent tension (SIMPLEPCPMEGA)/ initializing line of differentiated cells (SIMPLEPCPMEGA/PCPMEGA)
extern int is_voronoi_packing;								//!< If the packing is a voronoi packing, pressure of the central cell is evaluated differently.
extern int center_cell_type;								//!< Type of the central mega cell.
extern int outer_cell_type;									//!< Type of outer mega cell.


#endif