#include <string>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include<cassert> 
#include <algorithm> 
#include <stddef.h>

#include "EdgeClass.h"
#include "CellClass.h"
#include "VertexClass.h"
#include "NucleationSite.h"
#include "FiberEdgeClass.h"
#include "FiberCellClass.h"
#include "FiberConstants.h"

using namespace std;

//Finds the vector that bisects the angle formed by the edges
Vect NucleationSite::GetOutwardDirection(){
	Vect r;

	//Get the angles of the two edges 
	double a1 = (*edge1).GetAngle(parentvertex);
	double a2 = (*edge2).GetAngle(parentvertex);

	//Create two unit vectors in the direction of the edges 
	Vect unit1 = Vect(a1);
	Vect unit2 = Vect(a2);

	r.x = (unit1.x + unit2.x) / 2;
	r.y = (unit1.y + unit2.y) / 2;

	return r;
}

//Finds the vector that bisects the angle formed by the edges
double NucleationSite::GetOutwardAngle(){
	Vect r = (*this).GetOutwardDirection();
	//return r.Angle();
	return atan(r.y / r.x);
}

// cos 2n of outward angle is used to characterize fiber nucleation
double NucleationSite::Calculate_Cos2OutwardAngle(){
	double calc_angle = (*this).GetOutwardAngle();

	return cos(2 * calc_angle);
}

//Finds the angle between edge 1 and edge 2
double NucleationSite::GetAngle(){
	//Get the angles of the two edges 
	double a1 = (*edge1).GetAngle(parentvertex);
	double a2 = (*edge2).GetAngle(parentvertex);
	
	Vect u1 = Vect(a1); 
	Vect u2 = Vect(a2); 

	double ctheta = u1.x*u2.x + u1.y*u2.y; 

	return acos(ctheta); 
}

//Returns the probability of birthing a fiber in 1 min 
double NucleationSite::ProbFiberBirth(){
	double p = 0; 
	double cos2n_calculated = (*this).Calculate_Cos2OutwardAngle();
	double outward_angle_calc = (*this).GetOutwardAngle();
	outward_angle_calc = fabs(outward_angle_calc);
	double opening_angle_calc = (*this).GetAngle();

	if (outward_angle_calc <= fiber_nuc_cutoff_for_nucleation){
		p = NUCLEATIONRATE;
	}

	return p; 
}

//Returns the probability of fiber death (added by Hayden)
double NucleationSite::ProbFiberDeath(){
	double p = 0;
	double cos2n_calculated = (*this).Calculate_Cos2OutwardAngle();
	double outward_angle_calc = (*this).GetOutwardAngle();
	outward_angle_calc = fabs(outward_angle_calc);
	double opening_angle_calc = (*this).GetAngle();

	if (outward_angle_calc >= fiber_nuc_cutoff_for_death){
		p = PASSEDDEATHRATE;
	}

	return p;
}

//double NucleationSite::Calc_XDist_toDeathVertex(){
//	if ((*this).ProbFiberBirth() > 0){
//		vector<NucleationSite*> this_fiber_nuc_site_vect = (*cell).fibernucsites;
//		vector<double> corresponding_x_positions;
//
//		for (int nuc_site_index = 0; nuc_site_index != this_fiber_nuc_site_vect.size(); nuc_site_index++){
//			NucleationSite* this_fiber_nuc_site = this_fiber_nuc_site_vect[nuc_site_index];
//			Vertex* corresponding_vertex_ptr = (*this_fiber_nuc_site).parentvertex;
//			if ((*this_fiber_nuc_site).ProbFiberDeath() > 0){
//				corresponding_x_positions.push_back((*corresponding_vertex_ptr).r.x);
//			}
//		}
//		sort(corresponding_x_positions.begin(), corresponding_x_positions.end()); //sorts the vector min to max
//
//		// Get the list of ordered vertex positions
//		vector<vector<double>> xypos = (*cell).HelperOrderVerticesPos();
//		vector<double> xpos = xypos[0];
//		vector<double> ypos = xypos[1];
//		Vect geocent_this_cell = GetGeoCenter(xpos, ypos);
//
//		// calculate current position of vertex (parent)
//		Vect pos_nucleating_site = (*parentvertex).r;
//
//		// find next vertex that the fiber will run into along this edge
//		Vertex* nextv1 = (*edge1).GetOtherVertex(parentvertex);
//		Vertex* nextv2 = (*edge2).GetOtherVertex(parentvertex);
//		// store corresponding fiber nucleation sites
//		NucleationSite* nextv1_nucsite = NULL;
//		NucleationSite* nextv2_nucsite = NULL;
//		for (int nuc_site_index = 0; nuc_site_index != this_fiber_nuc_site_vect.size(); nuc_site_index++){
//			NucleationSite* this_fiber_nuc_site = this_fiber_nuc_site_vect[nuc_site_index];
//			Vertex* corresponding_vertex_ptr = (*this_fiber_nuc_site).parentvertex;
//
//			if ((*corresponding_vertex_ptr).id == (*nextv1).id){
//				nextv1_nucsite = this_fiber_nuc_site;
//			}
//
//			if ((*corresponding_vertex_ptr).id == (*nextv2).id){
//				nextv2_nucsite = this_fiber_nuc_site;
//			}
//		}
//
//		// find the positions of next vertices
//		Vect nextv1_pos = (*nextv1).r;
//		Vect nextv2_pos = (*nextv2).r;
//
//		// calculate distance along x to each of these vertices
//		double dist_to_v1 = (nextv1_pos - pos_nucleating_site).x;
//		double dist_to_v2 = (nextv2_pos - pos_nucleating_site).x;
//
//		// make this x distance signed (positive if death vertex is on the side closest to geometric center of cell).
//		if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v1 <= 0){
//			dist_to_v1 = -1 * fabs(dist_to_v1);
//		}
//		else{
//			dist_to_v1 = fabs(dist_to_v1);
//		}
//		if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v2 <= 0){
//			dist_to_v2 = -1 * fabs(dist_to_v2);
//		}
//		else{
//			dist_to_v2 = fabs(dist_to_v2);
//		}
//
//		if (dist_to_v1 <= dist_to_v2){
//			return dist_to_v2;
//		}
//		else{
//			return dist_to_v1;
//		}
//
//		//// check which vertex is farther away from birth vertex in x.
//		//if (dist_to_v1 <= dist_to_v2){
//		//	// check if that vertex is in the proper orientation for the fiber to die.
//		//	if ((*nextv2_nucsite).ProbFiberDeath() > 0){
//		//		return dist_to_v2;
//		//	}
//		//	else{
//		//		// if the original vertex is not in proper orientation for fiber death,
//		//		sort(xpos.begin(), xpos.end());		//sorts the vector min to max
//		//		int find_vert_ind = -1;
//		//		double mindist = (*cell).GetMaxHorizontalLength();
//
//		//		// if geometric center is to the left of the nucleating fiber,
//		//		if ((geocent_this_cell.x - pos_nucleating_site.x) <= 0){
//
//		//			double next_x_pos;
//		//			double dist_to_death_vert;
//		//			// Store the index which corresponds to nextv2_pos
//		//			for (int vertex_index1 = 0; vertex_index1 != corresponding_x_positions.size(); vertex_index1++){
//		//				if (corresponding_x_positions[vertex_index1] - nextv2_pos.x <= 0 && fabs(corresponding_x_positions[vertex_index1] - nextv2_pos.x) <= mindist){
//		//					mindist = fabs(corresponding_x_positions[vertex_index1] - nextv2_pos.x);
//		//					find_vert_ind = vertex_index1;
//		//					next_x_pos = corresponding_x_positions[vertex_index1];
//		//				}
//		//			}
//
//		//			if (find_vert_ind == -1){
//		//				next_x_pos = nextv2_pos.x;
//		//			}
//		//			dist_to_death_vert = next_x_pos - (pos_nucleating_site).x;
//
//		//			// make this x distance signed (positive if death vertex is on the side closest to geometric center of cell).
//		//			if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_death_vert <= 0){
//		//				dist_to_death_vert = -1 * fabs(dist_to_death_vert);
//		//			}
//		//			else{
//		//				dist_to_death_vert = fabs(dist_to_death_vert);
//		//			}
//
//		//			return dist_to_death_vert;
//		//		}
//		//		else{
//
//		//			double next_x_pos;
//		//			double dist_to_death_vert;
//		//			// Store the index which corresponds to nextv2_pos
//		//			for (int vertex_index1 = 0; vertex_index1 != corresponding_x_positions.size(); vertex_index1++){
//		//				if (corresponding_x_positions[vertex_index1] - nextv2_pos.x >= 0 && fabs(corresponding_x_positions[vertex_index1] - nextv2_pos.x) <= mindist){
//		//					mindist = fabs(corresponding_x_positions[vertex_index1] - nextv2_pos.x);
//		//					find_vert_ind = vertex_index1;
//		//					next_x_pos = corresponding_x_positions[vertex_index1];
//		//				}
//		//			}
//
//		//			if (find_vert_ind == -1){
//		//				next_x_pos = nextv2_pos.x;
//		//			}
//
//		//			dist_to_death_vert = next_x_pos - (pos_nucleating_site).x;
//
//		//			// make this x distance signed (positive if death vertex is on the side closest to geometric center of cell).
//		//			if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_death_vert <= 0){
//		//				dist_to_death_vert = -1 * fabs(dist_to_death_vert);
//		//			}
//		//			else{
//		//				dist_to_death_vert = fabs(dist_to_death_vert);
//		//			}
//
//		//			return dist_to_death_vert;
//		//		}
//		//	}
//		//}
//		//else{
//		//	// check if that vertex is in the proper orientation for the fiber to die.
//		//	if ((*nextv1_nucsite).ProbFiberDeath() > 0){
//		//		return dist_to_v1;
//		//	}
//		//	else{
//		//		// if the original vertex is not in proper orientation for fiber death,
//		//		sort(xpos.begin(), xpos.end()); //sorts the vector min to max
//		//		int find_vert_ind = -1;
//		//		double mindist = (*cell).GetMaxHorizontalLength();
//
//
//		//		// if geometric center is to the left of the nucleating fiber,
//		//		if ((geocent_this_cell.x - pos_nucleating_site.x) <= 0){
//
//		//			double next_x_pos;
//		//			double dist_to_death_vert;
//		//			// Store the index which corresponds to nextv2_pos
//		//			for (int vertex_index1 = 0; vertex_index1 != corresponding_x_positions.size(); vertex_index1++){
//		//				if (corresponding_x_positions[vertex_index1] - nextv1_pos.x <= 0 && fabs(corresponding_x_positions[vertex_index1] - nextv1_pos.x) <= mindist){
//		//					mindist = fabs(corresponding_x_positions[vertex_index1] - nextv1_pos.x);
//		//					find_vert_ind = vertex_index1;
//		//					next_x_pos = corresponding_x_positions[vertex_index1];
//		//				}
//		//			}
//
//		//			if (find_vert_ind == -1){
//		//				next_x_pos = nextv1_pos.x;
//		//			}
//
//		//			dist_to_death_vert = next_x_pos - (pos_nucleating_site).x;
//
//		//			// make this x distance signed (positive if death vertex is on the side closest to geometric center of cell).
//		//			if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_death_vert <= 0){
//		//				dist_to_death_vert = -1 * fabs(dist_to_death_vert);
//		//			}
//		//			else{
//		//				dist_to_death_vert = fabs(dist_to_death_vert);
//		//			}
//
//		//			return dist_to_death_vert;
//		//		}
//		//		else{
//
//		//			double next_x_pos;
//		//			double dist_to_death_vert;
//		//			// Store the index which corresponds to nextv2_pos
//		//			for (int vertex_index1 = 0; vertex_index1 != corresponding_x_positions.size(); vertex_index1++){
//		//				if (corresponding_x_positions[vertex_index1] - nextv1_pos.x >= 0 && fabs(corresponding_x_positions[vertex_index1] - nextv1_pos.x) <= mindist){
//		//					mindist = fabs(corresponding_x_positions[vertex_index1] - nextv1_pos.x);
//		//					find_vert_ind = vertex_index1;
//		//					next_x_pos = corresponding_x_positions[vertex_index1];
//		//				}
//		//			}
//
//		//			if (find_vert_ind == -1){
//		//				next_x_pos = nextv1_pos.x;
//		//			}
//
//
//		//			dist_to_death_vert = next_x_pos - (pos_nucleating_site).x;
//
//		//			// make this x distance signed (positive if death vertex is on the side closest to geometric center of cell).
//		//			if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_death_vert <= 0){
//		//				dist_to_death_vert = -1 * fabs(dist_to_death_vert);
//		//			}
//		//			else{
//		//				dist_to_death_vert = fabs(dist_to_death_vert);
//		//			}
//
//		//			return dist_to_death_vert;
//		//		}
//		//	}
//		//}
//	}
//	else{
//		return 0.0;
//	}
//}


double NucleationSite::Calc_XDist_toDeathVertex(){
	if ((*this).ProbFiberBirth() > 0){
		vector<NucleationSite*> this_fiber_nuc_site_vect = (*cell).fibernucsites;

		vector<int> is_fiber_birth_vertex_vect;
		is_fiber_birth_vertex_vect.resize(this_fiber_nuc_site_vect.size());

		vector<int> is_fiber_death_vertex_vect;
		is_fiber_death_vertex_vect.resize(this_fiber_nuc_site_vect.size());

		vector<double> corresponding_x_positions;
		vector<double> corresponding_y_positions;

		// Get the list of ordered vertex positions
		vector<vector<double>> xypos = (*cell).HelperOrderVerticesPos();
		vector<double> xpos = xypos[0];
		vector<double> ypos = xypos[1];
		Vect geocent_this_cell = GetGeoCenter(xpos, ypos);

		for (int nuc_site_index = 0; nuc_site_index != this_fiber_nuc_site_vect.size(); nuc_site_index++){

			NucleationSite* this_nuc_site_ptr = this_fiber_nuc_site_vect[nuc_site_index];
			Vertex* corresponding_vertex_ptr = (*this_nuc_site_ptr).parentvertex;
			Vect this_vertex_pos = (*corresponding_vertex_ptr).r;

			for (int vert_index = 0; vert_index != xpos.size(); vert_index++){
				if (xpos[vert_index] == this_vertex_pos.x && ypos[vert_index] == this_vertex_pos.y){

					is_fiber_birth_vertex_vect[vert_index] = 0;
					is_fiber_death_vertex_vect[vert_index] = 0;

					if ((*this_nuc_site_ptr).ProbFiberDeath() > 0){
						is_fiber_death_vertex_vect[vert_index] = 1;
						corresponding_x_positions.push_back((*corresponding_vertex_ptr).r.x);
						corresponding_y_positions.push_back((*corresponding_vertex_ptr).r.y);
					}

					if ((*this_nuc_site_ptr).ProbFiberBirth() > 0){
						is_fiber_birth_vertex_vect[vert_index] = 1;
					}
				}
			}
		}

		// calculate current position of vertex (parent)
		Vect pos_nucleating_site = (*parentvertex).r;

		// find next vertex that the fiber will run into along this edge
		Vertex* nextv1 = (*edge1).GetOtherVertex(parentvertex);
		Vertex* nextv2 = (*edge2).GetOtherVertex(parentvertex);

		// store corresponding fiber nucleation sites
		NucleationSite* nextv1_nucsite = NULL;
		NucleationSite* nextv2_nucsite = NULL;

		for (int nuc_site_index = 0; nuc_site_index != this_fiber_nuc_site_vect.size(); nuc_site_index++){
			NucleationSite* this_fiber_nuc_site = this_fiber_nuc_site_vect[nuc_site_index];
			Vertex* corresponding_vertex_ptr = (*this_fiber_nuc_site).parentvertex;

			if ((*corresponding_vertex_ptr).id == (*nextv1).id){
				nextv1_nucsite = this_fiber_nuc_site;
			}

			if ((*corresponding_vertex_ptr).id == (*nextv2).id){
				nextv2_nucsite = this_fiber_nuc_site;
			}
		}

		// find the positions of next vertices
		Vect nextv1_pos = (*nextv1).r;
		Vect nextv2_pos = (*nextv2).r;

		// calculate distance along x to each of these vertices
		double dist_to_v1 = (nextv1_pos - pos_nucleating_site).x;
		double dist_to_v2 = (nextv2_pos - pos_nucleating_site).x;

		// make this x distance signed (positive if death vertex is on the side closest to geometric center of cell).
		if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v1 <= 0){
			dist_to_v1 = -1 * fabs(dist_to_v1);
		}
		else{
			dist_to_v1 = fabs(dist_to_v1);
		}
		if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v2 <= 0){
			dist_to_v2 = -1 * fabs(dist_to_v2);
		}
		else{
			dist_to_v2 = fabs(dist_to_v2);
		}

		int number_of_death_vertices = corresponding_x_positions.size();
		double dist_to_death_vert;
		double next_x_pos;

		if (number_of_death_vertices == 0){
			if (dist_to_v1 <= dist_to_v2){
				dist_to_death_vert = dist_to_v2;
			}
			else{
				dist_to_death_vert = dist_to_v1;
			}
		}
		else if (number_of_death_vertices == 1){
			next_x_pos = corresponding_x_positions[0];

			dist_to_death_vert = next_x_pos - pos_nucleating_site.x;

			if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_death_vert <= 0){
				dist_to_death_vert = -1 * fabs(dist_to_death_vert);
			}
			else{
				dist_to_death_vert = fabs(dist_to_death_vert);
			}
		}
		else if (number_of_death_vertices == 2){
			nextv1_pos.x = corresponding_x_positions[0];
			nextv1_pos.y = corresponding_y_positions[0];
			nextv2_pos.x = corresponding_x_positions[1];
			nextv2_pos.y = corresponding_y_positions[1];

			dist_to_v1 = (nextv1_pos - pos_nucleating_site).x;
			dist_to_v2 = (nextv2_pos - pos_nucleating_site).x;

			if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v1 <= 0){
				dist_to_v1 = -1 * fabs(dist_to_v1);
			}
			else{
				dist_to_v1 = fabs(dist_to_v1);
			}
			if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v2 <= 0){
				dist_to_v2 = -1 * fabs(dist_to_v2);
			}
			else{
				dist_to_v2 = fabs(dist_to_v2);
			}

			if ((corresponding_y_positions[0] - (geocent_this_cell).y) * (corresponding_y_positions[1] - (geocent_this_cell).y) < 0){
				if (dist_to_v1 <= dist_to_v2){
					dist_to_death_vert = dist_to_v2;
				}
				else{
					dist_to_death_vert = dist_to_v1;
				}
			}
			else{
				if (dist_to_v1 <= dist_to_v2){
					dist_to_death_vert = dist_to_v1;
				}
				else{
					dist_to_death_vert = dist_to_v2;
				}
			}
		}
		else{
			vector<int> is_above_nucleating_vertex;
			is_above_nucleating_vertex.resize(corresponding_y_positions.size());
			double min_x_above = (*cell).GetMaxHorizontalLength() + 2;
			double min_x_below = (*cell).GetMaxHorizontalLength() + 2;

			for (int death_site_index = 0; death_site_index != corresponding_y_positions.size(); death_site_index++){
				is_above_nucleating_vertex[death_site_index] = 0;
				if ((corresponding_y_positions[death_site_index] - (geocent_this_cell).y) > 0){
					is_above_nucleating_vertex[death_site_index] = 1;

					dist_to_v1 = corresponding_x_positions[death_site_index] - pos_nucleating_site.x;

					if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v1 <= 0){
						dist_to_v1 = -1 * fabs(dist_to_v1);
					}
					else{
						dist_to_v1 = fabs(dist_to_v1);
					}

					if (min_x_above >= dist_to_v1){
						min_x_above = dist_to_v1;
					}
				}
				else{
					dist_to_v2 = corresponding_x_positions[death_site_index] - pos_nucleating_site.x;

					if ((geocent_this_cell.x - pos_nucleating_site.x) * dist_to_v2 <= 0){
						dist_to_v2 = -1 * fabs(dist_to_v2);
					}
					else{
						dist_to_v2 = fabs(dist_to_v2);
					}

					if (min_x_below >= dist_to_v2){
						min_x_below = dist_to_v2;
					}
				}
			}
			if (min_x_above == (*cell).GetMaxHorizontalLength() + 2){
				dist_to_death_vert = min_x_below;
			}
			else if (min_x_below == (*cell).GetMaxHorizontalLength() + 2){
				dist_to_death_vert = min_x_above;
			}
			else{
				if (min_x_above >= min_x_below){
					dist_to_death_vert = min_x_above;
				}
				else{
					dist_to_death_vert = min_x_below;
				}
			}
		}
		return dist_to_death_vert;
	}
	else{
		return 0.0;
	}
}