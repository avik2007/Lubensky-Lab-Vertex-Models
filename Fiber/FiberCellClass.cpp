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
#include <math.h>
#include <cmath> 
#include "Constants.h"
#include "FiberConstants.h"
#include "VertexClass.h"
#include "CellClass.h"
#include "EdgeClass.h"
#include "TissueClass.h"
#include "FiberCellClass.h"
#include "FiberClass.h"
#include "FiberEdgeClass.h"
#include "FiberTissueClass.h"
#include "NucleationSite.h"

using namespace std;

//! Removes fiber from the cell.
void FiberCell::RemoveFiber(Fiber* fptr){
	int i;
	for (i = 0; i<fibers.size(); i++){
		if ((*fibers[i]).id == (*fptr).id){
			fibers.erase(fibers.begin() + i);
		}
	}
}

//returns a 2 element vector where the first number is the left side of the cell and the second is the right side obaying pbc
vector<double> FiberCell::GetXExtent(){
	vector<Vertex*> verts;					//vertices in the cell
	vector<Vertex*> firsthalfverts;			//verties on one side of a cell which crosses the x bondry 
	vector<Vertex*> secondhalfverts;		//vertices on the other side of a cell which crosses the x bondry 
	vector<double> xpos;					//x positions of the vertices
	vector<double> firstxpos, secondxpos; 
	bool isflag = Cell::GetXFlag();			//true if the cell crosses the x boundary 
	bool edgeflag;
	int i, count;
	Edge* eptr;
	vector<double> result;

	//If the cell does not cross the X-boundary
	if (!isflag){
		verts = Cell::GetVertices();		//get the vertices
		//get the x pos of each vertex
		for (i = 0; i < verts.size(); i++){
			xpos.push_back((*verts[i]).r.x);
		}
		//find the min and max
		auto min = min_element(xpos.begin(), xpos.end());
		auto max = max_element(xpos.begin(), xpos.end());
		result.push_back(*min); 
		result.push_back(*max);
	}

	//If the cell does cross the x boundary
	else{
		verts = Cell::GetClockWiseVertices();
		edgeflag = 0;
		count = 0;
		//determine which vertices are on which side of the boundary
		while (!edgeflag){
			firsthalfverts.push_back(verts[count]); //add the vertex to this half
			eptr = (*verts[count]).GetEdgeBetween(verts[count + 1]); //find out if the next vertex crosses the boundary
			if ((*eptr).xflag){ edgeflag = 1; }
			else{ edgeflag = 0; }
			count++;
		}
		edgeflag = 0;
		while (!edgeflag){
			secondhalfverts.push_back(verts[count]); //add the vertex to this half
			eptr = (*verts[count]).GetEdgeBetween(verts[(count + 1) % verts.size()]); //find out if the next vertex crosses the boundary
			if ((*eptr).xflag){ edgeflag = 1; }
			else{ edgeflag = 0; }
			count++;
		}
		//now add any remaining verts to the first half
		if (count < verts.size()){
			for (i = count; i < verts.size(); i++){ firsthalfverts.push_back(verts[i]); }
		}
		//calculate the result
		for (i = 0; i < secondhalfverts.size(); i++){ secondxpos.push_back((*secondhalfverts[i]).r.x); }
		for (i = 0; i < firsthalfverts.size(); i++){ firstxpos.push_back((*firsthalfverts[i]).r.x); }
		if (firstxpos[0] < secondxpos[0]){
			auto min = min_element(secondxpos.begin(), secondxpos.end());
			auto max = max_element(firstxpos.begin(), firstxpos.end());
			result.push_back(*min);
			result.push_back(*max);
		}
		else{
			auto min = min_element(firstxpos.begin(), firstxpos.end());
			auto max = max_element(secondxpos.begin(), secondxpos.end());
			result.push_back(*min);
			result.push_back(*max);
		}
	}

	return result;
}

//returns a 2 element vector where the first number is the bottom edge of the cell and the second is the top edge obaying pbc
vector<double> FiberCell::GetYExtent(){
	vector<Vertex*> verts;					//vertices in the cell
	vector<Vertex*> firsthalfverts;			//verties on one side of a cell which crosses the y bondry 
	vector<Vertex*> secondhalfverts;		//vertices on the other side of a cell which crosses the y bondry 
	vector<double> ypos;					//y positions of the vertices
	vector<double> firstypos, secondypos;
	bool isflag = Cell::GetYFlag();			//true if the cell crosses the x boundary 
	bool edgeflag;
	int i, count;
	Edge* eptr;
	vector<double> result;

	//If the cell does not cross the Y-boundary
	if (!isflag){
		verts = Cell::GetVertices();		//get the vertices
		//get the y pos of each vertex
		for (i = 0; i < verts.size(); i++){
			ypos.push_back((*verts[i]).r.y);
		}
		//find the min and max
		auto min = min_element(ypos.begin(), ypos.end());
		auto max = max_element(ypos.begin(), ypos.end());
		result.push_back(*min);
		result.push_back(*max);
	}

	//If the cell does cross the y boundary
	else{
		verts = (*this).GetClockWiseVertices();
		edgeflag = 0;
		count = 0;
		//determine which vertices are on which side of the boundary
		while (!edgeflag){
			firsthalfverts.push_back(verts[count]); //add the vertex to this half
			eptr = (*verts[count]).GetEdgeBetween(verts[count + 1]); //find out if the next vertex crosses the boundary
			if ((*eptr).yflag){ edgeflag = 1; }
			else{ edgeflag = 0; }
			count++;
		}
		edgeflag = 0;
		while (!edgeflag){
			secondhalfverts.push_back(verts[count]); //add the vertex to this half
			eptr = (*verts[count]).GetEdgeBetween(verts[(count + 1) % verts.size()]); //find out if the next vertex crosses the boundary
			if ((*eptr).yflag){ edgeflag = 1; }
			else{ edgeflag = 0; }
			count++;
		}
		//now add any remaining verts to the first half
		if (count < verts.size()){
			for (i = count; i < verts.size(); i++){ firsthalfverts.push_back(verts[i]); }
		}
		//calculate the result
		for (i = 0; i < secondhalfverts.size(); i++){ secondypos.push_back((*secondhalfverts[i]).r.y); }
		for (i = 0; i < firsthalfverts.size(); i++){ firstypos.push_back((*firsthalfverts[i]).r.y); }
		if (firstypos[0] < secondypos[0]){
			auto min = min_element(secondypos.begin(), secondypos.end());
			auto max = max_element(firstypos.begin(), firstypos.end());
			result.push_back(*min);
			result.push_back(*max);
		}
		else{
			auto min = min_element(firstypos.begin(), firstypos.end());
			auto max = max_element(secondypos.begin(), secondypos.end());
			result.push_back(*min);
			result.push_back(*max);
		}
	}

	return result;
}

double inverseGaussian(double mu, double lambda, mt19937_64* generator){
	normal_distribution<> normaldist(0, 1);  //creates a normal distrabution with a mean of 0 and stdv of 1
	uniform_real_distribution<> uniformdist(0, 1); //creates a uniform distrabution on the interval 0 to 1

	double v = normaldist(*generator);			//sample from a normal distrabution with mu=0 and stdev=1
	double y = v*v;
	double x = mu + (mu*mu*y) / (2 * lambda) - (mu / (2 * lambda)) * sqrt(4 * mu*lambda*y + mu*mu*y*y);
	double test = uniformdist(*generator);  // sample from a uniform distribution between 0 and 1
	if (test <= (mu) / (mu + x)){ return x; }
	else { return (mu*mu) / x; }
}
//Sets the isCompleete flag
bool FiberCell::SetIsCompleeteFlag(){
	int edgenumber;
	Edge* e;
	Edge* enext = NULL;
	Vertex* v;
	Vertex* vnext;
	int firstv;
	int index = 0;

	edgenumber = (*this).edges.size();
	//If the cell has two or less edges then is is definatly not compleete 
	if (edgenumber < 3){
		(*this).isCompleete = false;
		return false; //set the is complete flag to false and exit
	}

	//If it has some edges
	if ((*this).edges[0] == NULL){
		return false;
	}
	e = (*this).edges[0]; //get the first edge
	v = (*e).vertex1;
	firstv = (*v).id; //get the first vertex and record its id

	//If it has a bad vertex 
	if ((*v).edges.size() > 4){
		(*this).isCompleete = false;
		return false; //set the is complete flag to false and exit
	}

	//Keep walking around the cell untill you either get back to the first vertex or run out of edges
	while (true){
		vnext = (*v).GetOtherVertex(e); //Get the next vertex
		if ((*vnext).id == firstv){		//if it is the same as the starting vertex the cell is closed and we are done
			(*this).isCompleete = true;
			return true;
		}
		//if not keep going
		if (GetNextEdge(e, vnext, &index)){ //if there is a next edge reset for the next run
			enext = (*vnext).edges[index];
			v = vnext;
			e = enext;
		}
		else{ //if there is no next edge then stop
			(*this).isCompleete = false;
			return false;
		}
	}
}

//takes and edge and a vertex on that edge and returns true if that vertex has another edge in this cell
bool FiberCell::GetNextEdge(Edge* e, Vertex* v, int* result){
	bool flag;
	Edge* teste;
	int cid = (*this).id;

	flag = false;
	//for every edge attached to this vertex (there must be at least one because you gave it e)
	for (int i = 0; i < (*v).edges.size(); i++){
		teste = (*v).edges[i];
		if ((*e).id != (*teste).id){		//if the edge is not the same one you put in 
			if ((*teste).cell1 != NULL && (*teste).cell2 != NULL){
				if (cid == (*(*teste).cell1).id || cid == (*(*teste).cell2).id){  //and one of its cells matches this cell
					flag = true;
					(*result) = i;
				}
			}
		}
	}

	return flag;
}


//returns the probablilty of a fiber neucleating from any vertex in the cell 
double FiberCell::CountNucleatingSites(){
	double count_of_nuc_sites = 0;
	double prob_birth = 0.0;

	for (int i = 0; i < fibernucsites.size(); i++){
		prob_birth = (*(fibernucsites[i])).ProbFiberBirth();

		if (prob_birth > 0.0){
			count_of_nuc_sites++;
		}
	}

	return count_of_nuc_sites;
}

double FiberCell::CountFiberDeathSites(){
	double count_of_death_sites = 0;
	double prob_death = 0.0;

	for (int i = 0; i < fibernucsites.size(); i++){
		prob_death = (*(fibernucsites[i])).ProbFiberDeath();

		if (prob_death > 0.0){
			count_of_death_sites++;
		}
	}

	return count_of_death_sites;
}

//returns the maximum horizontal extent of the cell
double FiberCell::GetMaxHorizontalLength(){
	if ((*this).isCompleete){
		double max_x_length = 0.0;
		vector<vector<double>> vert_pos_vect = (*this).HelperOrderVerticesPos();
		vector<double> xpos_vect = vert_pos_vect[0];
		vector<double> minmax_xpos_vect;

		vector<double> checkminmax = (*this).GetXExtent();

		//find the min and max
		auto min = min_element(xpos_vect.begin(), xpos_vect.end());
		auto max = max_element(xpos_vect.begin(), xpos_vect.end());
		minmax_xpos_vect.push_back(*min);
		minmax_xpos_vect.push_back(*max);
		

		//Return the maximum x extent of the cell
		max_x_length = minmax_xpos_vect[1] - minmax_xpos_vect[0];
		return max_x_length;
	}
	else{
		return 0.0;
	}
}

double FiberCell::GetMaxVerticalLength(){
	if ((*this).isCompleete){
		double max_y_length = 0.0;
		vector<vector<double>> vert_pos_vect = (*this).HelperOrderVerticesPos();
		vector<double> ypos_vect = vert_pos_vect[1];
		vector<double> minmax_ypos_vect;

		vector<double> checkminmax = (*this).GetYExtent();

		//find the min and max
		auto min = min_element(ypos_vect.begin(), ypos_vect.end());
		auto max = max_element(ypos_vect.begin(), ypos_vect.end());
		minmax_ypos_vect.push_back(*min);
		minmax_ypos_vect.push_back(*max);


		//Return the maximum x extent of the cell
		max_y_length = minmax_ypos_vect[1] - minmax_ypos_vect[0];
		return max_y_length;
	}
	else{
		return 0.0;
	}
}

//look at all nucleating sites (with appropriate opening angle) in a given cell, and calculate average distance along x direction to death TCJ
double FiberCell::CalcAverageXDist_toDeathVertex(){
	//will calculate mean distance along x to death TCJ using this.
	double mean_dist_x = 0.0;

	// if cell is not a weird edge cell,
	if ((*this).isCompleete){
		//count the total number of nucleating sites (= number of vertices)
		int number_nuc_sites = fibernucsites.size();
		//initialize a vector in which we will store the distance to the death vertex
		vector<double> distance_to_death_vertex;


		// loop through all vertices in this cell
		for (int nuc_sit_index = 0; nuc_sit_index < number_nuc_sites; nuc_sit_index++){
			// store pointer to nucleating site.
			NucleationSite* pointer_to_this_site = fibernucsites[nuc_sit_index];
			// check that nuc site has proper opening angle
			if ((*pointer_to_this_site).ProbFiberBirth() > 0){
				// push back the distance along x to death TCJ
				distance_to_death_vertex.push_back((*pointer_to_this_site).Calc_XDist_toDeathVertex());
			}
		}
		// store the size of the above vector (number of nucleating sites with proper opening angle).
		int distance_to_death_vertex_size = distance_to_death_vertex.size();
		// if there exists at least one nucleating site with proper opening angle,
		if (distance_to_death_vertex_size > 0){
			// calculate the mean distance along x to death TCJ
			for (int vect_index = 0; vect_index < distance_to_death_vertex_size; vect_index++){
				mean_dist_x += distance_to_death_vertex[vect_index] / distance_to_death_vertex_size;
			}
		}
	}
	// then return mean distance to death TCJ
	return mean_dist_x;
}