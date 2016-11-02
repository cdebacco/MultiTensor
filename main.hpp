
#ifndef MAIN_H
#define MAIN_H

#include "mlg.hpp"

#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>

extern boost::mt19937 gen;

/* ----------------GRAPH VARIABLES ----------------*/

int N=500;   // Number of nodes
int L=12; 	// Number of types of edges/layers
int K=4;	// Number of groups/communities
int N_real=1;   // Number of different initial realizations --> to tackle issue of getting stuck in local minima

/* ----------------CONVERGENCE VARIABLES ----------------*/

double tolerance =0.1;    // difference btw old and new values of the convergence variables
int decision = 10;            // number of times the convergence criterium is satisfied before stopping the simulation
int maxit=500;                

unsigned int rseed=0;  // seed for gen ---> used to pick a random real number using real01(). E.g. to initialize the variables u,v and w
bool out_adjacency=false;     // Flag to output single-layer adjacencies on file

double inf=1e10; 
double err_max=0.000001; // max err in the converge() routine
double err=0.1;              // magnitude of the noise added for the initialization form file input
int initialization=0;    // 0-> random initialization; 1-> u,v,w are initializes; 2 --> only w is initialized; 3 --> only u and v are initialized 
bool assortative=false;  // use the assortative model --> force W_matrix to be diagonal

/* ----------------I/O FILES ----------------*/

std::string folder="";   // input and output folder
std::string end_file="";    // output end of file
std::string adj="adjacency.dat";   // input adjacency matrix file
std::string w_file="w.dat";   // input file for the W matrix



#endif

