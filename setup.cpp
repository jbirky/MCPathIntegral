#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <math.h>
#include <cblas.h>
using namespace std;

typedef complex<double> cdouble;
typedef vector<cdouble> matrix;


// ===================================
// DEFINE PROBLEM CONSTANTS
// ===================================

// Physical constants
double  PI 	= 3.141592652;
double	h	= 1;				// Plancks constant
double	m 	= 1;				// Mass of particle
double	w	= 1; 				// Oscillation frequency
double  kb  = 1;				// Boltzmann constant
double  x0  = 0;				// x starting grid point
double  T   = .1;				// Total time	
double  Nd  = 100;				// Number data points
double  DEL_T = 1/Nd/T;			// Time step
// double  Tm  = 100;

// MCMC parameters
int		Ns = 100;				// Number of steps
int		Nw = 1; 				// Number of walkers
int		Nburn = pow(10,1);		// Burn in 
int 	Nskip = 1; 				// Number steps skipped per iteration
double	DEL_S = 1;				// Width of gaussian step
int 	Nsamp = round((Ns-Nburn)/Nskip); // Number of points kept per sample