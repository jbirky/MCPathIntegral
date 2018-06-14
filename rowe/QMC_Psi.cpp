// Final Project
// Quantum Harmonic Oscillator using Metropolis
/////////////////
// Trevor Rowe
// A11929779
/////////////////

#include <iostream>
#include <math.h>
#include <cmath>
#include <random>
using namespace std;

double normal_rand(double mean,double stdev);
int uniform_rand_int(int min,int max);
double uniform_rand(double min,double max);
double path_energy(std::vector<double> &vec, double delta_tau);
double energy_op(std::vector<double> &xvec);
double metropolis(std::vector<double> &vec,int index,double stdev,double delta_tau);


int N = 10;     // Lattice points

int main()
{
    int burn = 10;       // Throw away steps
    int M = 100;         // Number of path manipulations
    double init_rng = 4.0;   // range over which to pick initial values in path (+/-)
    double T = 0.1;          // Temperature to probe
    
    double beta = 1.0/T;        // Coldness
    double dtau = beta/N;       // Step size in imaginary time
    double delta = 2.0*sqrt(dtau);  // Sample range in metropolis

    std::vector<double> psisquare;    // Approximate ground state |psi|^2

    std::vector<double> x_coord;    // Coordinates of each quantum path

    // Initialize vector
    for(int i = 0; i < N; i++)
    {
        x_coord.push_back(0.00);
        psisquare.push_back(x_coord[i]);
    }

    // Manipulate paths using Quantum Path Integral Monte Carlo
    for(int m = 0; m < M; m++)
    {
        for (int i = 0; i < N; i++)
        {
            //int rand_ind = uniform_rand_int(0,9);
            x_coord[i] = metropolis(x_coord,i,delta,dtau);
            psisquare.push_back(x_coord[i]);
        }
        
    }

    // Output ground state wave function
    for(int i = burn*N; i < psisquare.size(); i++)
    {
        std::cout << psisquare[i] << std::endl;
    }

    cout<<"dtau "<<dtau<<endl;
    cout<<"size "<<psisquare.size()<<endl;

    
    return 0;
}


double path_energy(std::vector<double> &xvec, double delta_tau)
{
    double result = 0.00;
    
    for (int i = 1; i < N; i ++)
    {
        result += (0.5/delta_tau)*pow(xvec[i] - xvec[i - 1],2.0) + delta_tau*0.5*pow(0.5*(xvec[i] + xvec[i - 1]),2.0);
    }
    result += (0.5/delta_tau)*pow(xvec[0] - xvec[N - 1],2.0) + delta_tau*0.5*pow(0.5*(xvec[0] + xvec[N - 1]),2.0);
    
    return result;
}

double energy_op(std::vector<double> &xvec)
{
    double f = 0.00;
    for(int i = 0; i < N; i++)
    {
        f += pow(xvec[i],2.0);
    }
    return f/N;
}

// Implementation of Metropolis-Hastings algorithm
double metropolis(std::vector<double> &vec,int index,double range,double delta_tau)
{
    double init = vec[index];
    
    // Generate candidate sample from nomral distribution centered on init value
    double candidate = uniform_rand(init - range,init + range);
    
    // Create vector of all same coordinates except new candidate sample
    std::vector<double> cvec = vec;
    cvec[index] = candidate;
    
    // Accept if candidate occupies point of higher probability
    double boltz_candidate = exp(-path_energy(cvec,delta_tau));
    double boltz_current = exp(-path_energy(vec,delta_tau));
    
    if(boltz_candidate >= boltz_current)
    {
        return candidate;
    }
    
    // If candidate occupies point of lower probability
    // accept with P = ratio of boltzmann weights
    else
    {
        double test_number = uniform_rand(0.0,1.0);
        double ratio = boltz_candidate/boltz_current;
        if(test_number <= ratio)
        {
            return candidate;
        }
        else
        {
            return init;
        }
    }
}


double normal_rand(double mean,double stdev)
{
    static std::default_random_engine   generator;
    std::normal_distribution<double>    normal_dist(mean, stdev);
    
    return normal_dist(generator);
}

double uniform_rand(double min,double max)
{
    static std::default_random_engine         generator;
    std::uniform_real_distribution<double>    uni_dist(min, max);
    
    return uni_dist(generator);
}

int uniform_rand_int(int min,int max)
{
    static std::default_random_engine         generator;
    std::uniform_int_distribution<int>    uni_dist(min, max);
    
    return uni_dist(generator);
}
