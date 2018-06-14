// Final Project
// Quantum Harmonic Oscillator using Metropolis - <E>(T)
/////////////////
// Trevor Rowe
// A11929779
/////////////////

#include <iostream>
#include <math.h>
#include <cmath>
#include <random>
using namespace std;

double uniform_rand(double min,double max);
double path_energy(std::vector<double> &xvec, double delta_tau,int N_val);
double energy_op(std::vector<double> &xvec,int N_val);
double metropolis(std::vector<double> &vec,int index,double range,double delta_tau,int N_val);
double standard_error(int num_steps,double mean_energy,std::vector<double> &energies);


int main()
{
    double T0 = 0.05;        // Lower value for temperature
    double Tf = 4.0;        // Upper value for temperature
    double dT = 0.05;        // Change in temperature between steps
    int N = 10;            // Initial number of lattice points
    int nf = (int)round((Tf - T0)/dT + 1);     // Final time step
    
    
    std::cout << "Temp" << "\t" << "<E>" << "\t" << "error in <E>" << std::endl;
    
    for(int n = 0; n < nf; n++)
    {
        double T = T0 + n*dT;       // Temperature
        
        // int M;      // Number of measurements
        // int burn;   // Number of measurements to throw away
        // if(T <= 2.0)
        // {
        //     M = 800000;
        //     burn = 8000;
        // }
        // if(2.0 < T <= 2.5)
        // {
        //     M = 5000000;
        //     burn = 50000;
        // }
        // if(2.5 < T)
        // {
        //     M = 15000000;
        //     burn = 150000;
        // }

        long M = pow(10,3);
        long burn = pow(10,2);
        long steps = M - burn;       // Number of recorded measurments
        double beta = 1.0/T;        // Coldness
        double dtau = beta/N;       // Step size in imaginary time
        double delta = 2.0*sqrt(dtau);  // Sample range in metropolis
        
        std::vector<double> E;          // Used to calculate <E>
        std::vector<double> x_coord;    // Coordinates of each quantum path
        
        // Initialize vector
        for(int i = 0; i < N; i++)
        {
            x_coord.push_back(0.00);    // cold start
        }
        
        // Manipulate paths using Quantum Path Integral Monte Carlo
        for(int m = 0; m < M; m++)
        {
            for (int i = 0; i < N; i++)
            {
                x_coord[i] = metropolis(x_coord,i,delta,dtau,N);
            }
            
            if(m > burn)
            {
                E.push_back(energy_op(x_coord,N));
            }
        }
        
         double E_mean = 0.00;      // <E>
         
         for (int i = 0; i < steps; i++)
         {
             E_mean += E[i];
         }
         E_mean = E_mean/steps;
         
         double std_error = standard_error(steps,E_mean,E);
         
         std::cout << T << "\t" << E_mean <<  "\t" << std_error << std::endl;
    }
    
    return 0;
}

// Energy of quantum path
double path_energy(std::vector<double> &xvec, double delta_tau,int N_val)
{
    double result = 0.00;
    
    for (int i = 1; i < N_val; i ++)
    {
        result += (0.5/delta_tau)*pow(xvec[i] - xvec[i - 1],2.0) + delta_tau*0.5*pow(0.5*(xvec[i] + xvec[i - 1]),2.0);
    }
    result += (0.5/delta_tau)*pow(xvec[0] - xvec[N_val - 1],2.0) + delta_tau*0.5*pow(0.5*(xvec[0] + xvec[N_val - 1]),2.0);

    return result;
}

// <E>
double energy_op(std::vector<double> &xvec,int N_val)
{
    double f = 0.00;
    for(int i = 0; i < N_val; i++)
    {
        f += pow(xvec[i],2.0);
    }
    return f/N_val;
}

// Implementation of Metropolis-Hastings algorithm
double metropolis(std::vector<double> &vec,int index,double range,double delta_tau,int N_val)
{
    double init = vec[index];
    
    // Generate candidate sample from nomral distribution centered on init value
    double candidate = uniform_rand(init - range,init + range);
    
    // Create vector of all same coordinates except new candidate sample
    std::vector<double> cvec = vec;
    cvec[index] = candidate;
    
    // Accept if candidate occupies point of higher probability
    double boltz_candidate = exp(-path_energy(cvec,delta_tau,N_val));
    double boltz_current = exp(-path_energy(vec,delta_tau,N_val));
    
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

// Calculation of standard error at different bin sizes
double standard_error(int num_steps,double mean_energy,std::vector<double> &energies)
{
    std::vector<int> divisors;
    for(int i = 2; i < num_steps; i++)
    {
        if((num_steps % i) == 0)
        {
            divisors.push_back(i);
        }
    }
    int maxbin = divisors.back();
    
    double max_error = 0.00;
    int max_index = 0;
    
    for(int iter = 0; iter < divisors.size(); iter++)
    {
        double binsize = divisors[iter];
        
        // Calculate mean of each bin
        std::vector<double> bin_mean;
        for(int bin = 0; bin < num_steps/binsize; bin++)
        {
            double temp_mean = 0.00;
            for(int i = 0; i < binsize; i++)
            {
                temp_mean+= energies[bin*binsize + i];
            }
            bin_mean.push_back(temp_mean/binsize);
        }
        
        // Calculate error of binned means
        double binned_error = 0.00;
        for(int i = 0; i < bin_mean.size(); i++)
        {
            binned_error += pow(bin_mean[i] - mean_energy,2.0);
        }
        
        double standard_error = sqrt((binned_error/(bin_mean.size() - 1))/bin_mean.size());
        
        if(standard_error > max_error)
        {
            max_error = standard_error;
            max_index = iter;
        }
    }
    return max_error;
}

// Pseudo random number generator
// Generates a random real number between min and max
double uniform_rand(double min,double max)
{
    static std::default_random_engine         generator;
    std::uniform_real_distribution<double>    uni_dist(min, max);
    
    return uni_dist(generator);
}
