// ===================================
// JESSICA BIRKY (A13002163)
// ===================================

#include "setup.cpp"
#include "utils.cpp"
using namespace std;


// ===================================
// DECLARE FUNCTIONS
// ===================================

// double x(int i);
// double harmonicAction(vector<double> x);
// vector<double> sampleEnergyDist(double T);

double energy(vector<double> x);
vector<double> returnDist();

vector<double> estimateError(vector<double> vec);
vector<double> removeAutocorrelation(vector<double> vec);

// From utils.cpp
double getRandUniform(double min, double max);						// get random uniform number
vector<double> getRandGaussian(int n);								// get n random gaussian numbers

double printMatrixR(vector<double> k);
void   saveFileR(vector<double> vec, string save_name);

vector<double> vec_sum(vector<double> m1, vector<double> m2);
vector<double> vec_diff(vector<double> m1, vector<double> m2);
double vec_avg(vector<double> vec);
double vec_std(vector<double> vec);
double dot(vector<double> m1, vector<double> m2); 				


// ===================================
// RUN SIMULATION
// ===================================

int main() {

	clock_t begin = clock();
	srand(time(0));				// random number seed
	
	// ========================================= 
	// Question 1

	// vector<double> energies = sampleEnergyDist(T);
	// vector<double> sample = removeAutocorrelation(energies);
	// string e_save = "output/energies.dat";
	// saveFileR(sample, e_save);

	// string dist_type = "classical";
	// vector<double> exp = expected(dist_type);

	vector<double> dist = returnDist();
	vector<double> sample = removeAutocorrelation(dist);
	string save = "output/dist.dat";
	saveFileR(dist, save);


	// =========================================

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / pow(10,6);
	cout<<"Time elapsed: "<<elapsed_secs<<" sec"<<endl;

	return 0;
}


// ===================================
// MONTE CARLO
// ===================================

double energy(vector<double> x) {

	double energy = 0;
	for (int i=1; i<Nd; i++) {
		energy += m/(2*DEL_T) * pow((x[i] - x[i-1]),2.) + DEL_T*m*pow(w,2)/2*pow((x[i] +  x[i-1])/2, 2.);
	}
	energy += m/(2*DEL_T) * pow((x[0] - x[Nd-1]),2.) + DEL_T*m*pow(w,2)/2*pow((x[0] +  x[Nd-1])/2, 2.);

	return energy;
}

vector<double> returnDist() {

	int k;
	int km;
	int kp;
	vector<double> energies(Ns,0);
	vector<double> x(Nd,0);
	for (int i=0; i<Nd; i++) {x[i] = getRandUniform(-2,2);}
	x[0]  = x0;
	x[-1] = x0;

	int accept = 0;

	for (int i=0; i<Ns; i++) {
		k = (rand() % (int) Nd);

		km = k - 1;
		kp = k + 1;
		if (km == -1) {
			km = Nd;
		}
		vector<double> xprop = x;
		if ((k == 0) or (k == Nd)) {
			xprop[0] = x[0] + getRandGaussian(1)[0];
			xprop[Nd] = x[Nd] + getRandGaussian(1)[0];
		} else {
			xprop[k] = x[k] + getRandGaussian(1)[0];
		}

		double Ecurr = energy(x);
		double Eprop = energy(xprop);

		// cout<<Ecurr<<" "<<Eprop<<" "<<(Eprop-Ecurr)/T<<" "<<exp((Eprop-Ecurr)/T)<<endl;

		// Metropolis-Hastings algorithm
		if (Eprop < Ecurr) {
			x = xprop;
			energies[i] = Eprop;
			accept += 1;
		} else {
			if (getRandUniform(0,1) < exp(Eprop-Ecurr)) {
				x = xprop;
				energies[i] = Eprop;
				accept += 1;
			} else {
				energies[i] = Ecurr;
			}
		}
		string xsave = "output/xvec/x" + to_string(i) + ".csv";
		saveFileR(x, xsave);
	}

	double accept_rate = ((double) accept) / ((double) Ns);
	cout<<accept_rate<<" steps accepted"<<endl;

	return x;
}

vector<double> removeAutocorrelation(vector<double> vec) {

	
	vector<double> sample(Nsamp, 0);

	for (int i=0; i<Nsamp; i++) {
		sample[i] = vec[Nburn + i*Nskip];
	}

	return sample;
}


vector<double> estimateError(vector<double> vec) {

	vector<double> sample = removeAutocorrelation(vec);

	double avg = vec_avg(sample);
	double std = vec_std(sample) /sqrt((double) sample.size());

	vector<double> error = {avg, std};

	return error;
}

