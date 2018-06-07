
// ===================================
// RANDOM NUMBER FUNCTIONS
// ===================================

double getRandUniform(double min, double max) {

    // Random uniform distrubution sample

    return (max - min) * ((double)rand() / (double) RAND_MAX) + min;
}


vector<double> gaussianSample2D() {

    // Get 2D Gaussian sample using Box Muller method

    double U1 = getRandUniform(0,1);
    double U2 = getRandUniform(0,1);

    double theta = 2*PI*U1;
    double R = sqrt(-2*log(U2));

    vector<double> sample = {R*DEL_S*cos(theta), R*DEL_S*sin(theta)};

    return sample;
}


vector<double> getRandGaussian(int n) {

    vector<double> sample(n,0);
    for (int i=0; i<n; i++) {
        sample[i] = gaussianSample2D()[0];
    }

    return sample;
}


// ===================================
// MATRIX/VECTOR OPERATION FUNCTIONS
// ===================================

cdouble printMatrixC(matrix k) 
{
	int n = k.size();

	for (int i=0; i<n; i++) {
		cout<<k[i]<<endl;
	}
	return 0;
}


double printMatrixR(vector<double> k) 
{
	int n = k.size();

	for (int i=0; i<n; i++) {
		cout<<k[i]<<endl;
	}
	return 0;
}


void saveFileC(vector<cdouble> vec, string save_name) 
{
	ofstream outfile;
    outfile.open(save_name);

    int size = vec.size();
    for (int i=0; i<size; i++) {
    	outfile << vec[i] << endl;
    }

    outfile.close();
}


void saveFileR(vector<double> vec, string save_name) 
{
	ofstream outfile;
    outfile.open(save_name);

    int size = vec.size();
    for (int i=0; i<size; i++) {
    	outfile << vec[i] << endl;
    }

    outfile.close();
}


vector<cdouble> matrixMatrixMultiply(vector<cdouble> m1, vector<cdouble> m2) {

    int dim = sqrt(m1.size());

    //Create a length D^2 array of doubles, filled with the value 0.
    vector<cdouble> matrixnew(dim*dim,0);

    //cblas zgemm takes in three matrices: A,B,C. It stores in the value C
    //the matrix alpha*AB+beta*C. In this case beta=0 and alpha=1, so we just
    //do matrix multiplication AB.

    cdouble mult(1, 0.0);
    cdouble czero(0.0, 0.0);

    //use complex matrix matrix multiply.
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &mult, &m1[0], dim, &m2[0], dim, &czero, &matrixnew[0], dim);

    //return the vector object.
    return matrixnew;
}


vector<cdouble> matrixVectorMultiply(vector<cdouble> m1, vector<cdouble> m2) {

    int dim = m2.size();

	vector<cdouble> prod(dim,0);

    for (int i=0; i<dim; i++) {
    	prod[i] = 0;
    	for (int j=0; j<dim; j++) {
    		prod[i] += m1[i*dim + j] * m2[j];
    	}
    }

    return prod;
}


vector<cdouble> vectorVectorMultiply(vector<cdouble> m1, vector<cdouble> m2) {

    int dim = m1.size();

	vector<cdouble> prod(dim,0);

    for (int i=0; i<dim; i++) {
    	prod[i] += m1[i] * m2[i];
    }

    return prod;
}


double dot(vector<double> m1, vector<double> m2) {

    int dim = m1.size();

    double prod = 0;

    for (int i=0; i<dim; i++) {
        prod += m1[i] * m2[i];
    }

    return prod;
}


vector<double> vec_sum(vector<double> m1, vector<double> m2) {

    int dim = m1.size();

    vector<double> sum(dim,0);

    for (int i=0; i<dim; i++) {
        sum[i] += m1[i] - m2[i];
    }

    return sum;
}


vector<double> vec_diff(vector<double> m1, vector<double> m2) {

    int dim = m1.size();

    vector<double> diff(dim,0);

    for (int i=0; i<dim; i++) {
        diff[i] += m1[i] - m2[i];
    }

    return diff;
}


double vec_avg(vector<double> vec) {

    double sum = 0;
    for (int i=0; i<vec.size(); i++) {
        sum += vec[i];
    }
    double avg = sum/vec.size();

    return avg;
}


double vec_std(vector<double> vec) {

    double sum_sq;
    double sq_sum;
    int N = vec.size();

    for (int i=0; i<N; i++) {
        sq_sum += pow(vec[i],2)/N;
        sum_sq += vec[i]/N;
    }
    sum_sq = pow(sum_sq, 2);

    double var = sq_sum - sum_sq;

    return sqrt(var);

}