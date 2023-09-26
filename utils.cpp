#include "utils.h"

vector<double> CSVtoVector(string filename)
{
    vector<double> input_vec;
    // vector<string> string_vector;
    ifstream data(filename);
    string line;
    // int line_count = 0;
    // string cell;

    while (getline(data, line))
    {
        // string_vector.push_back(line);
        double entry = ::stof(line);
        input_vec.push_back(entry);
        // cout << weight << " ";
    }
    return input_vec;
}
vector<double> getRowFromMatrixFile(string& filename, int rowIndex, int numCol) {
    vector<double> rowVector;
    ifstream data(filename);
    string line;
    int currentRow = 0;
    while (getline(data, line)) {
        if (currentRow == 0) {
            // Skip the first row (header)
            currentRow++;
            continue;
        }
        if (currentRow-1 == rowIndex) {
            stringstream lineStream(line);
            string cell;
            int currentColumn = 0;

            while (getline(lineStream, cell, '\t')) 
            {
                if (currentColumn == 0) {
                    currentColumn++;
                    continue;
                }
                if (currentColumn > numCol){
                    break;
                }
                try {
                double entry = stod(cell);
                rowVector.push_back(entry);
                }
                catch (const exception& e) {
                    cerr << "Exception caught: " << e.what() << std::endl;
                }
                currentColumn++;
            }
        }
        currentRow++;
    }
    return rowVector;
}
void read_bedfile_row(vector<double>& rowData, string& geneID, const string& filename, int row, int skipcols, bool header) {
    // vector<double> rowData;
    ifstream data(filename);
    string line;
    // Skip the header if required
    if (header)
        getline(data, line);
    // Skip rows until the desired row
    for (int currentRow = 0; currentRow < row; ++currentRow) {
        if (!getline(data, line)) {
            cerr << "Desired row not found." << endl;
        }
    }
    stringstream lineStream(line);
    string cell;
    // Skip the first N columns
    for (int i = 0; i < skipcols; ++i) {
        if (!getline(lineStream, cell, '\t')) {
            cerr << "Not enough columns in the row." << endl;
        }
        if (i == skipcols-1)
            geneID = cell;

    }
    while (getline(lineStream, cell, '\t')) {
        try {
            double entry = stod(cell);
            rowData.push_back(entry);
        } catch (const exception& e) {
            cerr << "Exception caught: " << e.what() << endl;
        }
    }
    data.close();
}
vector<vector<double>> getTPMFromMatrixFile(const string& filename, vector<string>& geneID) {
    vector<vector<double>> rowsData;
    ifstream data(filename);
    string line;
    // int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line)) {
        stringstream lineStream(line);
        string cell;
        // int currentColumn = 0;

        // Skip the first two columns
        getline(lineStream, cell, '\t');
        geneID.push_back(cell);
        getline(lineStream, cell, '\t');

        vector<double> rowVector;
        while (getline(lineStream, cell, '\t')) {
            try {
                double entry = stod(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cerr << "Exception caught: " << e.what() << endl;
            }
            // currentColumn++;
        }

        rowsData.push_back(rowVector);
        // currentRow++;
    }

    // Close the file after reading
    data.close();

    return rowsData;
}
vector<vector<uint32_t>> getCountFromMatrixFile(const string& filename, vector<string>& geneID) {
    vector<vector<uint32_t>> rowsData;
    ifstream data(filename);
    string line;
    // int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line)) {
        stringstream lineStream(line);
        string cell;
        // int currentColumn = 0;

        // Skip the first two columns
        getline(lineStream, cell, '\t');
        geneID.push_back(cell);
        getline(lineStream, cell, '\t');

        vector<uint32_t> rowVector;
        while (getline(lineStream, cell, '\t')) {
            try {
                uint32_t entry = stoi(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cerr << "Exception caught: " << e.what() << endl;
            }
            // currentColumn++;
        }

        rowsData.push_back(rowVector);
        // currentRow++;
    }

    // Close the file after reading
    data.close();

    return rowsData;
}
vector<uint32_t> ScaleVector(vector<double> &v, int k)
{
    vector<uint32_t> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = v[i] * k;
    return intvec;
}
vector<vector<int32_t>> ScaleVector(vector<vector<double>>& v, int k)
{
    vector<vector<int32_t>> intvec(v.size(), vector<int32_t>(v[0].size()));
    for (int i=0; i<v.size(); i++)
    {
        for (int j=0; j<v[0].size(); j++)
        {
            intvec[i][j]=v[i][j]*k;
        }
    }
    return intvec;
}
vector<int32_t> ScaleVector_signed(vector<double> &v, int k)
{
    vector<int32_t> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = v[i] * k;
    return intvec;
}


vector<double> UnscaleVector_signed(vector<int32_t> &v, int k)
{
    vector<double> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = static_cast<double>(v[i]) / static_cast<double>(k);
    return intvec;
}
vector<uint32_t> ShiftVector(vector<int32_t>& vec, uint32_t number) 
{
    vector<uint32_t> shifted(vec.size());
    int32_t shiftedValue;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i]<0)
            shiftedValue = vec[i] + number;
        else
            shiftedValue = vec[i];
        shifted[i] = static_cast<uint32_t>(shiftedValue);
    }
    return shifted;
}

vector<vector<uint32_t>> ShiftVector(vector<vector<int32_t>>& vec, uint32_t number) 
{
    vector<vector<uint32_t>> shifted(vec.size(), vector<uint32_t>(vec[0].size()));
    int32_t shiftedValue;
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j=0; j< vec[0].size();j++)
        {
            if(vec[i][j]<0)
            {
                shiftedValue = vec[i][j] + number;
            }
            else{
                shiftedValue = vec[i][j];
            }
            shifted[i][j] = static_cast<uint32_t>(shiftedValue);
        }
        
    }
    return shifted;
}

vector<int32_t> UnshiftVector(vector<uint32_t>& shifted, uint32_t number)
{
    vector<int32_t> unshifted(shifted.size());
    for (size_t i = 0; i < shifted.size(); ++i) {
        if (shifted[i] > number/2)
            unshifted[i] = static_cast<int32_t>(shifted[i]) - number;
        else
            unshifted[i] = static_cast<int32_t>(shifted[i]);
        // unshifted[i] = static_cast<int32_t>(shiftedValue);
    }
    return unshifted;
}
uint32_t nearestPowerOf2(int N)
{
    uint32_t a = log2(N);
 
    // if (pow(2, a) == N)
    //     return N;
 
    return pow(2, a + 1);
}
double betaLogLikelihood(const gsl_vector *v, void *params) {
	double * p = (double *) params;
	double beta_shape1 = gsl_vector_get(v, 0);
	double beta_shape2 = gsl_vector_get(v, 1);
	return -1.0 * ((beta_shape1 - 1) * p[0] + (beta_shape2 - 1) * p[1] - p[2] * gsl_sf_lnbeta(beta_shape1, beta_shape2));
}

int mleBeta(vector < double > & pval, double & beta_shape1, double & beta_shape2) {

	//Set starting point to moment matching estimates
	gsl_vector * x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, beta_shape1);
	gsl_vector_set (x, 1, beta_shape2);

	//Set initial step sizes to shape1 and shape2 scales
	gsl_vector * ss = gsl_vector_alloc (2);
	gsl_vector_set (ss, 0, beta_shape1/10);
	gsl_vector_set (ss, 1, beta_shape2/10);

	//Initialize method and iterate
	double par [3];
	par[0] = 0.0;
	par[1] = 0.0;
	for (int e = 0 ; e < pval.size(); e ++) {
		if (pval[e] == 1.0) pval[e] = 0.99999999;
		par[0] += log (pval[e]);
		par[1] += log (1 - pval[e]);
	}
	par[2] = pval.size();
	gsl_multimin_function minex_func;
	minex_func.n = 2;
	minex_func.f = betaLogLikelihood;
	minex_func.params = par;

	//Initialize optimization machinery
	const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	//Optimization iteration
	size_t iter = 0;
	int status;
	double size;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 0.01);
	} while (status == GSL_CONTINUE && iter < 1000);

	//Output new beta shape values
	beta_shape1 = gsl_vector_get (s->x, 0);
	beta_shape2 = gsl_vector_get (s->x, 1);

	//Free allocated memory
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	return (status == GSL_SUCCESS);
}

double doublemean(vector < double > & X) {
		double mean = 0.0;
		for (int x = 0 ; x < X.size() ; x ++) mean += X[x];
		mean /= X.size();
		return mean;
	}

double doublevariance(vector < double > & X, double mean) {
    double variance = 0.0;
    for (int x = 0 ; x < X.size() ; x++) variance += (X[x] - mean) * (X[x] - mean);
    variance /= (X.size() - 1);
    return variance;
}
double degreeOfFreedom(const gsl_vector *v, void *params) {
    pair<vector<double>*, double>* p = static_cast<pair<vector<double>*, double>*>(params);	
    vector<double>& corr = *(p->first);
    vector < double > pval = vector < double >(corr.size(), 0.0);
    
	double mean = 0.0;
	for (int c = 0 ; c < corr.size() ; c++) {
		pval[c] = getPvalue(corr[c], gsl_vector_get(v, 0));
		mean += pval[c];
	}
	mean /= pval.size();
	double variance = 0.0;
	for (int c = 0 ; c < pval.size() ; c++) variance += (pval[c] - mean) * (pval[c] - mean);
	variance /= (pval.size() - 1);

	double shape2 = abs((mean * (mean * (1 - mean ) / variance - 1)) - 1.0);
	//cout << "O = " << mean << " " << shape2 << endl;
	return shape2;
}
int learnDF(vector < double > & corr, double & df) {

	//Set starting point to moment matching estimates
	gsl_vector * x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, df);

	//Set initial step sizes to shape1 and shape2 scales
	gsl_vector * ss = gsl_vector_alloc (1);
	gsl_vector_set (ss, 0, df * 0.1);

	pair<vector<double>*, double> params(&corr, df);

	gsl_multimin_function minex_func;
	minex_func.n = 1;
	minex_func.f = degreeOfFreedom;
	minex_func.params = &params;

	//Initialize optimization machinery
	const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 1);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	//Optimization iteration
	//cout << "\n ========================" << endl;
	size_t iter = 0;
	int status;
	double size;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 0.01);
		//printf ("%d %10.3e f() = %7.10f size = %.10f\n", iter, gsl_vector_get (s->x, 0), s->fval, size);

	} while (status == GSL_CONTINUE && iter < 20);

	//Output new beta shape values
	df = gsl_vector_get (s->x, 0);

	//Free allocated memory
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	return (status == GSL_SUCCESS);
}
double getTstat2(double corr, double df) {
    return df * corr * corr / (1 - corr * corr);
}
double getPvalueFromTstat2(double tstat2, double df) {
    return pf(tstat2, 1, df, 0, 0);
}
double getPvalue(double corr, double df) {
    return pf(df * corr * corr / (1 - corr * corr), 1, df,0,0);
}
double getSlope(double nominal_correlation, double psd, double gsd) {
    if (gsd < 1e-16 || psd < 1e-16) return 0;
    else return nominal_correlation * psd / gsd;
}
vector<double> center_normalize(vector<vector<double>>& M) {
    int rows = M.size();
    int cols = M[0].size();

    // Compute the mean and variance for each row
    vector<vector<double>> N(rows, vector<double>(cols, 0.0));
    vector<double> row_variances(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        // Calculate row mean
        double row_mean = 0.0;
        for (int j = 0; j < cols; ++j) {
            row_mean += M[i][j];
        }
        row_mean /= cols;

        // Center and normalize the row, and calculate variance
        for (int j = 0; j < cols; ++j) {
            N[i][j] = M[i][j] - row_mean;
            row_variances[i] += N[i][j] * N[i][j];
        }
    }

    // Normalize each row by its variance
    for (int i = 0; i < rows; ++i) {
        double norm = sqrt(row_variances[i]);
        for (int j = 0; j < cols; ++j) {
            N[i][j] /= norm;
        }
        row_variances[i] /= (cols-1); // Not the same
    }
    swap(N, M);
    return row_variances;
}
double center_normalize_vec(vector<double>& row) {
    int cols = row.size();

    // Calculate the mean of the row
    double row_mean = 0.0;
    for (int j = 0; j < cols; ++j) {
        row_mean += row[j];
    }
    row_mean /= cols;

    // Center and normalize the row, and calculate variance
    double row_variance = 0.0;
    for (int j = 0; j < cols; ++j) {
        row[j] -= row_mean;
        row_variance += row[j] * row[j];
    }

    // Normalize the row by its variance
    double norm = sqrt(row_variance);
    for (int j = 0; j < cols; ++j) {
        row[j] /= norm;
    }
    row_variance /= (cols - 1);

    return row_variance;
}