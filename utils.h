#ifndef UTILS_H
#define UTILS_H
#define MATHLIB_STANDALONE
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <Rmath.h>
#include <Eigen/Dense>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#define BETA_SHAPE1_MIN 0.1
#define BETA_SHAPE2_MIN 1
#define BETA_SHAPE1_MAX 10
#define BETA_SHAPE2_MAX 1000000	
using namespace std;
using namespace Eigen;
using namespace NTL;
vector<double> CSVtoVector(string filename);
vector<double> getRowFromMatrixFile(string& filename, int rowIndex);
void read_bedfile_row(vector<double>& rowData, string& geneID, const string& filename, int row, int skipcols, bool header);
vector<vector<double>> getTPMFromMatrixFile(const string& filename, vector<string>& geneID);
vector<vector<uint32_t>> getCountFromMatrixFile(const string& filename, vector<string>& geneID);
vector<uint32_t> ScaleVector(vector<double> &v, int k);
vector<vector<int32_t>> ScaleVector(vector<vector<double>>& v, int k);
vector<int32_t> ScaleVector_signed(vector<double> &v, int k);
vector<double> UnscaleVector_signed(vector<int32_t> &v, int k);
vector<uint32_t> ShiftVector(vector<int32_t>& vec, uint32_t number); 
vector<vector<uint32_t>> ShiftVector(vector<vector<int32_t>>& vec, uint32_t number);
vector<int32_t> UnshiftVector(vector<uint32_t>& shifted, uint32_t number);
uint32_t nearestPowerOf2(int N);
int mleBeta(vector < double > & pval, double & beta_shape1, double & beta_shape2);
double doublemean(vector < double > & X);
double doublevariance(vector < double > & X, double mean);
double degreeOfFreedom(const gsl_vector *v, void *params);
int learnDF(vector < double > & corr, double & df);
double getTstat2(double corr, double df);
double getPvalueFromTstat2(double tstat2, double df);
double getPvalue(double corr, double df);
double getSlope(double nominal_correlation, double psd, double gsd);
vector<double> center_normalize(vector<vector<double>>& M);
double center_normalize_vec(vector<double>& M);
// void ZZtoEigen(vector<vector<ZZ_p>>& v, MatrixXi& dest1, MatrixXi& dest2);
// void EigentoZZ(vector<uint32_t>& share1, vector<uint32_t>& share2, vector<vector<ZZ_p>>& dest);
#endif
