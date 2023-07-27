#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
using namespace std;

vector<double> CSVtoVector(string filename);
// vector<double> getRowFromMatrixFile(string& filename, int rowIndex, int numCol);

vector<double> getRowFromMatrixFile(string& filename, int rowIndex, int numCol);
vector<vector<double>> getTPMFromMatrixFile(const string& filename, int numRow, int numCol);
vector<vector<uint32_t>> getCountFromMatrixFile(const string& filename, int numRow, int numCol);
vector<uint32_t> ScaleVector(vector<double> &v, int k);
vector<vector<uint32_t>> ScaleVector(vector<vector<double>>& v, int k);
vector<int32_t> ScaleVector_signed(vector<double> &v, int k);
vector<double> UnscaleVector_signed(vector<uint32_t> &v, int k);
vector<double> ShiftVector(vector<double>& vec, double number); 
vector<double> UnshiftVector(vector<double>& shifted, double number);
uint32_t nearestPowerOf2(int N);
#endif
