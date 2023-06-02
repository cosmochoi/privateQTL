#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

std::vector<double> CSVtoVector(std::string filename);
std::vector<uint64_t> ScaleVector(std::vector<double> &v, int k);
uint32_t nearestPowerOf2(int N);
#endif
