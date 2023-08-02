#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <sstream>
using namespace std;

struct GeneData {
    std::string chr;
    uint32_t start;
    uint32_t end;
};

class prepareInput
{
public:
    uint32_t sampleN;
    uint32_t phenoN;
    uint32_t genoN;
    uint32_t ciswindow;
    unordered_map<string, GeneData> genePos;
    vector<uint32_t> snpPos;
    vector<string> snpChr;
    vector<string> snpIDs;
    vector<vector<int32_t>> geno;
    prepareInput(string& pheno_pos, string& geno_matrix, string& geno_pos, uint32_t window);
    string getCisRange(string geneID, vector<uint32_t>& positions);
    vector<uint32_t> getSNPrange(uint32_t start, uint32_t end, string chr);
    vector<vector<int32_t>> sliceGeno(vector<uint32_t> positions, string& chr);

private:
    // Add any private member functions or variables if needed
};

#endif 
