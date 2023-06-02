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

std::vector<double> getRowFromMatrixFile(const std::string& filename, int rowIndex) {
    std::vector<double> rowVector;

    std::ifstream data(filename);
    std::string line;
    int currentRow = 0;
    while (std::getline(data, line)) {
        if (currentRow == 0) {
            // Skip the first row (header)
            currentRow++;
            continue;
        }
        if (currentRow-1 == rowIndex) {
            std::stringstream lineStream(line);
            std::string cell;

            while (std::getline(lineStream, cell, '\t')) 
            {
                try {
                double entry = std::stod(cell);
                rowVector.push_back(entry);
                }
                catch (const std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }
            }
        }
        currentRow++;
    }
    return rowVector;
}

vector<uint64_t> ScaleVector(vector<double> &v, int k)
{
    vector<uint64_t> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = v[i] * k;
    return intvec;
}


uint32_t nearestPowerOf2(int N)
{
    uint32_t a = log2(N);
 
    // if (pow(2, a) == N)
    //     return N;
 
    return pow(2, a + 1);
}



