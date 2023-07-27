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

vector<vector<double>> getTPMFromMatrixFile(const string& filename, int numRow, int numCol) {
    vector<vector<double>> rowsData;
    ifstream data(filename);
    string line;
    int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line) && currentRow < numRow) {
        stringstream lineStream(line);
        string cell;
        int currentColumn = 0;

        // Skip the first two columns
        getline(lineStream, cell, '\t');
        getline(lineStream, cell, '\t');

        vector<double> rowVector;
        while (currentColumn < numCol && getline(lineStream, cell, '\t')) {
            try {
                double entry = stod(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cerr << "Exception caught: " << e.what() << endl;
            }
            currentColumn++;
        }

        rowsData.push_back(rowVector);
        currentRow++;
    }

    // Close the file after reading
    data.close();

    return rowsData;
}
vector<vector<uint32_t>> getCountFromMatrixFile(const string& filename, int numRow, int numCol) {
    vector<vector<uint32_t>> rowsData;
    ifstream data(filename);
    string line;
    int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line) && currentRow < numRow) {
        stringstream lineStream(line);
        string cell;
        int currentColumn = 0;

        // Skip the first two columns
        getline(lineStream, cell, '\t');
        getline(lineStream, cell, '\t');

        vector<uint32_t> rowVector;
        while (currentColumn < numCol && getline(lineStream, cell, '\t')) {
            try {
                uint32_t entry = stoi(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cerr << "Exception caught: " << e.what() << endl;
            }
            currentColumn++;
        }

        rowsData.push_back(rowVector);
        currentRow++;
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
vector<vector<uint32_t>> ScaleVector(vector<vector<double>>& v, int k)
{
    vector<vector<uint32_t>> intvec(v.size(), vector<uint32_t>(v[0].size()));
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


vector<double> UnscaleVector_signed(vector<uint32_t> &v, int k)
{
    vector<double> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = static_cast<double>(v[i]) / static_cast<double>(k);
    return intvec;
}
vector<double> ShiftVector(vector<double>& vec, double number) 
{
    vector<double> shifted(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        shifted[i] = vec[i] + number;
        // shifted[i] = static_cast<uint32_t>(shiftedValue);
    }
    return shifted;
}

vector<double> UnshiftVector(vector<double>& shifted, double number)
{
    vector<double> unshifted(shifted.size());
    for (size_t i = 0; i < shifted.size(); ++i) {
        unshifted[i] = shifted[i] - number;
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



