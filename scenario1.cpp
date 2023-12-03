
#include "mpc.h"
#include <cmath>
#include "utils.h"
#include <chrono>
#include "input.h"
#include <sys/resource.h>
struct rusage r_usage;

void bitdecompose(vector<uint32_t> &secrets, vector<BitVector> &bitInput)
{
    // vector<BitVector> bitInput;
    for (int i = 0; i < secrets.size(); i++)
    {
        string sinput = bitset<32>(secrets[i]).to_string();
        BitVector decomposed(sinput);
        bitInput.push_back(decomposed);
    }
}
double calculateVariance(const vector<double>& values) {
    double sum = 0.0;
    double mean = 0.0;
    size_t n = values.size();

    // Calculate the mean
    for (const double& value : values) {
        sum += value;
    }
    mean = sum / n;

    // Calculate the sum of squared differences
    double sumSquaredDiff = 0.0;
    for (const double& value : values) {
        double diff = value - mean;
        sumSquaredDiff += diff * diff;
    }

    // Calculate the variance
    double variance = sumSquaredDiff / n;
    return variance;
}
template <typename T>
void writeVectorToCSV(const vector<T>& data, string name)
{
    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name + ".csv";
    ofstream file(filename);
    if (file.is_open())
    {
        for (const T& value : data)
        {
            file << value << ",";
        }
        file.close();
        cout << string(name+ " vector successfully written to CSV file.") << endl;
    }
    else
    {
        cout << "Error opening the file." << endl;
    }
}
template <typename T>
void writematrixToTSV(const vector<vector<T>>& data, int startrow, int endrow, const string& name)
{
    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name + "_row" + std::to_string(startrow) + "_" + std::to_string(endrow) + ".tsv";
    ofstream file(filename);
    if (file.is_open())
    {
        for (int i = startrow; i <= endrow && i < data.size(); ++i)
        {
            for (const T& value : data[i])
            {
                file << value << "\t";
            }
            file << std::endl;
        }
        file.close();
        cout << string(name+" matrix successfully written to TSV file.") << endl;
    }
    else
    {
        cout << "Error opening the file." << endl;
    }
}
// vector<vector<double>> getMatrixFile(const string& filename, int startrow, int endrow, bool header, bool index) {
//     vector<vector<double>> rowsData;
//     ifstream data(filename);
//     string line;
//     int currentRow = 0;

//     if (header)// Skip the first row (header)
//         getline(data, line);

//     while (getline(data, line) && currentRow < endrow) {
//         if (currentRow >= startrow) { // Start reading from startrow
//             stringstream lineStream(line);
//             string cell;
//             int currentColumn = 0;

//             // Skip the first column
//             if (index)
//                 getline(lineStream, cell, '\t');

//             vector<double> rowVector;
//             while (getline(lineStream, cell, '\t')) {
//                 try {
//                     double entry = stod(cell);
//                     rowVector.push_back(entry);
//                 } catch (const exception& e) {
//                     cerr << "Exception caught: " << e.what() << endl;
//                 }
//                 currentColumn++;
//             }

//             rowsData.push_back(rowVector);
//         }

//         currentRow++;
//     }

//     // Close the file after reading
//     data.close();

//     return rowsData;
// }
vector<double> get_quantiles(vector<vector<double>>& phen_matrix, vector<vector<size_t>>& rank_matrix) {
    // vector<vector<size_t>> rank_matrix(phen_matrix[0].size(), vector<size_t>(phen_matrix.size()));
    for (size_t i = 0; i < phen_matrix[0].size(); i++) {
        vector<pair<double, size_t>> sorted_indices;
        for (size_t j = 0; j < phen_matrix.size(); j++) {
            sorted_indices.emplace_back(phen_matrix[j][i], j);
        }
        sort(sorted_indices.begin(), sorted_indices.end());
        for (size_t j = 0; j < phen_matrix.size(); j++) {
            rank_matrix[i][j] = sorted_indices[j].second;
        }
    }

    // vector<vector<double>> M = phen_matrix;
    // vector<vector<size_t>> Q = rank_matrix;
    size_t m = phen_matrix.size();
    size_t n = phen_matrix[0].size();

    vector<double> quantiles(m, 0.0);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            quantiles[j] += phen_matrix[rank_matrix[i][j]][i];
        }
    }

    return quantiles;
}

void sample_QN(vector<vector<double>>& phen_matrix, vector<vector<size_t>>& rank_matrix, vector<double>& total_quantiles) {
    // vector<vector<size_t>> rank_matrix(phen_matrix[0].size(), vector<size_t>(phen_matrix.size()));
    // for (size_t i = 0; i < phen_matrix[0].size(); i++) {
    //     vector<pair<double, size_t>> sorted_indices;
    //     for (size_t j = 0; j < phen_matrix.size(); j++) {
    //         sorted_indices.emplace_back(phen_matrix[j][i], j);
    //     }
    //     sort(sorted_indices.begin(), sorted_indices.end());
    //     for (size_t j = 0; j < phen_matrix.size(); j++) {
    //         rank_matrix[i][j] = sorted_indices[j].second;
    //     }
    // }

    // vector<vector<double>> M = phen_matrix;
    // vector<vector<size_t>> Q = rank_matrix;
    size_t m = phen_matrix.size();
    size_t n = phen_matrix[0].size();

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            phen_matrix[rank_matrix[i][j]][i] = total_quantiles[j]/n;
        }
    }
    // return phen_matrix;
}

vector<vector<double>> deseq2_cpm(vector<vector<uint32_t>>& counts_df) {
    int numGenes = counts_df.size();
    int numSamples = counts_df[0].size();

    vector<double> colSums(numSamples, 0.0);
    for (int j = 0; j < numSamples; ++j) {
        for (int i = 0; i < numGenes; ++i) {
            colSums[j] += static_cast<double>(counts_df[i][j]);
        }
    }
    vector<vector<double>> cpm_df(numGenes, vector<double>(numSamples, 0.0));
    for (int j = 0; j < numSamples; ++j) {
        for (int i = 0; i < numGenes; ++i) {
            cpm_df[i][j] = counts_df[i][j] / colSums[j] * 1e6;
        }
    }

    return cpm_df;
}

double findMedian(const vector<vector<double>>& matrix, size_t col) {
    vector<double> columnData;
    for (const auto& row : matrix) {
        columnData.push_back(row[col]);
    }

    sort(columnData.begin(), columnData.end());

    size_t size = columnData.size();
    if (size % 2 == 0) {
        return (columnData[size / 2 - 1] + columnData[size / 2]) / 2.0;
    } else {
        return columnData[size / 2];
    }
}

void dataclient(string norm_method, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend,string zscorefile, int permut)
{
    IOService ios;
    Endpoint epsend1(ios, address1, sendport1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendport2, EpMode::Server);
    Endpoint epsend3(ios, address3, sendport3, EpMode::Server);
    Endpoint eprecv1(ios, address1, recvport1, EpMode::Client);
    Endpoint eprecv2(ios, address2, recvport2, EpMode::Client);
    Endpoint eprecv3(ios, address3, recvport3, EpMode::Client);
    Channel owner_p1 = epsend1.addChannel();
    Channel owner_p2 = epsend2.addChannel();
    Channel owner_p3 = epsend3.addChannel();
    Channel p1_owner = eprecv1.addChannel();
    Channel p2_owner = eprecv2.addChannel();
    Channel p3_owner = eprecv3.addChannel();
    cout << "Owner established channels with the computing parties.\n";
    NTL::SetSeed((NTL::conv<NTL::ZZ>((long)27))); // Seed change
    vector<vector<double>> pheno;
    uint32_t p = pow(2,31);
    ZZ_p::init(to_ZZ(p)); 
    cout << "P: " << p << endl;
    owner_p1.send(p);
    owner_p2.send(p);
    owner_p3.send(p);
    auto pre_start = chrono::high_resolution_clock::now();

    vector<vector<int32_t>> geno_scaled;
    vector<double> geno_var, pheno_var;
    
    string pheno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/bed_template.tsv";
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_genotype.tsv";
    string geno_matrix = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_genotype_1kGresidualized.tsv";
    string geno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_variant.tsv";
    string pheno_cov = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/covariates/PC_covariate_df_670.csv";
    vector<vector<double>> covariates = getCovariates(pheno_cov);
    // cout << covariates.size() << "," << covariates[0].size() << endl;
    Residualizer res(covariates);   
    cout << "Loading Genotype matrix..." << flush;
    auto loadgeno = chrono::high_resolution_clock::now();
    prepareInput testinput(pheno_pos, geno_matrix,geno_pos,1000000);
    auto loadend = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = loadend - loadgeno;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    cout << durationInminutes << " minutes" << endl;

    for (int r=rowstart; r<rowend; r++)
    {
        cout <<"1"<< endl;
        vector<string> cisVariants;
        vector<double> std_ratio;
        string geneID;
        string pheno_file = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_qn_invCDF.bed";
        // string pheno_file = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/mis_acc.qn.bed";
        vector<double> norm_pheno;
        read_bedfile_row(norm_pheno,geneID, pheno_file, r,4,true);
        cout <<"2"<< endl;
        vector<double> pheno_res = res.transform(norm_pheno);
        // print_vector(pheno_res);
        double bed_var = center_normalize_vec(norm_pheno);
        vector<int32_t> pheno = ScaleVector_signed(norm_pheno, pow(10,5)); // integer version of secret

        // cout << string("bed file phenotype var: "+to_string(bed_var)+"\n");
        vector<uint32_t> range;
        string chromosome = testinput.getCisRange(geneID,range);
        // cout << string("gene "+geneID[r]+"/chr "+ chromosome+ "/start "+ to_string(range[0])+"/end "+ to_string(range[1])+"\n");
        vector<vector<double>> slicedgeno = testinput.sliceGeno(range, chromosome, -1,cisVariants);
        vector<vector<double>> geno_res = res.transform(slicedgeno);
        geno_var = center_normalize(geno_res);
        // writeVectorToCSV(geno_var, string("row"+to_string(r)+"_geno_var"));
        // cout << string("first geno var: "+to_string(geno_var[0])+"\n");
        geno_scaled = ScaleVector(geno_res, pow(10,4));
        
        string original = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/zscores.txt");
        vector<double> pre_centered = CSVtoVector(original);
        // cout << "phenotype variance: " << doublevariance(pre_centered, doublemean(pre_centered)) << endl;
        double pheno_var =doublevariance(pre_centered, doublemean(pre_centered));//0.9782648500530864;// 0.9596748533543238; //TODO::Don't put manual numbers here 
        // cout << "Phenotype_var: " << pheno_var << endl;
        // writeVectorToCSV(geno_var, "pQTL_geno_var");
        for (size_t j = 0; j < geno_var.size(); ++j) {
            std_ratio.push_back(sqrt(pheno_var / geno_var[j]));
        }
        // writeVectorToCSV(std_ratio, "pQTL_std_ratio");
        
        uint32_t inv = PowerMod(3, -1, p);
        vector<vector<uint32_t>> shares;
        for (int i = 0; i < 3; i++) {
            shares.push_back(vector<uint32_t>());
        }
        // making the shares: for each bit, loop through secrets
        for (int i=0; i<pheno.size(); i++)
        {
            ZZ_p x1 = random_ZZ_p();
            ZZ_p x2 = random_ZZ_p();
            ZZ_p x3;
            if (pheno[i]<0)
                x3 = (conv<ZZ_p>(pheno[i]+p)) - x1 - x2;
            else
                x3 = conv<ZZ_p>(pheno[i]) - x1 - x2;

            uint32_t send_x1 = conv<uint32_t>(x1);
            uint32_t send_x2 = conv<uint32_t>(x2);
            uint32_t send_x3 = conv<uint32_t>(x3);
            shares[0].push_back(send_x1);
            shares[0].push_back(send_x2);
            shares[1].push_back(send_x2);
            shares[1].push_back(send_x3);                    
            shares[2].push_back(send_x3);
            shares[2].push_back(send_x1);
        }
        owner_p1.send(shares[0]);
        owner_p2.send(shares[1]);
        owner_p3.send(shares[2]);
        
        // cout << "Secrets shared " << std::time(nullptr) << endl;
        cout << "Sent secret shared row values to parties.\n";
    
        vector<vector<uint32_t>> genoshares;
        for (int i = 0; i < 3; i++) {
            genoshares.push_back(vector<uint32_t>());
        }
        // sending geno shares
        for (int i=0; i< geno_scaled.size(); i++)
        {
            for (int j=0; j< geno_scaled[0].size(); j++)
            {
                ZZ_p i1 = random_ZZ_p();
                ZZ_p i2 = random_ZZ_p();
                ZZ_p i3;
                if (geno_scaled[i][j]<0)
                    i3 = (conv<ZZ_p>(geno_scaled[i][j]+p)) - i1 - i2;
                else
                    i3 = conv<ZZ_p>(geno_scaled[i][j]) - i1 - i2;

                uint32_t send_i1 = conv<uint32_t>(i1);
                uint32_t send_i2 = conv<uint32_t>(i2);
                uint32_t send_i3 = conv<uint32_t>(i3);
                genoshares[0].push_back(send_i1);
                genoshares[0].push_back(send_i2);
                genoshares[1].push_back(send_i2);
                genoshares[1].push_back(send_i3);
                genoshares[2].push_back(send_i3);
                genoshares[2].push_back(send_i1);
            }
        }
        uint32_t testg_row = geno_scaled.size();
        uint32_t testg_col = geno_scaled[0].size();
        vector<uint32_t> genoshape = {
            static_cast<uint32_t>(testg_row),
            static_cast<uint32_t>(testg_col),
        };

        owner_p1.send(genoshape);
        owner_p2.send(genoshape);
        owner_p3.send(genoshape);
        owner_p1.send(genoshares[0]);
        owner_p2.send(genoshares[1]);
        owner_p3.send(genoshares[2]);

        // sending std_ratio in plaintext to party 1
        owner_p1.send(std_ratio);
        string serializedvariants=string(geneID+";");
        for (const std::string& str : cisVariants) {
            serializedvariants += str + ";"; // Use a suitable delimiter
        }
        int estimatedSize = static_cast<int>(serializedvariants.size());
        // int safetyMargin = 128; // Adjust this based on your needs
        int bufferSize = estimatedSize;
        // Send the serialized data over the channel
        owner_p1.send(bufferSize);
        owner_p1.send(serializedvariants.data(), serializedvariants.size());
        // cout << "Geno shared " << std::time(nullptr) << endl;
        cout << "Sent secret shared geno to parties.\n";
        vector<double> agg_stat;
        p1_owner.recv(agg_stat);
    }

    owner_p1.close();
    owner_p2.close();
    owner_p3.close();
    p1_owner.close();
    p2_owner.close();
    p3_owner.close();
    epsend1.stop();
    epsend2.stop();
    epsend3.stop();
    eprecv1.stop();
    eprecv2.stop();
    eprecv3.stop();
    ios.stop();
}

vector<ZZ_p> convert(vector<int> input)
{
    vector<ZZ_p> output;
    for (int i=0; i<input.size(); i++)
    {
        output.push_back(conv<ZZ_p>(input[i]));
    }
    return output;
}
void runMPC(string norm_method, int pid,  string ownerIP, int ownerPort, int toOwnerPort,  string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2,int rowstart, int rowend, int permut, Logger& cisLogger, Logger& nominalLogger) //, atomic<int>& readyCounter, mutex& mtx, condition_variable& cv)
{    
    mpc testmpc;
    // std::promise<void> promise;
    // std::future<void> future = promise.get_future();
    // vector<vector<double>> resultVectors;
    testmpc.initialize(pid, ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    // if (norm_method == "deseq2")
    // {
    //     testmpc.receivePheno();
    // // vector<ZZ_p> random = testmpc.Frand(4);
    //     testmpc.logRatio();
    // }
    auto invcdf_start = chrono::high_resolution_clock::now();
    for (int row=rowstart; row<rowend; row++)
    {
        // cout << string("pid"+to_string(pid)+"_row"+to_string(row)+"\n");
        // vector<double> row_result;
        // cout << "start: " << std::time(nullptr) << endl;
        // testmpc.ready();
        // cout << "ready: " << std::time(nullptr) << endl;
        testmpc.permutPheno(permut);
        // cout << "receiveSecrets: " << std::time(nullptr) << endl;
        testmpc.receiveGeno();
        // cout << "receiveGeno: " << std::time(nullptr) << endl;
        // testmpc.genperm(row, numCol, norm_method, permut);
        // cout << "genperm: " << std::time(nullptr) << endl;
        testmpc.calc_corr(cisLogger, nominalLogger);
        // cout << "testMatrix: " << std::time(nullptr) << endl;
        // resultVectors.emplace_back(row_result);
        testmpc.clearVectors();
    }
    auto invcdf_end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = invcdf_end - invcdf_start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    // if (pid == 0)
    //     cout << string("InverseCDF Execution time: " + to_string(durationInminutes) + " minutes\n") << endl;
    testmpc.close();
    // swap(resultVector, resultVectors);
}
void setThreadAffinity(std::thread& thread, int coreId)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(coreId, &cpuset);
    pthread_setaffinity_np(thread.native_handle(), sizeof(cpu_set_t), &cpuset);
}
bool isPortOpen(int port) {
    int sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) {
        perror("socket");
        return false;
    }

    sockaddr_in addr{};
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    addr.sin_addr.s_addr = htonl(INADDR_ANY);

    int result = bind(sockfd, (struct sockaddr*)&addr, sizeof(addr));
    if (result < 0) {
        close(sockfd);
        return false;
    }

    close(sockfd);
    return true;
}
void startMPCset(string norm_method, vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, string zscorefile, int CPU_core, int permut, Logger& cislogger, Logger& nominalLogger)
{
    // cout << "startMpcset Start\n";
    vector<thread> threads;
    // vector<vector<double>> resultVectors1,resultVectors2,resultVectors3;
    thread dataClientThread(dataclient, norm_method, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], 
    availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow, zscorefile, permut);
    setThreadAffinity(dataClientThread,CPU_core+0);
    threads.emplace_back(move(dataClientThread));

    // int startidx1 = startidx;
    thread runMPC1([=, &cislogger, &nominalLogger]() {
        runMPC(norm_method, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], 
        startrow, endrow, permut,cislogger,nominalLogger);
    });
    // thread runMPC1(runMPC, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], startrow, endrow, resultVectors1);
    setThreadAffinity(runMPC1,CPU_core+1);
    threads.emplace_back(move(runMPC1));

    // // int startidx2 = startidx;
    thread runMPC2([=, &cislogger,&nominalLogger]() {
        runMPC(norm_method, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], 
        startrow, endrow,permut,cislogger,nominalLogger);
    });
    // thread runMPC2(runMPC, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], startrow, endrow, resultVectors2);
    setThreadAffinity(runMPC2,CPU_core+2);
    threads.emplace_back(move(runMPC2));

    thread runMPC3([=, &cislogger,&nominalLogger]() {
        runMPC(norm_method, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], 
        startrow, endrow, permut,cislogger,nominalLogger);
    });
    // thread runMPC3(runMPC, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], startrow, endrow, resultVectors3);
    setThreadAffinity(runMPC3,CPU_core+2);
    threads.emplace_back(move(runMPC3));
    for (auto& thread : threads) {
        thread.join();
    }
    // swap(finalVec, resultVectors1);
}

double cal_accuracy(vector<vector<double>>& predicted, string& filename, int startrow, int endrow)
{
    vector<vector<double>> actual = getMatrixFile(filename,startrow, endrow,true,true);
    // vector<double> actual = CSVtoVector(filename);
    // cout << string("actual matrix shape:" + to_string(actual.size())+ ","+to_string(actual[0].size())+"\n");
    // cout << string("predicted matrix shape:" + to_string(predicted.size())+ ","+to_string(predicted[0].size())+"\n");
    if (predicted.size() != actual.size()) {
        cout << "predicted size: " + to_string(predicted.size()) << endl;
        cout << "actual size: " + to_string(actual.size()) << endl;
        // throw runtime_error("Vector sizes do not match");
    }

    double sumSquaredDiff = 0.0;
    double max=0.0, min=0.0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        for (size_t j=0; j<predicted[0].size(); ++j)
        {
            double diff = predicted[i][j] - actual[i][j];
            sumSquaredDiff += diff * diff;
            if (actual[i][j]>=max)
                max=actual[i][j];
            if (actual[i][j]<=min)
                min=actual[i][j];
        }
        
    }

    double mse = sumSquaredDiff / (predicted.size()*predicted[0].size());
    double rmse = sqrt(mse);
    // cout << string("RMSE: " + to_string(rmse)) << endl;
    double normalized = max-min;
    // auto max = max_element(actual.begin(), actual.end());
    // auto min = min_element(actual.begin(), actual.end());
    // double normalized = *max - *min;
    return 1.0 - rmse/normalized;

}
int main(int argc, char* argv[])
{
    if (argc < 8)
    {
        cout << "Please provide at least four arg: low, mid, high, permutation, zscorefile name, norm method, cislog, nominallog.\n";
        return 1;
    }
    int lo_row = stoi(argv[1]);
    int mid_row = stoi(argv[2]);
    int hi_row = stoi(argv[3]);
    int permut = stoi(argv[4]);
    string zscorefile = argv[5];
    string norm_method = argv[6];
    string cis_log = argv[7];
    string nominal_log = argv[8];
    int startingPort = 12300;
    int assignedPorts=0;
    vector<int> openPorts;
    while (assignedPorts < 24)
    {
        int port = startingPort;
        if(isPortOpen(port))
        {
            openPorts.push_back(port);
            assignedPorts++;
        }
        startingPort++;
    }
    // print_vector(openPorts);
    vector<string> address{"localhost","localhost","localhost","localhost"};

    auto start = chrono::high_resolution_clock::now();
    // cout << "Start time: " << std::time(nullptr) << endl;
    int numCores = thread::hardware_concurrency();
    // mutex mtx;
    // condition_variable cv;
    // std::atomic<int> readyCounter(0);
    vector<thread> threads;
    // vector<vector<double>> resultVec1,resultVec2,resultVec3,resultVec4;
    string nominal = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/output/privateQTL_scenario1_"+nominal_log);
    string cis = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/output/privateQTL_scenario1_" +cis_log);
    Logger nominallogger(nominal+"_"+to_string(lo_row)+"_"+to_string(mid_row)+".tsv"), cislogger(cis+"_"+to_string(lo_row)+"_"+to_string(mid_row)+".tsv");
    Logger nominallogger2(nominal+"_"+to_string(mid_row)+"_"+to_string(hi_row)+".tsv"),cislogger2(cis+"_"+to_string(mid_row)+"_"+to_string(hi_row)+".tsv");
    nominallogger.log(string("phenID\tvarID\tvarIdx\tdof\tr_nom\tr2_nom\ttstat\tpval\tslope\tslope_se"));
    cislogger.log(string("phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\tpval_true_df\tr_nom\tr2_nom\ttstat\tpval_nominal\tslope\tslope_se\tpval_perm\tpval_beta"));
    nominallogger2.log(string("phenID\tvarID\tvarIdx\tdof\tr_nom\tr2_nom\ttstat\tpval\tslope\tslope_se"));
    cislogger2.log(string("phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\tpval_true_df\tr_nom\tr2_nom\ttstat\tpval_nominal\tslope\tslope_se\tpval_perm\tpval_beta"));
    thread thread1([&]() {
        startMPCset(norm_method, openPorts, 0, address, lo_row, mid_row, zscorefile, 0,permut,cislogger,nominallogger);
    });
    thread thread2([&]() {
        startMPCset(norm_method, openPorts, 12, address, mid_row, hi_row,zscorefile, 0, permut,cislogger2,nominallogger2);
    });

    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));

    for (auto& thread : threads) {
        thread.join();
    }

    // resultVec1.insert(resultVec1.end(), resultVec2.begin(), resultVec2.end());
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> totalduration = end - start;
    double totaldurationInSeconds = totalduration.count();
    double totaldurationInminutes = totaldurationInSeconds/60.0;
     if (getrusage(RUSAGE_SELF, &r_usage) == 0) {
        std::cout << "Memory usage: " << r_usage.ru_maxrss << " KB" << std::endl;
    } else {
        std::cerr << "Failed to get resource usage." << std::endl;
    }
    cout << "Execution time: " << totaldurationInminutes << " minutes" << endl;

    return 0;
}