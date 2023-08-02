
#include "mpc.h"
#include <cmath>
#include "utils.h"
#include <chrono>
#include "input.h"

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
vector<vector<double>> getMatrixFile(const string& filename, int startrow, int endrow, int numCol) {
    vector<vector<double>> rowsData;
    ifstream data(filename);
    string line;
    int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line) && currentRow < endrow) {
        if (currentRow >= startrow) { // Start reading from startrow
            stringstream lineStream(line);
            string cell;
            int currentColumn = 0;

            // Skip the first column
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
        }

        currentRow++;
    }

    // Close the file after reading
    data.close();

    return rowsData;
}
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

void dataclient(string norm_method, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend,int numCol, string zscorefile)
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
    uint32_t p = pow(2,20);
    ZZ_p::init(to_ZZ(p)); 
    owner_p1.send(p);
    owner_p2.send(p);
    owner_p3.send(p);
    auto pre_start = chrono::high_resolution_clock::now();
    vector<string> geneID;
    if (norm_method == "qn")
    {
        // string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/GTEx_Analysis_gene_tpm_matched_filtered.gct";
        string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_tpm_matched_filtered.tsv";
        pheno = getTPMFromMatrixFile(matched, geneID, 25949, 15253);
        // print_vector(testss);
        vector<vector<size_t>> rank(pheno[0].size(), vector<size_t>(pheno.size()));
        vector<double> quantiles =get_quantiles(pheno, rank);
        sample_QN(pheno, rank, quantiles);
        cout << "Quantile normalization completed.\n";
    }
    else if (norm_method == "deseq2")
    {
        vector<vector<uint32_t>> testg = {{1,2,3,4},{5,6,1,2},{3,4,5,6}};
        // string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/GTEx_Analysis_gene_counts_matched_filtered.gct";
        string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_reads_matched_filtered.tsv";
        vector<vector<uint32_t>> testp = getCountFromMatrixFile(matched,geneID, 25949,15253);
        vector<vector<double>> cpm_df = deseq2_cpm(testp);
        // print_vector(cpm_df);
        // vector<vector<uint32_t>> pheno = ScaleVector(cpm_df, pow(10,3));
        // vector<vector<uint32_t>> testp = {{1,0,3,3},{4,5,6,6},{7,8,9,6},{10,11,12,4}};
        vector<vector<uint32_t>> genoshares, phenoshares;
        for (int i = 0; i < 3; i++) {
            genoshares.push_back(vector<uint32_t>());
            phenoshares.push_back(vector<uint32_t>());
        }
        // sending geno pheno shares
        for (int i=0; i< testg.size(); i++)
        {
            for (int j=0; j< testg[0].size(); j++)
            {
                ZZ_p i1 = random_ZZ_p();
                ZZ_p i2 = random_ZZ_p();
                ZZ_p i3 = conv<ZZ_p>(testg[i][j]) - i1 - i2;
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
        uint32_t testp_row = cpm_df.size();
        uint32_t testp_col = cpm_df[0].size();
        uint32_t testg_row = testg.size();
        uint32_t testg_col = testg[0].size();
        vector<int> exclude;
        for (int i=0; i<cpm_df.size(); i++)
        {
            vector<uint32_t> share1_row, share2_row, share3_row;
            for (int j=0; j< testp[0].size(); j++)
            {
                if (cpm_df[i][j]<=0)
                {
                    share1_row.clear();
                    share2_row.clear();
                    share3_row.clear();
                    testp_row--;
                    exclude.push_back(i);
                    break;
                }
                double logcount = log(cpm_df[i][j]);
                uint32_t scaledcount = logcount*pow(10,3);
                // cout << scaledcount << endl;
                ZZ_p i1 = random_ZZ_p();
                ZZ_p i2 = random_ZZ_p();
                ZZ_p i3 = conv<ZZ_p>(scaledcount) - i1 - i2;
                uint32_t send_i1 = conv<uint32_t>(i1);
                uint32_t send_i2 = conv<uint32_t>(i2);
                uint32_t send_i3 = conv<uint32_t>(i3);
                share1_row.push_back(send_i1);
                share1_row.push_back(send_i2);
                share2_row.push_back(send_i2);
                share2_row.push_back(send_i3);
                share3_row.push_back(send_i3);
                share3_row.push_back(send_i1);
            }
            phenoshares[0].insert(phenoshares[0].end(), share1_row.begin(), share1_row.end());
            phenoshares[1].insert(phenoshares[1].end(), share2_row.begin(), share2_row.end());
            phenoshares[2].insert(phenoshares[2].end(), share3_row.begin(), share3_row.end());
        }
        vector<uint32_t> matshape = {
            static_cast<uint32_t>(testg.size()),
            static_cast<uint32_t>(testg[0].size()),
            static_cast<uint32_t>(testp_row),
            static_cast<uint32_t>(testp_col)
        };

        owner_p1.send(matshape);
        owner_p2.send(matshape);
        owner_p3.send(matshape);
        owner_p1.send(genoshares[0]);
        owner_p2.send(genoshares[1]);
        owner_p3.send(genoshares[2]);
        owner_p1.send(phenoshares[0]);
        owner_p2.send(phenoshares[1]);
        owner_p3.send(phenoshares[2]);
        cout << "Sent secret shared geno pheno to parties.\n";
        vector<double> ref1, ref2, ref3;
        // cout << testp.size()-exclude.size() << endl;
        vector<vector<double>> finalratio(testp[0].size(), vector<double>(testp.size()-exclude.size()));
        p1_owner.recv(ref1);
        p2_owner.recv(ref2);
        p3_owner.recv(ref3);
        // print_vector(ref1);
        int skipped = 0;
        for (int i=0; i<cpm_df.size(); i++)
        {
            auto it = find(exclude.begin(), exclude.end(), i);
            if (it != exclude.end())
            {
                // cout <<string("Gene "+to_string(i)+" excluded;\n");
                skipped++;
                continue;
            }
            for(int j=0; j<cpm_df[0].size();j++)
            {
                // double logcount = log(testp[i][j]);
                finalratio[j][i-skipped] = log(cpm_df[i][j]) - ref1[i-skipped]/pow(10,3);
            }
        }
        vector<double> sizefactors;
        for (int j=0; j<finalratio.size(); j++)
        {
            sort(finalratio[j].begin(), finalratio[j].end());
            size_t size = finalratio[j].size();
            if (size % 2 == 0) {
                double medval = (finalratio[j][size / 2 - 1] + finalratio[j][size / 2]) / 2.0;
                sizefactors.push_back(exp(medval));
            } else {
                sizefactors.push_back(exp(finalratio[j][size / 2]));
            }
        }
        for (int i=0; i<cpm_df.size(); i++)
        {
            for (int j=0; j<cpm_df[0].size(); j++)
            {
                cpm_df[i][j] = cpm_df[i][j]/sizefactors[j];
            }
        }
        swap(cpm_df, pheno);
        // print_vector(pheno);
    }
    else 
    {
        throw invalid_argument("Please choose normalization method between qn and deseq2.\n");
    }
    auto pre_end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = pre_end - pre_start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    cout << "Normalization Execution time: " << durationInminutes << " minutes" << endl;
    for (int r=rowstart; r<rowend; r++)
    {
        int ready1, ready2, ready3;
        p1_owner.recv(ready1);
        p2_owner.recv(ready2);
        p3_owner.recv(ready3);
        if ((ready1 == 1) && (ready2 ==1) && (ready3==1))
        {
            cout << "Client will send secrets.\n";
            vector<BitVector> bitInput;
            // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/samplewise_QN.tsv";
            // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/correct_sampleQN.tsv";
            // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/sample_thres125_fp1.tsv";
            // vector<double> input = getRowFromMatrixFile(filename,r, numCol);

            // string filename ="/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/testdata2.txt";
            // vector<double> input = pheno[r][:numCol];
            vector<double> input(pheno[r].begin(), pheno[r].begin() + numCol);
            cout << string("gene "+ geneID[r]+"\n");
            string pheno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/bed_template.tsv";
            string geno_matrix = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_centered.tsv";
            string geno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_variantdf.tsv";
            prepareInput testinput(pheno_pos, geno_matrix,geno_pos,1000000);
            vector<uint32_t> range;
            string chromosome = testinput.getCisRange(geneID[r],range);
            cout << string("gene "+geneID[r]+"/chr "+ chromosome+ "/start "+ to_string(range[0])+"/end "+ to_string(range[1])+"\n");
            vector<vector<int32_t>> slicedgeno = testinput.sliceGeno(range, chromosome);
            cout << string("genotype shape: "+ to_string(slicedgeno.size())+"/ "+to_string(slicedgeno[0].size())+"\n");
            vector<uint32_t> secrets = ScaleVector(input, pow(10,3)); // integer version of secret
            auto smax = max_element(secrets.begin(), secrets.end());
            // uint32_t absmax = abs(*smax);
            uint32_t maxsecret = nearestPowerOf2(*smax);
            bitdecompose(secrets, bitInput);

            
            /// ZSCORE FILE
            string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/"+zscorefile+".txt");
            // string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/" + zscorefile + ".txt");
            vector<double> zscore_input = CSVtoVector(zscore_filename);
            auto min = min_element(zscore_input.begin(), zscore_input.end());
            double shiftsize = abs(*min);
            // cout << "shift size: " << shiftsize << endl;
            vector<double> zscores_shifted = ShiftVector(zscore_input, shiftsize);
            vector<uint32_t> zscores_scaled = ScaleVector(zscores_shifted, pow(10,2));
            auto max = max_element(zscores_scaled.begin(), zscores_scaled.end());
            uint32_t maxzscore = *max;
            uint32_t maxzsecret = nearestPowerOf2(maxzscore);
            uint32_t secretsize = nearestPowerOf2(bitInput.size());
            // uint32_t p = 32768;
            // uint32_t p = (maxsecret > maxzsecret) ? ((maxsecret > secretsize) ? maxsecret : secretsize) : ((maxzsecret > secretsize) ? maxzsecret : secretsize);
            cout << "P: " << p << endl;
            // ZZ_p::init(to_ZZ(p)); 
            uint32_t inv = PowerMod(3, -1, p);
            vector<uint32_t> vectorsize{(uint32_t) bitInput.size(), (uint32_t) bitInput[0].size(), p, inv};
            owner_p1.send(vectorsize);
            owner_p2.send(vectorsize);
            owner_p3.send(vectorsize);
            cout << "vector size sent.\n";
            // Send the original input for verification
            // owner_p1.asyncSend(input);
            // owner_p2.asyncSend(input);
            // owner_p3.asyncSend(input);
            vector<vector<uint32_t>> shares;
            for (int i = 0; i < 3; i++) {
                shares.push_back(vector<uint32_t>());
            }
            // making the shares: for each bit, loop through secrets
            for (int j=bitInput[0].size()-1; j >=0; j--)
            {
                for (int i=0; i<bitInput.size(); i++)
                {
                    ZZ_p x1 = random_ZZ_p();
                    ZZ_p x2 = random_ZZ_p();
                    ZZ_p x3 = to_ZZ_p(bitInput[i][j]) - x1 - x2;
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
            }
            owner_p1.send(shares[0]);
            owner_p2.send(shares[1]);
            owner_p3.send(shares[2]);
            int prelim1, prelim2, prelim3;
            p1_owner.recv(prelim1);
            p2_owner.recv(prelim2);
            p3_owner.recv(prelim3);
            if ((prelim1 == 1) && (prelim1 ==1) && (prelim1==1))
            {
                cout << "Client will share zscore and identity.\n";
                // string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/"+zscorefile+".txt");
                // string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/" + zscorefile + ".txt");
                // vector<double> zscore_input = CSVtoVector(zscore_filename);
                // print_vector(zscore_input);
                // vector<int32_t> zscores_int = ScaleVector_signed(zscore_input, pow(10,5));
                // auto min = min_element(zscore_input.begin(), zscore_input.end());
                // double shiftsize = abs(*min);
                // cout << "shift size: " << shiftsize << endl;
                // vector<double> zscores_shifted = ShiftVector(zscore_input, shiftsize);
                // vector<uint32_t> zscores_scaled = ScaleVector(zscores_shifted, pow(10,2));
                // cout << "zscores size: " << zscores_scaled.size() << endl;
                vector<vector<uint32_t>> identity_shares;
                for (int i = 0; i < 3; i++) {
                    identity_shares.push_back(vector<uint32_t>());
                }
                //sending identity shares
                for (int j=0; j< bitInput.size(); j++)
                {
                    ZZ_p i1 = random_ZZ_p();
                    ZZ_p i2 = random_ZZ_p();
                    ZZ_p i3 = conv<ZZ_p>(j+1) - i1 - i2;
                    uint32_t send_i1 = conv<uint32_t>(i1);
                    uint32_t send_i2 = conv<uint32_t>(i2);
                    uint32_t send_i3 = conv<uint32_t>(i3);
                    identity_shares[0].push_back(send_i1);
                    identity_shares[0].push_back(send_i2);
                    identity_shares[1].push_back(send_i2);
                    identity_shares[1].push_back(send_i3);
                    identity_shares[2].push_back(send_i3);
                    identity_shares[2].push_back(send_i1);
                }
                owner_p1.send(identity_shares[0]);
                owner_p2.send(identity_shares[1]);
                owner_p3.send(identity_shares[2]);
                owner_p1.send(shiftsize);
                owner_p2.send(shiftsize);
                owner_p3.send(shiftsize);
                vector<vector<uint32_t>> zscore_shares;
                for (int i = 0; i < 3; i++) {
                    zscore_shares.push_back(vector<uint32_t>());
                }
                for (int i=0; i<zscores_scaled.size(); i++)
                {
                    ZZ_p z1 = random_ZZ_p();
                    ZZ_p z2 = random_ZZ_p();
                    ZZ_p z3 = to_ZZ_p(zscores_scaled[i]) - z1 - z2;
                    uint32_t send_z1 = conv<uint32_t>(z1);
                    uint32_t send_z2 = conv<uint32_t>(z2);
                    uint32_t send_z3 = conv<uint32_t>(z3);
                    zscore_shares[0].push_back(send_z1);
                    zscore_shares[0].push_back(send_z2);
                    zscore_shares[1].push_back(send_z2);
                    zscore_shares[1].push_back(send_z3);
                    zscore_shares[2].push_back(send_z3);
                    zscore_shares[2].push_back(send_z1);
                }
                owner_p1.send(zscore_shares[0]);
                owner_p2.send(zscore_shares[1]);
                owner_p3.send(zscore_shares[2]);
            }
            
            cout << "Sent secret shared row values to parties.\n";
        }
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
void runMPC(string norm_method, int pid,  string ownerIP, int ownerPort, int toOwnerPort,  string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2,int rowstart, int rowend, int numCol, vector<vector<double>>& resultVector) //, atomic<int>& readyCounter, mutex& mtx, condition_variable& cv)
{
    // string to_print = "runMPC: ownerPort: " + to_string(ownerPort);
    // to_print.append(" pid " + to_string(pid));
    // to_print.append(" recPort1 " + to_string(recPort1));
    // to_print.append(" recPort2 " + to_string(recPort2));
    // to_print.append("\n");
    // cout << to_print;
    // mpc testmpc(readyCounter, mtx, cv);
    mpc testmpc;
    // std::promise<void> promise;
    // std::future<void> future = promise.get_future();
    vector<vector<double>> resultVectors;
    testmpc.initialize(pid, ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    if (norm_method == "deseq2")
    {
        testmpc.receivePheno();
    // vector<ZZ_p> random = testmpc.Frand(4);
        testmpc.logRatio();
    }
    auto invcdf_start = chrono::high_resolution_clock::now();
    for (int row=rowstart; row<rowend; row++)
    {
        // cout << string("pid"+to_string(pid)+"_row"+to_string(row)+"\n");
        vector<double> row_result;
        testmpc.ready();
        testmpc.receiveSecrets();
        row_result=testmpc.genperm(row, numCol, norm_method);
        resultVectors.emplace_back(row_result);
        testmpc.clearVectors();
    }
    auto invcdf_end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = invcdf_end - invcdf_start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    cout << string("pid "+to_string(pid)+" inverseCDF Execution time: " + to_string(durationInminutes) + " minutes\n") << endl;
    testmpc.close();
    swap(resultVector, resultVectors);
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
void startMPCset(string norm_method, vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, int numCol, string zscorefile, int CPU_core, vector<vector<double>>& finalVec)
{
    // cout << "startMpcset Start\n";
    vector<thread> threads;
    vector<vector<double>> resultVectors1,resultVectors2,resultVectors3;
    thread dataClientThread(dataclient, norm_method, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], 
    availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow, numCol, zscorefile);
    setThreadAffinity(dataClientThread,CPU_core+0);
    threads.emplace_back(move(dataClientThread));

    // int startidx1 = startidx;
    thread runMPC1([=, &resultVectors1]() {
        runMPC(norm_method, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], 
        startrow, endrow, numCol, resultVectors1);
    });
    // thread runMPC1(runMPC, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], startrow, endrow, resultVectors1);
    setThreadAffinity(runMPC1,CPU_core+1);
    threads.emplace_back(move(runMPC1));

    // // int startidx2 = startidx;
    thread runMPC2([=, &resultVectors2]() {
        runMPC(norm_method, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], 
        startrow, endrow,numCol,resultVectors2);
    });
    // thread runMPC2(runMPC, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], startrow, endrow, resultVectors2);
    setThreadAffinity(runMPC2,CPU_core+2);
    threads.emplace_back(move(runMPC2));

    thread runMPC3([=, &resultVectors3]() {
        runMPC(norm_method, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], 
        startrow, endrow, numCol,resultVectors3);
    });
    // thread runMPC3(runMPC, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], startrow, endrow, resultVectors3);
    setThreadAffinity(runMPC3,CPU_core+2);
    threads.emplace_back(move(runMPC3));
    for (auto& thread : threads) {
        thread.join();
    }
    swap(finalVec, resultVectors1);
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
        cout << "Vector of vectors successfully written to TSV file." << endl;
    }
    else
    {
        cout << "Error opening the file." << endl;
    }
}
double cal_accuracy(vector<vector<double>>& predicted, string& filename, int startrow, int endrow, int numCol)
{
   vector<vector<double>> actual = getMatrixFile(filename,startrow, endrow, numCol);
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
    if (argc < 5)
    {
        cout << "Please provide at least four arg: low, mid, high, samplecount, zscorefile.\n";
        return 1;
    }
    int lo_row = stoi(argv[1]);
    int mid_row = stoi(argv[2]);
    int hi_row = stoi(argv[3]);
    int samplecount = stoi(argv[4]);
    string zscorefile = argv[5];
    string norm_method = argv[6];
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
    int numCores = thread::hardware_concurrency();
    // mutex mtx;
    // condition_variable cv;
    // std::atomic<int> readyCounter(0);
    vector<thread> threads;
    vector<vector<double>> resultVec1,resultVec2,resultVec3,resultVec4;
    thread thread1([&]() {
        startMPCset(norm_method, openPorts, 0, address, lo_row, mid_row, samplecount, zscorefile, 0, resultVec1);
    });
    thread thread2([&]() {
        startMPCset(norm_method, openPorts, 12, address, mid_row, hi_row,samplecount,zscorefile, 0, resultVec2);
    });
    // thread thread3([&]() {
    //     startMPCset(openPorts, 24, address, 2, 3,samplecount,zscorefile, 0, resultVec3);
    // });
    // thread thread4([&]() {
    //     startMPCset(openPorts, 36, address, 3, 4,samplecount,zscorefile, 0, resultVec4);
    // });
    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));
    // threads.emplace_back(move(thread3));
    // threads.emplace_back(move(thread4));
    for (auto& thread : threads) {
        thread.join();
    }
    // cout << resultVec1.size() << endl;
    // cout << resultVec2[0].size() << endl;
    resultVec1.insert(resultVec1.end(), resultVec2.begin(), resultVec2.end());
    string actual_filename;
    if (norm_method == "qn")
    {
        actual_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/qn.tsv";
    }
    else if (norm_method == "deseq2")
    {
        actual_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/deseq2.tsv";
    }
    else 
    {
        throw invalid_argument("either QN or deseq2 normalization please.");
    }
    double accuracy = cal_accuracy(resultVec1, actual_filename, lo_row, hi_row, samplecount);
    cout << string("row "+to_string(lo_row)+" to "+to_string(hi_row)+" accuracy: "+ to_string(accuracy)+"\n");
    writematrixToTSV(resultVec1,lo_row,hi_row,"0728");
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> totalduration = end - start;
    double totaldurationInSeconds = totalduration.count();
    double totaldurationInminutes = totaldurationInSeconds/60.0;
    cout << "Execution time: " << totaldurationInminutes << " minutes" << endl;

    return 0;
}