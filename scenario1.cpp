
#include "mpc.h"
#include <cmath>
#include "utils.h"
#include <chrono>
#include "input.h"
#include <sys/resource.h>
struct rusage r_usage;


void dataclient(string norm_method, string split_set, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend, int permut)
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
    uint64_t p = pow(2,50);
    ZZ_p::init(to_ZZ(p)); 
    cout << "P: " << p << endl;
    owner_p1.send(p);
    owner_p2.send(p);
    owner_p3.send(p);
    auto pre_start = chrono::high_resolution_clock::now();

    vector<vector<int64_t>> geno_scaled;
    vector<double> geno_var;
    
    string pheno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/bed_template.tsv";
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_genotype.tsv";
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_genotype_1kGresidualized.tsv";
    string geno_matrix = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/genotype/re/genotype_imputed_projected_residualized_" + split_set +"_concatenated.tsv";
    string geno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_variant.tsv";
    // string pheno_cov = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/pheno_blood_18PCs.tsv";//qn PCA
    // string pheno_cov = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/phenotype/covariate_deseq_pca.csv";
    // vector<vector<double>> covariates = getCovariates(pheno_cov);

    // Residualizer res(covariates);   
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
        vector<string> cisVariants;
        vector<double> std_ratio;
        string geneID;
        string pheno_file;
        if (norm_method == "qn"){
            // pheno_file = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_qn_invCDF.bed";
            pheno_file = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/phenotype/residualized_qn.bed"; 
        }
        else if (norm_method == "tmm"){
            // pheno_file = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/phenotype/sorted_blood_rnaseq_tmm_invCDF.bed";
            pheno_file = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/phenotype/residualized_tmm.bed"; 
        }
        else if (norm_method == "deseq2"){
            // pheno_file = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/phenotype/sorted_blood_rnaseq_deseq_invCDF.bed";
            // pheno_file = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/phenotype/concatenated_residualized_deseq.bed"; //scenario1
            pheno_file = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/private_deseq2_invcdf_"+split_set+"_residualized.bed"; //scenario2_matmult
        }
        // string pheno_file = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/mis_acc.qn.bed";
        vector<double> norm_pheno;
        read_bedfile_row(norm_pheno,geneID, pheno_file, r,4,true);
        // vector<double> pheno_res = res.transform(norm_pheno);
        double bed_var = center_normalize_vec(norm_pheno);
        cout << "pheno var: " << bed_var << endl;
        vector<int64_t> pheno = ScaleVector_signed(norm_pheno, pow(10,6)); // integer version of secret
        
        cout << string("bed file phenotype var: "+to_string(bed_var)+"\n");
        vector<uint64_t> range;
        string chromosome = testinput.getCisRange(geneID,range);

        vector<vector<double>> slicedgeno = testinput.sliceGeno(range, chromosome, -1,cisVariants);
        // vector<vector<double>> geno_res = res.transform(slicedgeno);
        geno_var = center_normalize(slicedgeno);
        
        // vector<vector<double>> geno_res = getMatrixFile("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/pQTL_row_0_geno_centered_resd.tsv", 0, 3000, false, false);
        // geno_var = TSVtoDoubleVector("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/pQTL_row_0_geno_var.tsv");
        geno_scaled = ScaleVector(slicedgeno, pow(10,6));
        // writeVectorToTSV(geno_var, string("pQTL_row_"+to_string(r)+"_geno_var"));
        // writematrixToTSV(geno_res,  string("pQTL_row_"+to_string(r)+"_geno_centered_resd"));
        // writeVectorToTSV(cisVariants, string("pQTL_row_"+to_string(r)+"_cisVar"));
        
        // vector<string> cisVariants = TSVtoVector("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/pQTL_row_0_cisVar.tsv");
        
        for (size_t j = 0; j < geno_var.size(); ++j) {
            std_ratio.push_back(sqrt(bed_var / geno_var[j]));
        }
        // writeVectorToTSV(std_ratio, "pQTL_std_ratio");
        // std_ratio = TSVtoDoubleVector("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/pQTL_std_ratio.tsv");
        uint64_t inv = PowerMod(3, -1, p);
        vector<vector<uint64_t>> shares;
        for (int i = 0; i < 3; i++) {
            shares.push_back(vector<uint64_t>());
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

            uint64_t send_x1 = conv<uint64_t>(x1);
            uint64_t send_x2 = conv<uint64_t>(x2);
            uint64_t send_x3 = conv<uint64_t>(x3);
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
        
        cout << "Sent secret shared gene pheno values to parties.\n";
        vector<vector<uint64_t>> genoshares;
        for (int i = 0; i < 3; i++) {
            genoshares.push_back(vector<uint64_t>());
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

                uint64_t send_i1 = conv<uint64_t>(i1);
                uint64_t send_i2 = conv<uint64_t>(i2);
                uint64_t send_i3 = conv<uint64_t>(i3);
                genoshares[0].push_back(send_i1);
                genoshares[0].push_back(send_i2);
                genoshares[1].push_back(send_i2);
                genoshares[1].push_back(send_i3);
                genoshares[2].push_back(send_i3);
                genoshares[2].push_back(send_i1);
            }
        }
        uint64_t testg_row = geno_scaled.size();
        uint64_t testg_col = geno_scaled[0].size();
        vector<uint64_t> genoshape = {
            static_cast<uint64_t>(testg_row),
            static_cast<uint64_t>(testg_col),
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
        int complete;
        p1_owner.recv(complete);
        if (complete != 1)
        {
            cout << "Didn't receive confirmation. Did row finish?" << endl;
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
void runMPC(string norm_method, int pid,  string ownerIP, int ownerPort, int toOwnerPort,  string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2,int rowstart, int rowend, int permut, Logger& cisLogger, Logger& nominalLogger) //, atomic<int>& readyCounter, mutex& mtx, condition_variable& cv)
{    
    mpc testmpc;
    testmpc.initialize(pid, ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    auto invcdf_start = chrono::high_resolution_clock::now();
    for (int row=rowstart; row<rowend; row++)
    {
        testmpc.permutPheno(permut);
        testmpc.receiveGeno();
        testmpc.calc_corr(cisLogger, nominalLogger);
        testmpc.clearVectors();
    }
    auto invcdf_end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = invcdf_end - invcdf_start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    testmpc.close();
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
void startMPCset(string norm_method, string split_set, vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, int CPU_core, int permut, Logger& cislogger, Logger& nominalLogger)
{
    vector<thread> threads;
    thread dataClientThread(dataclient, norm_method, split_set, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], 
    availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow, permut);
    setThreadAffinity(dataClientThread,CPU_core+0);
    threads.emplace_back(move(dataClientThread));

    // int startidx1 = startidx;
    thread runMPC1([=, &cislogger, &nominalLogger]() {
        runMPC(norm_method, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], 
        startrow, endrow, permut,cislogger,nominalLogger);
    });
    setThreadAffinity(runMPC1,CPU_core+1);
    threads.emplace_back(move(runMPC1));

    thread runMPC2([=, &cislogger,&nominalLogger]() {
        runMPC(norm_method, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], 
        startrow, endrow,permut,cislogger,nominalLogger);
    });
    setThreadAffinity(runMPC2,CPU_core+2);
    threads.emplace_back(move(runMPC2));

    thread runMPC3([=, &cislogger,&nominalLogger]() {
        runMPC(norm_method, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], 
        startrow, endrow, permut,cislogger,nominalLogger);
    });
    setThreadAffinity(runMPC3,CPU_core+2);
    threads.emplace_back(move(runMPC3));
    for (auto& thread : threads) {
        thread.join();
    }
}

int main(int argc, char* argv[])
{
    if (argc < 8)
    {
        cout << "Please provide at least four arg: low, mid, high, permutation, norm method, cislog, nominallog.\n";
        return 1;
    }
    int lo_row = stoi(argv[1]);
    int mid_row = stoi(argv[2]);
    int hi_row = stoi(argv[3]);
    int permut = stoi(argv[4]);
    // string zscorefile = argv[5]; // don't need zscores for scenario1
    string norm_method = argv[5];
    string split_set = argv[6];
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
    vector<string> address{"localhost","localhost","localhost","localhost"};

    auto start = chrono::high_resolution_clock::now();
    int numCores = thread::hardware_concurrency();
    vector<thread> threads;
    string nominal = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/output/" + split_set + "/privateQTL_scenario1_"+split_set+"_"+nominal_log);
    string cis = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/output/" + split_set + "/privateQTL_scenario1_"+split_set+"_"+cis_log);
    Logger nominallogger(nominal+"_"+to_string(lo_row)+"_"+to_string(mid_row)+".tsv"), cislogger(cis+"_"+to_string(lo_row)+"_"+to_string(mid_row)+".tsv");
    Logger nominallogger2(nominal+"_"+to_string(mid_row)+"_"+to_string(hi_row)+".tsv"),cislogger2(cis+"_"+to_string(mid_row)+"_"+to_string(hi_row)+".tsv");
    nominallogger.log(string("phenID\tvarID\tvarIdx\tdof\tr_nom\tr2_nom\ttstat\tpval\tslope\tslope_se"));
    cislogger.log(string("phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\tpval_true_df\tr_nom\tr2_nom\ttstat\tpval_nominal\tslope\tslope_se\tpval_perm\tpval_beta"));
    nominallogger2.log(string("phenID\tvarID\tvarIdx\tdof\tr_nom\tr2_nom\ttstat\tpval\tslope\tslope_se"));
    cislogger2.log(string("phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\tpval_true_df\tr_nom\tr2_nom\ttstat\tpval_nominal\tslope\tslope_se\tpval_perm\tpval_beta"));
    thread thread1([&]() {
        startMPCset(norm_method, split_set, openPorts, 0, address, lo_row, mid_row, 0,permut,cislogger,nominallogger);
    });
    thread thread2([&]() {
        startMPCset(norm_method, split_set, openPorts, 12, address, mid_row, hi_row, 0, permut,cislogger2,nominallogger2);
    });

    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));

    for (auto& thread : threads) {
        thread.join();
    }

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