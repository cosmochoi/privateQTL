
#include "mpc.h"
#include <cmath>
#include "utils.h"
#include <chrono>


void bitdecompose(vector<uint32_t> &secrets, vector<BitVector> *bitInput)
{
    // vector<BitVector> bitInput;
    for (int i = 0; i < secrets.size(); i++)
    {
        string sinput = bitset<32>(secrets[i]).to_string();
        BitVector decomposed(sinput);
        bitInput->push_back(decomposed);
    }
}

void dataclient(int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend,int numCol, string zscorefile)
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
    // PRNG ownerprng;
    // ownerprng.SetSeed(toBlock(27));
    for (int r=rowstart; r<rowend; r++)
    {
        int ready1, ready2, ready3;
        p1_owner.recv(ready1);
        p2_owner.recv(ready2);
        p3_owner.recv(ready3);
        if ((ready1 == 1) && (ready2 ==1) && (ready3==1))
        {
            vector<BitVector> bitInput;
            // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/samplewise_QN.tsv";
            string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/correct_sampleQN.tsv";
            // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/sample_thres125_fp1.tsv";
            vector<double> input = getRowFromMatrixFile(filename,r, numCol);

            // string filename ="/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/testdata2.txt";
            // vector<double> input = CSVtoVector(filename);
            // bitvector representation of input
            
            vector<uint32_t> secrets = ScaleVector(input, pow(10,2)); // integer version of secret
            bitdecompose(secrets, &bitInput);

            NTL::SetSeed((NTL::conv<NTL::ZZ>((long)27))); // Seed change
            cout << "size of input: " << bitInput.size() << endl;
            cout << "size of each vector: " << bitInput[0].size() << endl;

            uint32_t p = nearestPowerOf2(bitInput.size());
            // uint32_t p=32768;
            cout << "P: " << p << endl;
            ZZ_p::init(to_ZZ(p)); 
            uint32_t inv = PowerMod(3, -1, p);
            vector<uint32_t> vectorsize{(uint32_t) bitInput.size(), (uint32_t) bitInput[0].size(), p, inv};
            owner_p1.send(vectorsize);
            owner_p2.send(vectorsize);
            owner_p3.send(vectorsize);

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
            string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/"+zscorefile+".txt");
            // string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/" + zscorefile + ".txt");
            vector<double> zscore_input = CSVtoVector(zscore_filename);
            // print_vector(zscore_input);
            // vector<int32_t> zscores_int = ScaleVector_signed(zscore_input, pow(10,5));
            auto min = min_element(zscore_input.begin(), zscore_input.end());
            double shiftsize = abs(*min);
            // cout << "shift size: " << shiftsize << endl;
            vector<double> zscores_shifted = ShiftVector(zscore_input, shiftsize);
            vector<uint32_t> zscores_scaled = ScaleVector(zscores_shifted, pow(10,2));
            cout << "zscores size: " << zscores_scaled.size() << endl;
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
            cout << "Sent secret shared values to parties.\n";
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
void runMPC(int pid,  string ownerIP, int ownerPort, int toOwnerPort,  string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2,int rowstart, int rowend, int numCol, vector<vector<double>>& resultVector) //, atomic<int>& readyCounter, mutex& mtx, condition_variable& cv)
{
    string to_print = "runMPC: ownerPort: " + to_string(ownerPort);
    to_print.append(" pid " + to_string(pid));
    to_print.append(" recPort1 " + to_string(recPort1));
    to_print.append(" recPort2 " + to_string(recPort2));
    to_print.append("\n");
    cout << to_print;
    // mpc testmpc(readyCounter, mtx, cv);
    mpc testmpc;
    // std::promise<void> promise;
    // std::future<void> future = promise.get_future();
    vector<vector<double>> resultVectors;
    testmpc.initialize(pid, ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    
    for (int row=rowstart; row<rowend; row++)
    {
        cout << string("pid"+to_string(pid)+"_row"+to_string(row)+"\n");
        vector<double> row_result;
        testmpc.ready();
        testmpc.receiveSecrets();
        row_result=testmpc.genperm(row, numCol);
        resultVectors.emplace_back(row_result);
        testmpc.clearVectors();
    }
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
void startMPCset(vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, int numCol, string zscorefile, int CPU_core, vector<vector<double>>& finalVec)
{
    cout << "startMpcset Start\n";
    vector<thread> threads;
    vector<vector<double>> resultVectors1,resultVectors2,resultVectors3;
    thread dataClientThread(dataclient, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], 
    availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow, numCol, zscorefile);
    setThreadAffinity(dataClientThread,CPU_core+0);
    threads.emplace_back(move(dataClientThread));

    string to_print = "1 setup Channels: ownerPort: " + to_string(availPorts[startidx + 6]);
    to_print.append(" recPort1 " + to_string(availPorts[startidx + 0]));
    to_print.append(" recPort2 " + to_string(availPorts[startidx + 1]));
    to_print.append("\n");
    cout << to_print;

    // int startidx1 = startidx;
    thread runMPC1([=, &resultVectors1]() {
        runMPC(0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], 
        startrow, endrow, numCol, resultVectors1);
    });
    // thread runMPC1(runMPC, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], startrow, endrow, resultVectors1);
    setThreadAffinity(runMPC1,CPU_core+1);
    threads.emplace_back(move(runMPC1));

    to_print = "2 setup Channels: ownerPort: " + to_string(availPorts[startidx + 7]);
    to_print.append(" recPort1 " + to_string(availPorts[startidx + 3]));
    to_print.append(" recPort2 " + to_string(availPorts[startidx + 2]));
    to_print.append("\n");
    cout << to_print;

    // // int startidx2 = startidx;
    thread runMPC2([=, &resultVectors2]() {
        runMPC(1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], 
        startrow, endrow,numCol,resultVectors2);
    });
    // thread runMPC2(runMPC, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], startrow, endrow, resultVectors2);
    setThreadAffinity(runMPC2,CPU_core+2);
    threads.emplace_back(move(runMPC2));

    to_print = "3 setup Channels: ownerPort: " + to_string(availPorts[startidx + 8]);
    to_print.append(" recPort1 " + to_string(availPorts[startidx + 4]));
    to_print.append(" recPort2 " + to_string(availPorts[startidx + 5]));
    to_print.append("\n");
    cout << to_print;

    thread runMPC3([=, &resultVectors3]() {
        runMPC(2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], 
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
    print_vector(openPorts);
    vector<string> address{"localhost","localhost","localhost","localhost"};

    auto start = chrono::high_resolution_clock::now();
    int numCores = thread::hardware_concurrency();
    // mutex mtx;
    // condition_variable cv;
    // std::atomic<int> readyCounter(0);
    vector<thread> threads;
    vector<vector<double>> resultVec1,resultVec2;
    thread thread1([&]() {
        startMPCset(openPorts, 0, address, lo_row, mid_row, samplecount, zscorefile, 0, resultVec1);
    });
    thread thread2([&]() {
        startMPCset(openPorts, 12, address, mid_row, hi_row,samplecount,zscorefile, 0, resultVec2);
    });
    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));
    for (auto& thread : threads) {
        thread.join();
    }
    cout << resultVec1.size() << endl;
    // cout << resultVec2[0].size() << endl;
    resultVec1.insert(resultVec1.end(), resultVec2.begin(), resultVec2.end());
    cout << resultVec1.size() << endl;
    // startMPCset(openPorts, 0, address, lo_row, mid_row, samplecount, zscorefile, 0, resultVec1);
    // startMPCset(openPorts, 12, address, mid_row, hi_row,samplecount,zscorefile, 0, resultVec2);
    // startMPCset(openPorts, 24, address, 2, 3,threads, 4);
    // startMPCset(openPorts, 36, address, 3, 4,threads, 4);
        // threads.emplace_back(dataclient, owner1, toOwner1, address1, owner2, toOwner2, address2, owner3, toOwner3,address3);
        // threads.emplace_back(runMPC, 0, "localhost", owner1, toOwner1, address2, port12, port21, address3, port13, port31);
        // threads.emplace_back(runMPC, 1, "localhost", owner2, toOwner2, address3, port23, port32, address1, port21, port12);
        // threads.emplace_back(runMPC, 2, "localhost", owner3, toOwner3, address1, port31, port13, address2, port32, port23);

    // thread dataClientThread(dataclient, owner1, toOwner1, address1, owner2, toOwner2, address2, owner3, toOwner3, address3, lo_row, mid_row);
    // setThreadAffinity(dataClientThread,0);
    // threads.emplace_back(move(dataClientThread));
    // thread runMPC1(runMPC, 0, data_address, owner1, toOwner1, address2, port12, port21, address3, port13, port31, lo_row, mid_row);
    // setThreadAffinity(runMPC1,1);
    // threads.emplace_back(move(runMPC1));
    // thread runMPC2(runMPC, 1, data_address, owner2, toOwner2, address3, port23, port32, address1, port21, port12, lo_row, mid_row);
    // setThreadAffinity(runMPC2,2);
    // threads.emplace_back(move(runMPC2));
    // thread runMPC3(runMPC, 2, data_address, owner3, toOwner3, address1, port31, port13, address2, port32, port23, lo_row, mid_row);
    // setThreadAffinity(runMPC3,2);
    // threads.emplace_back(move(runMPC3));

    // thread dataClientThread2(dataclient, owner1_2, toOwner1_2, address1, owner2_2, toOwner2_2, address2, owner3_2, toOwner3_2, address3, 2, hi_row);
    // setThreadAffinity(dataClientThread2,0);
    // threads.emplace_back(move(dataClientThread2));
    // thread runMPC1_2(runMPC, 0, data_address, owner1_2, toOwner1_2, address2, port12_2, port21_2, address3, port13_2, port31_2,2, hi_row);
    // setThreadAffinity(runMPC1_2,1);
    // threads.emplace_back(move(runMPC1_2));
    // thread runMPC2_2(runMPC, 1, data_address, owner2_2, toOwner2_2, address3, port23_2, port32_2, address1, port21_2, port12_2,2, hi_row);
    // setThreadAffinity(runMPC2_2,2);
    // threads.emplace_back(move(runMPC2_2));
    // thread runMPC3_2(runMPC, 2, data_address, owner3_2, toOwner3_2, address1, port31_2, port13_2, address2, port32_2, port23_2,2, hi_row);
    // setThreadAffinity(runMPC3_2,2);
    // threads.emplace_back(move(runMPC3_2));    

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    cout << "Execution time: " << durationInminutes << " minutes" << endl;

    return 0;
}