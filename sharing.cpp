
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

void dataclient(int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend)
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
            string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/samplewise_QN.tsv";
            // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/correct_sampleQN.tsv";
            // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/sample_thres125_fp1.tsv";
            vector<double> input = getRowFromMatrixFile(filename,r);

            // string filename ="/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/testdata2.txt";
            // vector<double> input = CSVtoVector(filename);
            // bitvector representation of input
            
            vector<uint32_t> secrets = ScaleVector(input, pow(10,5)); // integer version of secret
            bitdecompose(secrets, &bitInput);

            NTL::SetSeed((NTL::conv<NTL::ZZ>((long)27))); // Seed change
            cout << "size of input: " << bitInput.size() << endl;
            cout << "size of each vector: " << bitInput[0].size() << endl;

            // uint32_t p = nearestPowerOf2(bitInput.size());
            uint32_t p=32768;
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
            // string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/zscores.txt";
            string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/zscores.txt";
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
void runMPC(int pid, const string ownerIP, int ownerPort, int toOwnerPort, const string address1, int recPort1, int sendPort1, 
const string address2, int recPort2, int sendPort2,int rowstart, int rowend)//, atomic<int>& readyCounter, mutex& mtx, condition_variable& cv)
{
    // mpc testmpc(readyCounter, mtx, cv);
    mpc testmpc;
    // std::promise<void> promise;
    // std::future<void> future = promise.get_future();
    // ZZ_p::init(to_ZZ(16));
    
    // cout << conv<int>(testvec[0]) <<endl;
    testmpc.initialize(pid, ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    
    for (int row=rowstart; row<rowend; row++)
    {
        testmpc.ready();
        testmpc.receiveSecrets();
        testmpc.genperm(row);
        testmpc.clearVectors();
    }
    
    // vector<int> input, rho, sigma;
    // vector<ZZ_p> testvec, rhoZZ, sigmaZZ;
    // if (pid == 0) {
    //     // input={6,4,8,8,3,2,9,2,7,9,5,3,4,6,3,10,1,13,12,2,7,8,5,2};
    //     // testvec = convert(input);
    //     rho={3,13,6,4,8,7,5,3};
    //     rhoZZ=convert(rho);
    //     // sigma={2,8,8,6,9,5,11,3};
    //     // sigmaZZ=convert(sigma);
    // }
    // else if (pid ==1) {
    //     // input={4,6,8,0,2,12,2,6,9,0,3,9,6,6,10,4,13,3,2,3,8,2,2,9};
    //     // testvec = convert(input);
    //     rho={13,2,4,9,7,5,3,9};
    //     rhoZZ=convert(rho);
    //     // sigma={8,7,6,5,5,4,3,6};
    //     // sigmaZZ=convert(sigma);
    // } 
    // else {
    //     // input={6,6,0,8,12,3,6,9,0,7,9,5,6,4,4,3,3,1,3,12,2,7,9,5};
    //     // testvec = convert(input);
    //     rho={2,3,9,6,5,8,9,5};
    //     rhoZZ=convert(rho);
    //     // sigma={7,2,5,8,4,9,6,11};
    //     // sigmaZZ=convert(sigma);
    // }
    // vector<ZZ_p> rho2(rhoZZ);
    // testmpc.compose(rhoZZ, rho2);
    // print_vector(testmpc.reveal(rho2, false));


       

    // vector<ZZ_p> og{to_ZZ_p(4),to_ZZ_p(2),to_ZZ_p(1),to_ZZ_p(3)};
    //         for (int i = 0; i < 100; i++) {

    // vector<ZZ_p> randpi=testmpc.Frand(4);
    // vector<ZZ_p> randrho=testmpc.Frand(4);

    // vector<ZZ_p> reconstructedPi = testmpc.reveal(randpi, true);
    // cout << "pi: " << endl;
    // print_vector(randpi);
    // vector<ZZ_p> reconstructedRho = testmpc.reveal(rhoZZ, false);
    // cout <<"original: "<<endl;
    // print_vector(reconstructedRho);
    // vector<ZZ_p> reconstructedRho = testmpc.reveal(testvec, false);
    // cout << "original: "<<endl;
    // print_vector(reconstructedRho);
    // for(int i = 0; i < 0; i++) {
    //     cout << "REVEAL" + i << endl;
    //     assert_equal(reconstructed, testmpc.reveal(randpi), "reveal" + i);
    // }

    // for(int j = 0; j < 1000; j++) {
    //     testmpc.reshare(randpi, 1);
    //     assert_equal(reconstructed, testmpc.reveal(randpi), "reshare" + j);
    // }

    // testmpc.shuffle(randpi, rhoZZ);
    // testmpc.shuffle(randpi, testvec);
    // cout << "after shuffle: " << endl;
    // print_vector(testmpc.reveal(testvec, false));

    // vector<ZZ_p> afterRho = testmpc.reveal(randrho, false);
    // cout <<"After shuffle: " << endl;
    // print_vector(afterRho);
    // testmpc.unshuffle(randpi, testvec);
    // testmpc.unshuffle(randpi, rhoZZ);

    // cout <<"After unshuffle: " << endl;
    // vector<ZZ_p> unshuffledRho = testmpc.reveal(rhoZZ, false);
    // vector<ZZ_p> unshuffledRho = testmpc.reveal(testvec, false);
    // print_vector(unshuffledRho);
    // assert_equal(reconstructedRho, unshuffledRho, "unshuffle");
    //         }
    // } catch (const std::exception &e)
    //         {
    //             cerr << "Exception caught: " << e.what() << std::endl;
    //         }

    // print_vector(testmpc.reveal(randrho));
    // vector<ZZ_p> testresult = testmpc.genbitperm(testvec);
    // vector<ZZ_p> reconstructed = testmpc.reveal(testresult);
    // testmpc.apply_shared_perm(testresult, testvec);
    // vector<ZZ_p> sigma
    // testmpc.compose(sigmaZZ,rhoZZ);
    // vector<ZZ_p> test = testmpc.reveal(rhoZZ);
    // testmpc.apply_perm_local(og, rho);
    // print_vector(test);
    // promise.set_value(); // signal that initialization is done
    // future.get();        // wait for all threads to finish
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
void startMPCset(vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, vector<thread>& threads)
{
    thread dataClientThread(dataclient, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow);
    setThreadAffinity(dataClientThread,0);
    threads.emplace_back(move(dataClientThread));
    thread runMPC1(runMPC, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], startrow, endrow);
    setThreadAffinity(runMPC1,1);
    threads.emplace_back(move(runMPC1));
    thread runMPC2(runMPC, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], startrow, endrow);
    setThreadAffinity(runMPC2,2);
    threads.emplace_back(move(runMPC2));
    thread runMPC3(runMPC, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], startrow, endrow);
    setThreadAffinity(runMPC3,2);
    threads.emplace_back(move(runMPC3));
}
int main()
{
    int lo_row = 0;
    int mid_row = 1;
    int hi_row = 3;
    int startingPort = 12300;
    int assignedPorts=0;
    vector<int> openPorts;
    while (assignedPorts < 48)
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
    // int port12 = openPorts[0];//12345;            // p2 to p1
    // int port13 = openPorts[1]; //12346;            // p3 to p1
    // int port21 = openPorts[2];//12355;            // p1 to p2
    // int port23 = openPorts[3];//12356;            // p3 to p2
    // int port31 = openPorts[4];//12365;            // p1 to p3
    // int port32 = openPorts[5];//12366;            // p2 to p3
    // int owner1 = openPorts[6];//12375;            // owner to p1
    // int owner2 = openPorts[7];//12376;            // owner to p2
    // int owner3 = openPorts[8];//12377;            // owner to p3
    // int toOwner1 = openPorts[9];//12378;           // p1 to owner
    // int toOwner2 = openPorts[10];//12379;           //p2 to owner
    // int toOwner3 = openPorts[11];//12380;           //p3 to owner
    vector<string> address{"localhost","localhost","localhost","localhost"};

    auto start = std::chrono::high_resolution_clock::now();
    int numCores = std::thread::hardware_concurrency();
    // mutex mtx;
    // condition_variable cv;
    // std::atomic<int> readyCounter(0);
    std::vector<std::thread> threads;
    startMPCset(openPorts, 0, address, 0, 1,threads);
    startMPCset(openPorts, 12, address, 1, 2,threads);
    startMPCset(openPorts, 24, address, 2, 3,threads);
    startMPCset(openPorts, 36, address, 3, 4,threads);
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


    // int port12_2 = 12381;            // p2 to p1
    // int port13_2 = 12382;            // p3 to p1
    // int port21_2 = 12383;            // p1 to p2
    // int port23_2 = 12384;            // p3 to p2
    // int port31_2 = 12385;            // p1 to p3
    // int port32_2 = 12386;            // p2 to p3
    // int owner1_2 = 12387;            // owner to p1
    // int owner2_2 = 12388;            // owner to p2
    // int owner3_2 = 12389;            // owner to p3
    // int toOwner1_2 = 12390;           // p1 to owner
    // int toOwner2_2 = 12391;           //p2 to owner
    // int toOwner3_2 = 12392;           //p3 to owner
    // // // string address1 = "localhost"; // address of the remote ssh server for party 1
    // // // string address2 = "localhost"; // address of the remote ssh server for party 2
    // // // string address3 = "localhost"; // address of the remote ssh server for party 3
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
    for (auto& thread : threads)
    {
        thread.join();
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    cout << "Execution time: " << durationInminutes << " minutes" << endl;

    return 0;
}