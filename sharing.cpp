
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

void dataclient(int sendport1, std::string address1, int sendport2, std::string address2, int sendport3, std::string address3, int r)
{
    IOService ios;
    Endpoint epsend1(ios, address1, sendport1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendport2, EpMode::Server);
    Endpoint epsend3(ios, address3, sendport3, EpMode::Server);
    Channel owner_p1 = epsend1.addChannel();
    Channel owner_p2 = epsend2.addChannel();
    Channel owner_p3 = epsend3.addChannel();
    cout << "Owner established channels with the computing parties.\n";
    // PRNG ownerprng;
    // ownerprng.SetSeed(toBlock(27));
    vector<BitVector> bitInput;
    // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/samplewise_QN.tsv";
    // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/correct_sampleQN.tsv";
    // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/sample_thres125_fp0.tsv";
    // vector<double> input = getRowFromMatrixFile(filename,r);

    string filename ="/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/testdata2.txt";
    vector<double> input = CSVtoVector(filename);
    // bitvector representation of input
    
    
    vector<uint32_t> secrets = ScaleVector(input, pow(10,5)); // integer version of secret
    // print_vector(secrets);
    bitdecompose(secrets, &bitInput);

    NTL::SetSeed((NTL::conv<NTL::ZZ>((long)27))); // Seed change
    uint32_t p = nearestPowerOf2(bitInput.size());
    // uint32_t p = 32768;
    cout << "P: " << p << endl;
    ZZ_p::init(to_ZZ(p)); 
    // ZZ_p::init(to_ZZ(32768)); 
    uint32_t inv = PowerMod(3, -1, p);
    vector<uint32_t> vectorsize{(uint32_t) bitInput.size(), (uint32_t) bitInput[0].size(), p, inv};
    // cout << "sent size: " << bitInput.size() << endl;
    owner_p1.send(vectorsize);
    owner_p2.send(vectorsize);
    owner_p3.send(vectorsize);
    // print_vector(bitInput);
    // making the shares: for each bit, loop through secrets
    for (int j=bitInput[0].size()-1; j >=0; j--)
    {
        for (int i=0; i<bitInput.size(); i++)
        {
            ZZ_p x1 = random_ZZ_p();
            ZZ_p x2 = random_ZZ_p();
            ZZ_p x3 = to_ZZ_p(bitInput[i][j]) - x1 - x2;
            owner_p1.send(conv<uint32_t>(x1));
            owner_p1.send(conv<uint32_t>(x2));
            owner_p2.send(conv<uint32_t>(x2));
            owner_p2.send(conv<uint32_t>(x3));
            owner_p3.send(conv<uint32_t>(x3));
            owner_p3.send(conv<uint32_t>(x1));
        }
    }
    // string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/zscores.txt";
    string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/zscores2.txt";
    vector<double> zscore_input = CSVtoVector(zscore_filename);
    // print_vector(zscore_input);
    // vector<int32_t> zscores_int = ScaleVector_signed(zscore_input, pow(10,5));
    auto min = min_element(zscore_input.begin(), zscore_input.end());
    double shiftsize = abs(*min);
    cout << "shift size: " << shiftsize << endl;
    vector<double> zscores_shifted = ShiftVector(zscore_input, shiftsize);
    // print_vector(zscores_shifted);
    vector<uint32_t> zscores_scaled = ScaleVector(zscores_shifted, pow(10,3));
    // print_vector(zscores_scaled);
    //sending identity shares
    for (int j=0; j< bitInput.size(); j++)
    {
        ZZ_p i1 = random_ZZ_p();
        ZZ_p i2 = random_ZZ_p();
        ZZ_p i3 = conv<ZZ_p>(j+1) - i1 - i2;
        owner_p1.send(conv<uint32_t>(i1));
        owner_p1.send(conv<uint32_t>(i2));
        owner_p2.send(conv<uint32_t>(i2));
        owner_p2.send(conv<uint32_t>(i3));
        owner_p3.send(conv<uint32_t>(i3));
        owner_p3.send(conv<uint32_t>(i1));
    }
    owner_p1.send(shiftsize);
    owner_p2.send(shiftsize);
    owner_p3.send(shiftsize);
    for (int i=0; i<zscores_scaled.size(); i++)
    {
        ZZ_p z1 = random_ZZ_p();
        ZZ_p z2 = random_ZZ_p();
        ZZ_p z3 = to_ZZ_p(zscores_scaled[i]) - z1 - z2;
        owner_p1.send(conv<uint32_t>(z1));
        owner_p1.send(conv<uint32_t>(z2));
        owner_p2.send(conv<uint32_t>(z2));
        owner_p2.send(conv<uint32_t>(z3));
        owner_p3.send(conv<uint32_t>(z3));
        owner_p3.send(conv<uint32_t>(z1));
    }


    owner_p1.close();
    owner_p2.close();
    owner_p3.close();
    epsend1.stop();
    epsend2.stop();
    epsend3.stop();
    ios.stop();
    cout << "Sent secret shared values to parties.\n";
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
void runMPC(int pid, const string ownerIP, int ownerPort, const string address1, int recPort1, int sendPort1, 
const string address2, int recPort2, int sendPort2, int row)//, atomic<int>& readyCounter, mutex& mtx, condition_variable& cv)
{
    // mpc testmpc(readyCounter, mtx, cv);
    mpc testmpc;
    // std::promise<void> promise;
    // std::future<void> future = promise.get_future();
    // ZZ_p::init(to_ZZ(16));
    
    // cout << conv<int>(testvec[0]) <<endl;
    testmpc.initialize(pid, ownerIP, ownerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    testmpc.receiveSecrets();
    // testmpc.receiveOriginalData();
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


    testmpc.genperm(row);   

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

int main()
{
    int port12 = 12345;            // p2 to p1
    int port13 = 12346;            // p3 to p1
    int port21 = 12355;            // p1 to p2
    int port23 = 12356;            // p3 to p2
    int port31 = 12365;            // p1 to p3
    int port32 = 12366;            // p2 to p3
    int owner1 = 12375;            // owner to p1
    int owner2 = 12376;            // owner to p2
    int owner3 = 12377;            // owner to p3
    string address1 = "localhost"; // address of the remote ssh server for party 1
    string address2 = "localhost"; // address of the remote ssh server for party 2
    string address3 = "localhost"; // address of the remote ssh server for party 3

    auto start = std::chrono::high_resolution_clock::now();
    int numCores = std::thread::hardware_concurrency();
    // mutex mtx;
    // condition_variable cv;
    // std::atomic<int> readyCounter(0);
    std::vector<std::thread> threads;
    for (int row = 0; row < 1; row++)
    {
        threads.emplace_back(dataclient, owner1, address1, owner2, address2, owner3, address3, row);
        // threads.emplace_back(runMPC, 0, "localhost", owner1, address2, port12, port21, address3, port13, port31,ref(readyCounter),ref(mtx), ref(cv));
        // threads.emplace_back(runMPC, 1, "localhost", owner2, address3, port23, port32, address1, port21, port12,ref(readyCounter),ref(mtx), ref(cv));
        // threads.emplace_back(runMPC, 2, "localhost", owner3, address1, port31, port13, address2, port32, port23,ref(readyCounter),ref(mtx), ref(cv));
        threads.emplace_back(runMPC, 0, "localhost", owner1, address2, port12, port21, address3, port13, port31, row);
        threads.emplace_back(runMPC, 1, "localhost", owner2, address3, port23, port32, address1, port21, port12, row);
        threads.emplace_back(runMPC, 2, "localhost", owner3, address1, port31, port13, address2, port32, port23, row);
        // Set affinity for each thread to a separate CPU core
        // for (int i = 0; i < 4; i++)
        // {
        //     int coreId = i % numCores;
        //     setThreadAffinity(threads[i], coreId);
        // }

        // Wait for all threads to finish
        for (auto& thread : threads)
        {
            thread.join();
        }

        threads.clear();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    double durationInSeconds = duration.count();
    std::cout << "Execution time: " << durationInSeconds << " seconds" << std::endl;

    return 0;
}