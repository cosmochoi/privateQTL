#include <thread>
#include "mpc.h"
#include <cmath>
#include "utils.h"
#include <chrono>

// Define the function to be computed
bool computeFunction(bool x1, bool x2, bool x3)
{
    // Example function: output 1 if exactly two inputs are true, 0 otherwise
    return (x1 + x2 + x3) == 2;
}

void bitdecompose(vector<uint64_t> &secrets, vector<BitVector> *bitInput)
{
    // vector<BitVector> bitInput;
    for (int i = 0; i < secrets.size(); i++)
    {
        string sinput = bitset<64>(secrets[i]).to_string();
        BitVector decomposed = BitVector(sinput);
        bitInput->push_back(decomposed);
    }
}

void dataclient(int sendport1, std::string address1, int sendport2, std::string address2, int sendport3, std::string address3)
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

    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/samplewise_QN.tsv";
    vector<BitVector> bitInput; // bitvector representation of input
    vector<double> input = getRowFromMatrixFile(filename,0);
    vector<uint64_t> secrets = ScaleVector(input, pow(10,5)); // integer version of secret
    // print_vector(secrets);
    bitdecompose(secrets, &bitInput);

    NTL::SetSeed((NTL::conv<NTL::ZZ>((long)27))); // Seed change
    uint64_t p = nearestPowerOf2(bitInput.size());
    ZZ_p::init(to_ZZ(p)); 
    uint64_t inv = PowerMod(3, -1, p);
    vector<uint64_t> vectorsize{bitInput.size(), bitInput[0].size(), p, inv};
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
            owner_p1.send(conv<uint64_t>(x1));
            owner_p1.send(conv<uint64_t>(x2));
            owner_p2.send(conv<uint64_t>(x2));
            owner_p2.send(conv<uint64_t>(x3));
            owner_p3.send(conv<uint64_t>(x3));
            owner_p3.send(conv<uint64_t>(x1));
        }
    }


    owner_p1.close();
    owner_p2.close();
    owner_p3.close();
    cout << "Sent secret shared values to parties.\n";
}

// void party(int id, int ownerport, int recPort1, int recPort2, bool input, std::string address1, int sendport1, std::string address2, int sendport2)
// {
//     // Initialize the IOService and create channels to the other parties
//     IOService ios;
//     Endpoint p_owner(ios, "localhost", ownerport, EpMode::Client);
//     Endpoint eprec1(ios, "localhost", recPort1, EpMode::Client);
//     Endpoint eprec2(ios, "localhost", recPort2, EpMode::Client);
//     Endpoint epsend1(ios, address1, sendport1, EpMode::Server);
//     Endpoint epsend2(ios, address2, sendport2, EpMode::Server);
//     Channel dataowner = p_owner.addChannel();
//     Channel chl = eprec1.addChannel();
//     Channel chl1 = eprec2.addChannel();
//     Channel chl2 = epsend1.addChannel();
//     Channel chl3 = epsend2.addChannel();
//     try
//     {
//         // Receive input from data owner
//         // std::vector<u8> share1;
//         // std::vector<u8> share2;
//         BitVector share1 = BitVector();
//         BitVector share2 = BitVector();
//         cout << "ready to receive:\n";
//         dataowner.recv(share1);
//         cout << "received:\n";
//         dataowner.recv(share2);
//         vector<BitVector> recshares{share1, share2};
//         cout << "Party " << id << " received shares:\n";
//         print_vector(recshares);
//         // BitVector shared1 = BitVector(share1.data(), share1.size());
//         // BitVector shared2 = BitVector(share2.data(), share2.size());
//         // print_vector(shares);
//     }
//     catch (const exception &e)
//     {
//         cout << "except\n";
//         cout << e.what() << endl;
//     }

//     // // Share the input with the other parties
//     // uint8_t inputBuffer[1];
//     // inputBuffer[0] = input;
//     // chl2.asyncSend(inputBuffer, 1);
//     // chl3.asyncSend(inputBuffer, 1);

//     // // Receive the inputs from the other parties
//     // uint8_t inputBuffers[2][1];
//     // chl.recv(inputBuffers[0], 1);
//     // chl1.recv(inputBuffers[1], 1);

//     // // Compute the function on the inputs
//     // bool output = computeFunction(input, inputBuffers[0][0], inputBuffers[1][0]);
//     // std::cout << output << "ffff" << std::endl;
//     // // Share the output with the other parties
//     // uint8_t outputBuffer[1];
//     // outputBuffer[0] = output;
//     // chl2.send(outputBuffer, 1);
//     // chl3.send(outputBuffer, 1);

//     // Close the channels
//     chl.close();
//     chl1.close();
//     chl2.close();
//     chl3.close();
// }
vector<ZZ_p> convert(vector<int> input)
{
    vector<ZZ_p> output;
    for (int i=0; i<input.size(); i++)
    {
        output.push_back(conv<ZZ_p>(input[i]));
    }
    return output;
}
void runMPC(int pid, const string ownerIP, int ownerPort, const string address1, int recPort1, int sendPort1, const string address2, int recPort2, int sendPort2)
{
    mpc testmpc;
    std::promise<void> promise;
    std::future<void> future = promise.get_future();
    ZZ_p::init(to_ZZ(16));
    
    // cout << conv<int>(testvec[0]) <<endl;
    testmpc.initialize(pid, ownerIP, ownerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    testmpc.receiveSecrets();
    vector<int> input, rho, sigma;
    vector<ZZ_p> testvec, rhoZZ, sigmaZZ;
    if (pid == 0) {
        input={6,4,8,8,3,2,9,2,7,9,5,3,4,6,3,10,1,13,12,2,7,8,5,2};
        testvec = convert(input);
        rho={3,13,6,4,8,7,5,3};
        rhoZZ=convert(rho);
        sigma={2,8,8,6,9,5,11,3};
        sigmaZZ=convert(sigma);
    }
    else if (pid ==1) {
        input={4,6,8,0,2,12,2,6,9,0,3,9,6,6,10,4,13,3,2,3,8,2,2,9};
        testvec = convert(input);
        rho={13,2,4,9,7,5,3,9};
        rhoZZ=convert(rho);
        sigma={8,7,6,5,5,4,3,6};
        sigmaZZ=convert(sigma);
    } 
    else {
        input={6,6,0,8,12,3,6,9,0,7,9,5,6,4,4,3,3,1,3,12,2,7,9,5};
        testvec = convert(input);
        rho={2,3,9,6,5,8,9,5};
        rhoZZ=convert(rho);
        sigma={7,2,5,8,4,9,6,11};
        sigmaZZ=convert(sigma);
    }
    testmpc.genperm(testvec);   
    // vector<ZZ_p> og{to_ZZ_p(4),to_ZZ_p(2),to_ZZ_p(1),to_ZZ_p(3)};
    // vector<ZZ_p> randpi=testmpc.Frand(4);
    // testmpc.shuffle(randpi, testvec);
    // testmpc.unshuffle(randpi, testvec);
    // vector<ZZ_p> testresult = testmpc.genbitperm(testvec);
    // vector<ZZ_p> reconstructed = testmpc.reveal(testresult);
    // testmpc.apply_shared_perm(testresult, testvec);
    // vector<ZZ_p> sigma
    // testmpc.compose(sigmaZZ,rhoZZ);
    // vector<ZZ_p> test = testmpc.reveal(rhoZZ);
    // testmpc.apply_perm_local(og, rho);
    // print_vector(test);
    promise.set_value(); // signal that initialization is done
    future.get();        // wait for all threads to finish
    testmpc.close();
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

    // // Run each party as a separate thread
    // std::thread t0(dataclient, 4, owner1, address1, owner2, address2, owner3, address3);
    // std::thread t2(party, 0, owner1, port12, port13, false, address2, port21, address3, port31);
    // std::thread t3(party, 1, owner2, port21, port23, false, address3, port12, address1, port32);
    // std::thread t4(party, 2, owner3, port31, port32, true, address1, port13, address2, port23);
    // string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/testdata.txt";
    // vector<BitVector> bitInput;
    // vector<double> input = CSVtoVector(filename);
    // vector<uint64_t> secrets = ScaleVector(input, pow(10,3));
    // bitdecompose(secrets, &bitInput);
    // print_vector(secrets);
    // print_vector(bitInput);
    auto start = std::chrono::high_resolution_clock::now();
    vector<std::future<void>> futures;
    futures.push_back(async(launch::async, dataclient, owner1, address1, owner2, address2, owner3, address3));
    futures.push_back(async(launch::async, runMPC, 0, "localhost", owner1, address2, port12, port21, address3, port13, port31));
    futures.push_back(async(launch::async, runMPC, 1, "localhost", owner2, address3, port23, port32, address1, port21, port12));
    futures.push_back(async(launch::async, runMPC, 2, "localhost", owner3, address1, port31, port13, address2, port32, port23));

    for (auto &future : futures)
    {
        try
        {
            future.get();
        }
        catch (const std::exception &e)
        {
            cerr << "Exception caught: " << e.what() << std::endl;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    double durationInSeconds = duration.count();
    std::cout << "Execution time: " << durationInSeconds << " seconds" << std::endl;

    return 0;
}
