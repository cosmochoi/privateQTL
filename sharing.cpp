#include <thread>
#include "mpc.h"
#include <cmath>
#include "utils.h"

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

// Flatten a vector<vector<BitVector>> object into a one-dimensional byte vector
vector<uint8_t> flatten(const vector<vector<BitVector>>& data) {
  vector<uint8_t> result;
  for (const auto& row : data) {
    for (const auto& elem : row) {
      const uint8_t* bytes = elem.data();
      size_t num_bytes = elem.sizeBytes();
      result.insert(result.end(), bytes, bytes + num_bytes);
    }
  }
  return result;
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
    PRNG ownerprng;
    ownerprng.SetSeed(toBlock(27));

    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/testdata.txt";
    vector<BitVector> bitInput; // bitvector representation of input
    vector<double> input = CSVtoVector(filename);
    vector<uint64_t> secrets = ScaleVector(input, pow(10,3)); // integer version of secret
    bitdecompose(secrets, &bitInput);

    // print_vector(secrets);
    // print_vector(bitInput);

    vector<size_t> vectorsize{bitInput.size()};
    owner_p1.send(vectorsize);
    owner_p2.send(vectorsize);
    owner_p3.send(vectorsize);
    uint32_t n = nearestPowerOf2(bitInput.size());
    cout << "N: " << n << endl;
    // int p = 11;
    uint32_t p = PowerMod(3, -1, n);
    cout << "P: " << p << endl;
    // // making the shares
    // for (int i=0; i < bitInput.size(); i++){
    //     vector<BitVector> triples(3);
    //     triples[2].assign(bitInput[i]);
    //     for (int j = 0; j < 2; j++)
    //     {
    //         triples[j].assign(bitInput[i]);
    //         triples[j].randomize(ownerprng);
    //         triples[2] ^= triples[j];
    //     }

    //     owner_p1.send(triples[0]);
    //     owner_p1.send(triples[1]);
    //     owner_p2.send(triples[1]);
    //     owner_p2.send(triples[2]);
    //     owner_p3.send(triples[2]);
    //     owner_p3.send(triples[0]);
    // }
    for (int i=0; i< 3; i++)
    {
        BitVector bitSecret = bitInput[i];
        int randomNumbers[3];
        int sum =0;
        cout << bitSecret << endl;
        for (int j=3; j>=0; --j) //starting from the LSB
        {
            for (int m=0; i<3; m++)
            {
                randomNumbers[m] = ownerprng.get<uint32_t>();
                sum += randomNumbers[m];
            }
            int result = sum % 11;
            for (int m = 0; m < 3; m++) 
            {
            randomNumbers[m] = (randomNumbers[m] + (result - bitSecret[j]) + 11) % 11;
            }
        }
        
    }

    // auto share1_bytes = flatten(share1);
    // auto share2_bytes = flatten(share2);
    // auto share3_bytes = flatten(share3);
    // cout << "Data owner sending shares to parties:\n";
    // owner_p1.send(span<uint8_t>(share1_bytes.data(), share1_bytes.size()));
    // owner_p2.send(span<uint8_t>(share2_bytes.data(), share2_bytes.size()));
    // owner_p3.send(span<uint8_t>(share3_bytes.data(), share3_bytes.size()));

    // // sharing shares to each party
    // for (size_t i = 0; i < triples.size(); ++i)
    // {
    //     // Determine which element to send
    //     BitVector current = triples[i];
    //     BitVector next = triples[(i + 1) % triples.size()];

    //     // Create a vector of size 2 containing the elements
    //     // std::vector<BitVector> v{current, next};
    //     // Send the vector to the corresponding party
    //     if (i == 0)
    //     {
    //         // Send vector[0], vector[1] to party1
    //         owner_p1.send(current);
    //         owner_p1.send(next);
    //     }
    //     else if (i == 1)
    //     {
    //         // Send vector[1], vector[2] to party2
    //         owner_p2.send(current);
    //         owner_p2.send(next);
    //     }
    //     else if (i == 2)
    //     {
    //         // Send vector[2], vector[0] to party3
    //         owner_p3.send(current);
    //         owner_p3.send(next);
    //     }
    // }
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
void runMPC(int pid, const string ownerIP, int ownerPort, const string address1, int recPort1, int sendPort1, const string address2, int recPort2, int sendPort2)
{
    mpc testmpc;
    std::promise<void> promise;
    std::future<void> future = promise.get_future();
    testmpc.initialize(pid, ownerIP, ownerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    testmpc.receiveSecrets();
    testmpc.Frand(5);
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
    // t1.join();
    // t2.join();
    // t3.join();

    return 0;
}
