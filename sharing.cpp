#include <iostream>
#include <thread>
#include <bitset>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Defines.h>

using namespace osuCrypto;
using namespace std;

//<<move to h file>>
template <typename T>
void print_vector(vector<T> printme)
{
    for (int j = 0; j < printme.size(); j++)
    {
        std::cout << printme[j] << " ";
    }
    std::cout << "\n";
}

// Define the function to be computed
bool computeFunction(bool x1, bool x2, bool x3)
{
    // Example function: output 1 if exactly two inputs are true, 0 otherwise
    return (x1 + x2 + x3) == 2;
}
void dataclient(uint64_t input, int sendport1, std::string address1, int sendport2, std::string address2, int sendport3, std::string address3)
{
    IOService ios;
    Endpoint epsend1(ios, address1, sendport1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendport2, EpMode::Server);
    Endpoint epsend3(ios, address3, sendport3, EpMode::Server);
    Channel owner_p1 = epsend1.addChannel();
    Channel owner_p2 = epsend2.addChannel();
    Channel owner_p3 = epsend3.addChannel();
    std::cout << "Established channels with the computing parties.\n";
    PRNG ownerprng;
    ownerprng.SetSeed(toBlock(27));
    std::string sinput = std::bitset<64>(input).to_string();
    BitVector decomposed = BitVector(sinput);
    // making the shares
    // The last share should be xor'd with all other bitvector
    std::vector<BitVector> triples(3);
    triples[2].assign(decomposed);
    for (int i = 0; i < 2; i++)
    {
        triples[i].assign(sinput);
        triples[i].randomize(ownerprng);
        triples[2] ^= triples[i];
    }

    // sharing shares to each party
    for (size_t i = 0; i < triples.size(); ++i)
    {
        // Determine which element to send
        BitVector current = triples[i];
        BitVector next = triples[(i + 1) % triples.size()];

        // Create a vector of size 2 containing the elements
        // std::vector<BitVector> v{current, next};
        // Send the vector to the corresponding party
        if (i == 0)
        {
            // Send vector[0], vector[1] to party1
            owner_p1.send(current);
            owner_p1.send(next);
        }
        else if (i == 1)
        {
            // Send vector[1], vector[2] to party2
            owner_p2.send(current);
            owner_p2.send(next);
        }
        else if (i == 2)
        {
            // Send vector[2], vector[0] to party3
            owner_p3.send(current);
            owner_p3.send(next);
        }
    }
    owner_p1.close();
    owner_p2.close();
    owner_p3.close();
    std::cout << "Sent secret shared values to parties.\n";
}

void party(int id, int ownerport, int recPort1, int recPort2, bool input, std::string address1, int sendport1, std::string address2, int sendport2)
{
    // Initialize the IOService and create channels to the other parties
    IOService ios;
    Endpoint p_owner(ios, "localhost", ownerport, EpMode::Client);
    Endpoint eprec1(ios, "localhost", recPort1, EpMode::Client);
    Endpoint eprec2(ios, "localhost", recPort2, EpMode::Client);
    Endpoint epsend1(ios, address1, sendport1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendport2, EpMode::Server);
    Channel dataowner = p_owner.addChannel();
    Channel chl = eprec1.addChannel();
    Channel chl1 = eprec2.addChannel();
    Channel chl2 = epsend1.addChannel();
    Channel chl3 = epsend2.addChannel();
    try
    {
        // Receive input from data owner
        // std::vector<u8> share1;
        // std::vector<u8> share2;
        BitVector share1 = BitVector();
        BitVector share2 = BitVector();
        std::cout << "ready to receive:\n";
        dataowner.recv(share1);
        std::cout << "received:\n";
        dataowner.recv(share2);
        vector<BitVector> recshares{share1, share2};
        std::cout << "Party " << id << " received shares:\n";
        print_vector(recshares);
        // BitVector shared1 = BitVector(share1.data(), share1.size());
        // BitVector shared2 = BitVector(share2.data(), share2.size());
        // print_vector(shares);
    }
    catch (const std::exception &e)
    {
        std::cout << "except\n";
        std::cout << e.what() << std::endl;
    }

    // // Share the input with the other parties
    // uint8_t inputBuffer[1];
    // inputBuffer[0] = input;
    // chl2.asyncSend(inputBuffer, 1);
    // chl3.asyncSend(inputBuffer, 1);

    // // Receive the inputs from the other parties
    // uint8_t inputBuffers[2][1];
    // chl.recv(inputBuffers[0], 1);
    // chl1.recv(inputBuffers[1], 1);

    // // Compute the function on the inputs
    // bool output = computeFunction(input, inputBuffers[0][0], inputBuffers[1][0]);
    // std::cout << output << "ffff" << std::endl;
    // // Share the output with the other parties
    // uint8_t outputBuffer[1];
    // outputBuffer[0] = output;
    // chl2.send(outputBuffer, 1);
    // chl3.send(outputBuffer, 1);

    // Close the channels
    chl.close();
    chl1.close();
    chl2.close();
    chl3.close();
}

int main()
{
    int port12 = 12345; // port number for party 1
    int port13 = 12346; // port number for party 2
    int port21 = 12355; // port number for party 3
    int port23 = 12356; // port number for party 1
    int port31 = 12365; // port number for party 2
    int port32 = 12366; // port number for party 3
    int owner1 = 12375; // port number for dataowner
    int owner2 = 12376;
    int owner3 = 12377;
    std::string address1 = "localhost"; // address of the remote ssh server for party 1
    std::string address2 = "localhost"; // address of the remote ssh server for party 2
    std::string address3 = "localhost"; // address of the remote ssh server for party 3

    // Run each party as a separate thread
    std::thread t1(dataclient, 4, owner1, address1, owner2, address2, owner3, address3);
    std::thread t2(party, 0, owner1, port12, port13, false, address2, port21, address3, port31);
    std::thread t3(party, 1, owner2, port21, port23, false, address3, port12, address1, port32);
    std::thread t4(party, 2, owner3, port31, port32, true, address1, port13, address2, port23);

    t1.join();
    t2.join();
    t3.join();
    t4.join();

    return 0;
}
