#ifndef MPC_H
#define MPC_H
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Defines.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <map>
#include <set>
#include <bitset>
#include <fstream>

using namespace osuCrypto;
using namespace std;
using namespace NTL;

// struct Share;
// struct Party;

// struct Share
// {
//     Party *owner;
//     int secret_id;
//     vector<BitVector> replShares; //(a1, a2)
//     Share(Party *o, int secret_id)
//     {
//         owner = o;
//         this->secret_id = secret_id;
//     }
// };

// struct Party
// {
//     int pid;
//     set<Share *> shares;
//     int numShares = 0;
//     map<int, PRNG *> seedpair;
//     explicit Party(int PID)
//     {
//         this->pid = PID;
//     };
//     void addShare(Share *share)
//     {
//         this->shares.insert(share);
//         numShares++;
//     }
//     Share retrieveShares(int sID)
//     {
//         auto it = find_if(shares.begin(), shares.end(),
//                           [sID](const Share *s)
//                           { return s->secret_id == sID; });
//     }
// };

// extern vector<double> CSVtoVector(string filename);
// extern vector<uint64_t> ScaleVector(vector<double> &v, int k);

class mpc
{
public:
    IOService ios;
    PRNG globalprng;
    int pid;
    uint32_t inv;
    mpc(){};
    bool initialize(int pid, string ownerIP, int ownerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupChannels(string ownerIP, int ownerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupSeeds(); // shared seeds for PRNG
    void receiveSecrets();
    vector<vector<uint64_t>> Frand(uint64_t bufferSize);
    vector<vector<uint64_t>> reveal(vector<vector<uint64_t>> pi);
    void reshare(vector<BitVector>& shares, int reshareID);
    void reshare(vector<BitVector>& shares);
    vector<uint64_t> Fmult(vector<uint64_t> k_i, vector<uint64_t> s_i);
    // vector<vector<uint64_t>> genbitperm(vector<vector<uint64_t>> keybit);
    void close();

private:
    block commonSeed = oc::toBlock(27);
    map<int, PRNG *> seedpair;
    vector<vector<BitVector>> shares;
    Channel dataowner;
    Channel fromPlus;
    Channel fromMinus;
    Channel toPlus;
    Channel toMinus;
};

template <typename T>
void print_vector(vector<T> printme)
{
    for (int j = 0; j < printme.size(); j++)
    {
        std::cout << printme[j] << " ";
    }
    std::cout << "\n";
}

template <typename T>
void print_vector(vector<vector<T>> printme)
{
    for (int i =0; i<printme.size(); i++)
    {
        print_vector(printme[i]);
    }
    std::cout << "\n";
}

#endif