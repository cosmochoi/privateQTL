#ifndef MPC_H
#define MPC_H
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Defines.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <map>
#include <set>
#include <bitset>
#include <fstream>

using namespace osuCrypto;
using namespace std;
using namespace NTL;

class mpc
{
public:
    IOService ios;
    PRNG globalprng;
    int pid;
    uint64_t inv;
    uint64_t lk;
    uint64_t n;
    mpc(){};
    bool initialize(int pid, string ownerIP, int ownerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupChannels(string ownerIP, int ownerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupSeeds(); // shared seeds for PRNG
    void receiveSecrets();
    vector<ZZ_p> Frand(uint64_t bufferSize);
    vector<ZZ_p> reveal(vector<ZZ_p> &pi);
    void reshare(vector<ZZ_p> &shares, int reshareID);
    void reshare(vector<ZZ_p> &shares);
    vector<ZZ_p> Fmult(vector<ZZ_p> k_i, vector<ZZ_p> s_i);
    vector<ZZ_p> genbitperm(vector<ZZ_p> &keybit);
    void apply_perm_local(vector<ZZ_p> &v, vector<ZZ_p> &pi);
    void shuffle(vector<ZZ_p> &pi, vector<ZZ_p> &a);
    void unshuffle(vector<ZZ_p> &pi, vector<ZZ_p> &b);
    void apply_shared_perm(vector<ZZ_p> &rho, vector<ZZ_p> &k);
    void compose(vector<ZZ_p> &sigma, vector<ZZ_p> &rho);
    void genperm(vector<ZZ_p> k);
    void close();

private:
    block commonSeed = oc::toBlock(27);
    map<int, PRNG *> seedpair;
    vector<ZZ_p> shares;
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
template <typename T>
void print_vector(Vec<T> printme)
{
    for (int i =0; i<printme.length(); i++)
    {
        std::cout << printme.get(i) << " ";
    }
    std::cout << "\n";
}


#endif