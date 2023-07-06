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
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm> //for debugging
#include <sys/socket.h>
#include <netinet/in.h>
#include <cstring>

using namespace osuCrypto;
using namespace std;
using namespace NTL;

class mpc
{

public:
    IOService ios;
    PRNG globalprng;
    int pid;
    uint32_t inv;
    uint32_t lk;
    uint32_t n;
    uint32_t p;
    double shiftsize;
    vector<double> originalData; // TODO: Remove this!
    // reference_wrapper<std::atomic<int>> readyCounter;
    // reference_wrapper<std::mutex> mtx;
    // reference_wrapper<std::condition_variable> cv;
    // mpc(atomic<int>& counter, mutex& mtx, condition_variable& cv)
    //     : readyCounter(counter), mtx(mtx), cv(cv) {};
    mpc(){};
    bool initialize(int pid, string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupChannels(string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupSeeds(); // shared seeds for PRNG
    bool areVectorsEmpty() const {
        return shares.empty() && identity.empty() && zscores.empty();
    }
    void ready();
    void receiveSecrets();
    void assertSize(vector<ZZ_p> pi, string tag);
    void assertValid(vector<ZZ_p> sigma, string tag);
    vector<uint32_t> apply_plaintext_perm(vector<uint32_t> rho, vector<uint32_t> sigma);
    vector<ZZ_p> Frand(uint32_t bufferSize);
    vector<ZZ_p> reveal(vector<ZZ_p> &pi, bool isperm);
    void reshare(vector<ZZ_p> &shares, int reshareID);
    // void reshare(vector<ZZ_p> &shares);
    vector<ZZ_p> Fmult(vector<ZZ_p> k_i, vector<ZZ_p> s_i);
    vector<ZZ_p> genbitperm(vector<ZZ_p> &keybit);
    vector<ZZ_p> inversePerm(vector<ZZ_p> pi);
    void apply_perm_local(bool participate,vector<ZZ_p> &v, vector<ZZ_p> &pi);
    void shuffle(vector<ZZ_p> &pi, vector<ZZ_p> &a);
    void unshuffle(vector<ZZ_p> &pi, vector<ZZ_p> &b);
    void apply_shared_perm(vector<ZZ_p> &rho, vector<ZZ_p> &k);
    void compose(vector<ZZ_p> &sigma, vector<ZZ_p> &rho);
    // vector<ZZ_p> get_shared_inverse(vector<ZZ_p> sigma);
    vector<double> genperm(int row, int numCol);
    void clearVectors() {
        shares.clear();
        identity.clear();
        zscores.clear();
    }
    void close();

private:
    block commonSeed = oc::toBlock(27);
    map<int, PRNG *> seedpair;
    vector<ZZ_p> shares;
    vector<ZZ_p> identity;
    vector<ZZ_p> zscores;
    // atomic<int>& readyCounter;
    // atomic<bool> allPartiesReady;
    Channel dataowner;
    Channel toOwner;
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

template <typename T>
void assert_equal(vector<T> one, vector<T> two, string tag)
{
    if (one.size() != two.size()) {
        throw logic_error("[" + tag + "]: vectors are not equal length");
    }
    for (int i = 0; i < one.size(); i++)
    {
        if (one[i] != two[i]) {
            throw logic_error("[" + tag + "]: vectors are not equal");
        }
    }
}

// vector<ZZ_p> convVec(vector<uint32_t> v){
//     vector<ZZ_p> converted(v.size());
//     for (int i=0; i<v.size(); i++)
//     {
//         converted[i] = conv<ZZ_p>(v[i]);
//     }
//     return converted;
// }
// vector<uint32_t> convVec(vector<ZZ_p> v){
//     vector<uint32_t> converted(v.size());
//     for (int i=0; i<v.size(); i++)
//     {
//         converted[i] = conv<uint32_t>(v[i]);
//     }
//     return converted;
// }
#endif