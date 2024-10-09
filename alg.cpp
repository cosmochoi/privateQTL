/*
Sorting algorithm implementation based on
"Efficient Secure Three-Party Sorting with Applications to Data Analysis and Heavy Hitters"
by Gilad Asharov et al.
*/

#include <iostream>
#include <bitset>
#include "alg.h"
#include <string>
#include "mpc.h"

// reshare: Given shares of secret, mask shares with random value and return shares to receiver
void reshare(Party *senderA, Party *senderB, int secretID, Party *receiver)
{
    vector<BitVector> senderAShare = senderA->retrieveShares(secretID).replShares;
    vector<BitVector> senderBShare = senderB->retrieveShares(secretID).replShares;
    vector<BitVector> receiverShares(2);
    auto it = senderA->seedpair.find(senderB->pid);
    if (it != senderA->seedpair.end())
    {
        vector<BitVector> r(3);
        PRNG &sharedPRNG = *(it->second);
        size_t size = senderAShare[0].size();
        r[0].resize(size);
        r[1].resize(size);
        r[2].resize(size);
        sharedPRNG.get((u8 *)r[0].data(), size);
        sharedPRNG.get((u8 *)r[1].data(), size);
        r[2] ^= r[0];
        r[2] ^= r[1];
        receiverShares[0] = senderAShare[1] ^ r[(senderA->pid + 1) % 3];
        receiverShares[1] = senderBShare[0] ^ r[senderB->pid];
        receiver->retrieveShares(secretID).replShares = receiverShares;
    }
    else
    {
        throw invalid_argument("partner and sender do not share PRNG");
    }
}

vector<int> Fmult(vector<int> keyvector, vector<vector<int>> sumvector)
{
    vector<int> permutation(keyvector.size());
    for (int i = 0; i < keyvector.size(); i++)
    {
        permutation[i] = sumvector[0][i] + keyvector[i] * (sumvector[1][i] - sumvector[0][i]);
    }
    return permutation;
}

// Frand: generate random permutation
vector<uint64_t> Frand(PRNG *prng, uint64_t bufferSize)
{
    if (bufferSize <= 1)
    {
        throw std::invalid_argument("bufferSize must be positive");
    }

    vector<uint64_t> pi(bufferSize);
    for (uint64_t i = 0; i < bufferSize; i++)
    {
        pi[i] = i + 1;
    }
    for (uint64_t i = bufferSize - 1; i >= 1; i--)
    {
        uint64_t j = prng->get<uint64_t>() % (i + 1);
        uint64_t temp = pi[i];
        pi[i] = pi[j];
        pi[j] = temp;
    }
    return pi;
}

// GENBITPERM: sorting algorithm for one bit vector
vector<int> genbitperm(vector<int> unsorted)
{
    vector<int> permutation(unsorted.size());
    vector<vector<int>> f(2, vector<int>(unsorted.size()));
    vector<vector<int>> s(2, vector<int>(unsorted.size()));
    for (int i = 0; i < unsorted.size(); i++)
    {
        f[0][i] = 1 - unsorted[i];
        f[1][i] = unsorted[i];
    }
    char prefix_sum = 0;
    for (int j = 0; j < 2; j++)
    {
        for (int i = 0; i < unsorted.size(); i++)
        {
            prefix_sum = prefix_sum + f[j][i];
            s[j][i] = prefix_sum;
        }
    }
    permutation = Fmult(unsorted, s);
    return permutation;
}

// Apply permutation
template <typename T>
void apply_permutation(std::vector<T> &v, const std::vector<int> &indices)
{
    std::vector<T> v2;
    v2.reserve(v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        v2.push_back(v[indices[i] - 1]);
    }
    v = std::move(v2);
}

int main()
{
    int test[4] = {1, 2, 3, 4};

    vector<string> bitdecomp;
    vector<int> perm;

    // vector<char> binary;
    for (int i = 0; i < 4; i++)
    {
        bitdecomp.push_back(bitset<8>(test[i]).to_string());
    }
    // vector<int> bittest{1, 1, 0, 0};
    // perm = genbitperm(bittest);
    // apply_permutation(bittest, perm);
    // print_vector(bittest);

    vector<BitVector> keys(4);
    for (int j = 0; j < 4; j++)
    {
        keys[j] = BitVector(bitdecomp[j]);
        // std::cout << keys[j] << std::endl;
    }
    mpc testmpc(3);
    vector<BitVector> testshare(3);
    testshare = testmpc.createShares(keys[1]);

    cout << "secret: " << keys[1] << endl;
    print_vector(testshare);
    cout << (testshare[0] ^ testshare[1] ^ testshare[2]) << endl;
    testmpc.distributeShares(testshare, testmpc.parties);
    // for (int i = 0; i < testmpc.parties.size(); i++)
    // {
    //     cout << "Party " << testmpc.parties[i].pid << " shares : " << testmpc.parties[i].numShares << endl;
    // }
    auto ittt = testmpc.parties[0].seedpair.find(testmpc.parties[1].pid);
    if (ittt != testmpc.parties[0].seedpair.end())
    {

        PRNG &testPRNG = *(ittt->second);
        cout << keys.size() << endl;
        if (&testPRNG == nullptr)
        {
            cout << "testPRNG null" << endl;
        }
        else
        {
            vector<uint64_t> randomperm = Frand(&testPRNG, 4);
            print_vector(randomperm);
        }
    }
    else
    {
        cout << "couldn't find seedpair" << endl;
    }
    return 0;
}