#include "mpc.h"
#include "utils.h"

bool mpc::initialize(int pid, string ownerIP, int ownerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    this->pid = pid;
    globalprng.SetSeed(commonSeed);
    if (!setupChannels(ownerIP, ownerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2))
    {
        cout << "Network channel failed.\n";
    }
    if (!setupSeeds())
    {
        cout << "PRNG setup failed.\n";
    }
    return true;
}
bool mpc::setupChannels(string ownerIP, int ownerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    Endpoint p_owner(ios, ownerIP, ownerPort, EpMode::Client);
    Endpoint eprec1(ios, address1, recPort1, EpMode::Client);
    Endpoint eprec2(ios, address2, recPort2, EpMode::Client);
    Endpoint epsend1(ios, address1, sendPort1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendPort2, EpMode::Server);
    this->dataowner = p_owner.addChannel();
    this->fromPlus = eprec1.addChannel();
    this->fromMinus = eprec2.addChannel();
    this->toPlus = epsend1.addChannel();
    this->toMinus = epsend2.addChannel();
    cout << "Established channels with the computing parties and data owner.\n";
    p_owner.stop();
    eprec1.stop();
    eprec2.stop();
    epsend1.stop();
    epsend2.stop();
    return true;
}
bool mpc::setupSeeds()
{
    // cout << "setting up shared seeds within parties..\n";
    // lower number pid party sends prng seed to higher number pid party
    //  Check if they already share a prng
    auto it = this->seedpair.find((this->pid + 1) % 3);
    auto it2 = this->seedpair.find((this->pid + 2) % 3);
    if (it != this->seedpair.end())
    {
        cout << "Party " << this->pid << " already share a prng with party " << (this->pid + 1) % 3 << endl;
        return false;
    }
    if (it2 != this->seedpair.end())
    {
        cout << "Party " << this->pid << " already share a prng with party " << (this->pid + 2) % 3 << endl;
        return false;
    }
    try
    {
        PRNG prng(oc::toBlock(this->pid)); //initializing private prng with pid..?
        uint64_t plusSeed = prng.get<uint64_t>();
        vector<uint64_t> buffer(sizeof(uint64_t));
        memcpy(buffer.data(), &plusSeed, sizeof(uint64_t));
        this->toPlus.asyncSend(buffer.data(), buffer.size());

        vector<uint64_t> dest;
        this->fromMinus.recv(dest);
        // retrieve uint64_t from buffer
        uint64_t minusSeed = dest[0];

        // PRNG wPlus(oc::toBlock(plusSeed));
        // PRNG wMinus(oc::toBlock(minusSeed));
        PRNG *wPlus = new PRNG(oc::toBlock(plusSeed));
        PRNG *wMinus = new PRNG(oc::toBlock(minusSeed));

        /// add prng to each party map
        this->seedpair.insert({(this->pid + 1) % 3, wPlus});
        this->seedpair.insert({(this->pid + 2) % 3, wMinus});
    }
    catch (const std::exception &e)
    {
        cout << "except\n";
        cout << e.what() << std::endl;
    }

    cout << "Party " << (this->pid + 1) % 3 << " and Party " << (this->pid + 2) % 3 << " seed setup complete." << endl;
    return true;
}

// Convert a one-dimensional byte vector into a vector<vector<BitVector>> object
//might have to DELETE; double check
std::vector<std::vector<BitVector>> unflatten(const std::vector<uint8_t>& data, size_t num_rows, size_t num_cols) {
    std::vector<std::vector<BitVector>> result(num_rows, std::vector<BitVector>(num_cols));
    size_t pos = 0;
    for (size_t i = 0; i < num_rows; i++) {
        for (size_t j = 0; j < num_cols; j++) {
            size_t bitvec_size = result[i][j].sizeBytes();
            std::memcpy(result[i][j].data(), data.data() + pos, bitvec_size);
            pos += bitvec_size;
        }
    }
    return result;
}

void mpc::receiveSecrets()
{
    vector<int> dest;
    this->dataowner.recv(dest);
    // cout << "Size of input: " << dest[0] << endl;
    vector<vector<BitVector>> result;

    try {
    for (int i = 0; i < dest[0]; i++) {
        vector<BitVector> resultRow;
        for (int j = 0; j < 2; j++) {
            BitVector received_data;
            this->dataowner.recv(received_data);
            resultRow.push_back(received_data);
        }
        result.push_back(resultRow);
    }

        this->shares = result;
        // cout << "Final size of vector " << this->shares.size() << endl;
        // print_vector(result);
    }
    catch (const exception &e)
    {
        cout << "except\n";
        cout << e.what() << endl;
    }
    cout << this->shares[0][0][0] << endl;
    cout << this->shares[0][1][0] << endl;
    vector<uint64_t> ki{this->shares[0][0][0],this->shares[0][1][0]};
    vector<uint64_t> mult;
    mult = Fmult(ki, ki);
    print_vector(mult);
}

vector<vector<uint64_t>> mpc::Frand(uint64_t bufferSize) 
{
    if (bufferSize <= 1)
    {
        throw std::invalid_argument("bufferSize must be positive");
    }
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &plus_PRNG = *(plus_it->second);
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    vector<vector<uint64_t>> pi(2,vector<uint64_t>(bufferSize));
    for (uint64_t i = 0; i < bufferSize; i++)
    {
        pi[0][i] = i + 1;
    }
    for (uint64_t i = bufferSize - 1; i >= 1; i--)
    {
        uint64_t j = minus_PRNG.get<uint64_t>() % (i + 1);
        uint64_t temp = pi[0][i];
        pi[0][i] = pi[0][j];
        pi[0][j] = temp;
    }
    for (uint64_t i = 0; i < bufferSize; i++)
    {
        pi[1][i] = i + 1;
    }
    for (uint64_t i = bufferSize - 1; i >= 1; i--)
    {
        uint64_t j = plus_PRNG.get<uint64_t>() % (i + 1);
        uint64_t temp = pi[1][i];
        pi[1][i] = pi[1][j];
        pi[1][j] = temp;
    }
    // print_vector(pi);
    return pi;
}

vector<vector<uint64_t>> mpc::reveal(vector<vector<uint64_t>> pi)
{
    this->toPlus.send(pi[0]);
    this->toMinus.send(pi[1]);
    vector<uint64_t> receivedMinus, receivedPlus;
    this->fromMinus.recv(receivedMinus);
    this->fromPlus.recv(receivedPlus);
    vector<vector<uint64_t>> reconstructed;
    if (receivedMinus != receivedPlus) {
        throw string("received inputs do not match.");
    }
    else {
        reconstructed.push_back(receivedMinus);
        reconstructed.push_back(pi[0]);
        reconstructed.push_back(pi[1]);
    }
    return reconstructed;
}

void mpc::reshare(vector<BitVector>& shares, int reshareID)
{
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &plus_PRNG = *(plus_it->second);
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    BitVector r_i, r_iplus;
    for (uint32_t i = 0; i < shares[0].size(); i++) {
        bool randomBit = minus_PRNG.get<bool>();  
        r_i[i] = randomBit;     
    }
    for (uint32_t i = 0; i < shares[1].size(); i++) {
        bool randomBit = plus_PRNG.get<bool>(); 
        r_iplus[i] = randomBit;     
    }
    vector<BitVector> newshare{shares[0]^r_i, shares[1]^r_iplus};
    if (reshareID == (this->pid + 1) % 3)
    {
        this->toPlus.send(shares[1]^r_iplus);
    }
    else 
    {
        this->toMinus.send(shares[0]^r_i);
    }
    shares = newshare;
}

void mpc::reshare(vector<BitVector>& shares)
{
    // vector<BitVector> newshare(2);
    this->fromMinus.recv(shares[0]);
    this->fromPlus.recv(shares[1]);
    // shares = newshare;
}

//is adding alpha beta gamma necessary? 
vector<uint64_t> mpc::Fmult(vector<uint64_t> k_i, vector<uint64_t> s_i)
{
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    vector<uint64_t> t_i(2);
    uint64_t ri = k_i[0]*s_i[0] + k_i[1]*s_i[0] + k_i[0]*s_i[1];
    t_i[0]=ri;
    this->toPlus.send(ri);
    this->fromMinus.recv(t_i[1]);
    return t_i;
}

vector<vector<uint64_t>> mpc::genbitperm(vector<vector<uint64_t>> keybit)
{
    vector<vector<uint64_t>> f0(keybit.size(), vector<uint64_t>(2));
    vector<vector<uint64_t>> f1(keybit.size(), vector<uint64_t>(2));
    // vector<vector<uint64_t>> s(keybit.size(), vector<uint64_t>(2));
    for (int i=0; i<keybit.size(); i++) 
    {
        vector<uint64_t> ones(keybit[i].size(), 1);
        for (int j=0; j<keybit[i].size(); j++) 
        {
            f0[i][j] = ones[j]-keybit[i][j];
            f1[i][j] = keybit[i][j];
        }
        
    }

    return keybit;
}

void mpc::close()
{
    cout << "Closing channels.\n";
    this->dataowner.close();
    this->fromPlus.close();
    this->fromMinus.close();
    this->toPlus.close();
    this->toMinus.close();
    this->ios.stop();
    cout << "channels closed.\n";
}
