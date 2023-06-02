#include "mpc.h"
#include "utils.h"

bool mpc::initialize(int pid, string ownerIP, int ownerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    // ZZ_p::init(to_ZZ(4));
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

    return true;
}

void mpc::receiveSecrets()
{
    vector<uint64_t> dest;
    this->dataowner.recv(dest);
    ZZ_p::init(to_ZZ(dest[2])); ///p
    this->n = dest[0]; //number of secrets
    this->lk = dest[1]; //number of bits
    this->inv = dest[3]; //inv 
    vector<uint64_t> result;
    try 
    {
        for (int i = 0; i < dest[0]*dest[1]*2; i++) 
        {
            uint64_t received_data;
            this->dataowner.recv(received_data);
            result.push_back(received_data);
        }
        vector<ZZ_p> temp(result.begin(), result.end());
        swap(this->shares, temp);
    }
    catch (const exception &e)
    {
        cout << "except\n";
        cout << e.what() << endl;
    }
}

vector<ZZ_p> mpc::Frand(uint64_t bufferSize) 
{
    if (bufferSize <= 1)
    {
        throw std::invalid_argument("bufferSize must be positive");
    }
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &plus_PRNG = *(plus_it->second);
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    vector<uint64_t> pi(2*bufferSize);
    for (uint64_t i = 0; i < bufferSize; i++)
    {
        pi[i*2+0] = i + 1;
        pi[i*2+1] = i + 1;
    }
    for (uint64_t i = bufferSize - 1; i >= 1; i--)
    {
        uint64_t j = minus_PRNG.get<uint64_t>() % (i + 1);
        uint64_t temp = pi[i*2+0];
        pi[i*2+0] = pi[j*2+0];
        pi[j*2+0] = temp;

        uint64_t k = plus_PRNG.get<uint64_t>() % (i + 1);
        uint64_t temp2 = pi[i*2+1];
        pi[i*2+1] = pi[k*2+1];
        pi[k*2+1] = temp2;
    }

    vector<ZZ_p> ZZ_pi(pi.begin(), pi.end());

    return ZZ_pi;
}

vector<ZZ_p> mpc::reveal(vector<ZZ_p>& pi)
{
    vector<uint64_t> share1;
    vector<uint64_t> share2;
    for (int i=0; i<pi.size(); i+=2)
    {
        share1.push_back(conv<uint64_t>(pi[i]));
        share2.push_back(conv<uint64_t>(pi[i+1]));
    }
    this->toPlus.send(share1);
    this->toMinus.send(share2);
    vector<uint64_t> receivedMinus, receivedPlus;
    this->fromMinus.recv(receivedMinus);
    this->fromPlus.recv(receivedPlus);
    vector<ZZ_p> reconstructed;
    if (receivedMinus != receivedPlus) 
    {
        throw string("received inputs do not match.");
    }
    else 
    {
        for (int i=0; i<receivedMinus.size(); i++)
        {
            reconstructed.push_back(conv<ZZ_p>(receivedMinus[i])+conv<ZZ_p>(share1[i])+conv<ZZ_p>(share2[i]));
        }
    }
    return reconstructed;
}

void mpc::reshare(vector<ZZ_p>& shares, int reshareID)
{
    vector<ZZ_p> randoms(3);
    if (reshareID == (this->pid + 1) % 3) // share to plus
    {
        auto minus_it = this->seedpair.find((this->pid + 2) % 3);
        PRNG &minus_PRNG = *(minus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(minus_PRNG.get<uint32_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int j=0; j<shares.size()/2; j++)
        {
            shares[2*j]+=randoms[this->pid];
            shares[2*j+1]+=randoms[(this->pid +1)%3];
            this->toPlus.send(conv<uint64_t>(shares[2*j+1]));
        }
    }
    else  //share to minus
    {
        auto plus_it = this->seedpair.find((this->pid + 1) % 3);
        PRNG &plus_PRNG = *(plus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(plus_PRNG.get<uint32_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int j=0; j<shares.size()/2; j++)
        {
            shares[2*j]+=randoms[this->pid];
            shares[2*j+1]+=randoms[(this->pid +1)%3];
            this->toMinus.send(conv<uint64_t>(shares[2*j]));
        }
    }
}

void mpc::reshare(vector<ZZ_p> &shares)
{
    vector<uint64_t> newshare(shares.size());
    for (int i=0; i<shares.size()/2; i++)
    {
        this->fromMinus.recv(newshare[2*i]);
        this->fromPlus.recv(newshare[2*i+1]);
    }
    vector<ZZ_p> temp(newshare.begin(), newshare.end());
    shares.swap(temp);
}

//is adding alpha beta gamma necessary? 
vector<ZZ_p> mpc::Fmult(vector<ZZ_p> k_i, vector<ZZ_p> s_i)
{
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    vector<ZZ_p> t_i(2);
    ZZ_p ri = k_i[0]*s_i[0] + k_i[1]*s_i[0] + k_i[0]*s_i[1];
    t_i[0]=ri;
    this->toMinus.send(conv<int>(ri));
    int received_data;
    this->fromPlus.recv(received_data);
    conv(t_i[1], received_data);
    return t_i;
}

vector<ZZ_p> mpc::genbitperm(vector<ZZ_p> &keybit)
{
    vector<ZZ_p> f0(keybit.size()), f1(keybit.size()), s0(keybit.size()), s1(keybit.size());
    vector<ZZ_p> s(2, to_ZZ_p(0));
    vector<ZZ_p> p;
    f1 = keybit;
    for (int i=0; i<keybit.size()/2; i++)
    {
        f0[2*i] = this->inv-keybit[2*i];
        f0[2*i+1] = this->inv-keybit[2*i+1];
        s[0]+=f0[2*i];
        s[1]+=f0[2*i+1];
        s0[2*i] = s[0];
        s0[2*i+1] = s[1];
    }
    for (int i=0; i<keybit.size()/2; i++)
    {
        s[0]+=f1[2*i];
        s[1]+=f1[2*i+1];
        s1[2*i] = s[0];
        s1[2*i+1] = s[1];
    }
    vector<ZZ_p> t;
    for (int i=0; i <keybit.size()/2; i++)
    {
        vector<ZZ_p> sminus(2);
        vector<ZZ_p> ti(2);
        sminus[0] = s1[2*i] - s0[2*i];
        sminus[1] = s1[2*i+1] - s0[2*i+1];
        vector<ZZ_p> ki(keybit.begin()+2*i, keybit.begin()+2*(i+1));
        ti = Fmult(ki, sminus);
        t.insert(t.end(), ti.begin(), ti.end());
    }
    for (int m=0; m<keybit.size(); m++)
    {
        p.push_back(s0[m]+t[m]);
    }
    return p;
}

vector<ZZ_p> inversePerm(vector<ZZ_p> pi)
{
    vector<uint64_t> arr2(pi.size());
    vector<ZZ_p> inverse(pi.size());
    for (int i = 0; i < pi.size(); i++)
        arr2[conv<int>(pi[i]) - 1] = i + 1;
 
    for (int i = 0; i < pi.size(); i++)
    {
        inverse[i] = conv<ZZ_p>(arr2[i]);
    }
    return inverse;
}

void mpc::apply_perm_local(vector<ZZ_p> &v, vector<ZZ_p> &pi)
{
    if (pi.size() != v.size()/2)
        throw std::invalid_argument("local permutation should be half the size of shared vector");
    vector<ZZ_p> v2(v.size());
    for (size_t i = 0; i < v.size()/2; i++)
    {
        uint64_t pi_i = conv<uint64_t>(pi[i]);
        // uint64_t pi2 = conv<uint64_t>(pi[2*i+1]);
        v2[(pi_i-1)*2] =v[2*i];
        v2[(pi_i-1)*2+1]=v[2*i+1];
    }

    vector<ZZ_p> testing;
    testing = inversePerm(pi);
    v = std::move(v2);
}
void mpc::shuffle(vector<ZZ_p> &pi, vector<ZZ_p> &a)
{
    if (pi.size() != a.size())
        throw std::invalid_argument("Your shares and pi size are different");
    vector<ZZ_p> pi_m(pi.size()/2), pi_p(pi.size()/2);
    for (int i=0; i<pi.size()/2; i++)
    {
        pi_m[i] = pi[2*i];
        pi_p[i] = pi[2*i+1];
    }

    if (this->pid == 0)
    {
        apply_perm_local(a, pi_m);
        reshare(a, 1);
        apply_perm_local(a, pi_p);
        reshare(a, 2);
        reshare(a);
    }
    else if (this->pid == 1)
    {
        reshare(a);
        apply_perm_local(a, pi_m);
        reshare(a, 2);
        apply_perm_local(a, pi_p);
        reshare(a, 0);
    }
    else
    {
        apply_perm_local(a, pi_p);
        reshare(a, 1);
        reshare(a);
        apply_perm_local(a, pi_m);
        reshare(a, 0);
    }
}

void mpc::unshuffle(vector<ZZ_p> &pi, vector<ZZ_p> &b)
{
    vector<ZZ_p> pi_m, pi_p;
    for (size_t i=0; i<pi.size()/2; i++)
    {
        pi_m.push_back(pi[2*i]);
        pi_p.push_back(pi[2*i+1]);
    }
    vector<ZZ_p> inv_m(pi_m.size()), inv_p(pi_p.size());
    inv_m = inversePerm(pi_m);
    inv_p = inversePerm(pi_p); 
    if (this->pid == 0)
    {
        reshare(b);
        apply_perm_local(b,inv_p);
        reshare(b,2);
        apply_perm_local(b,inv_m);
        reshare(b,1);
    }
    else if (this->pid == 1)
    {
        apply_perm_local(b,inv_p);
        reshare(b, 0);
        apply_perm_local(b,inv_m);
        reshare(b, 2);
        reshare(b);
    }
    else
    {
        apply_perm_local(b,inv_m);
        reshare(b, 0);
        reshare(b);
        apply_perm_local(b,inv_p);
        reshare(b, 1);
    }
}

void mpc::apply_shared_perm(vector<ZZ_p> &rho, vector<ZZ_p> &k)
{
    
    if (rho.size() != k.size())
        throw std::invalid_argument("rho and k size don't match");
    vector<ZZ_p> pi=Frand(k.size()/2);
    
    shuffle(pi, rho);
    shuffle(pi, k);
    vector<ZZ_p> reconstructed = reveal(rho);
    apply_perm_local(k, reconstructed);
    // print_vector(k);
}

void mpc::compose(vector<ZZ_p> &sigma, vector<ZZ_p> &rho)
{
    vector<ZZ_p> pi=Frand(sigma.size()/2);
    shuffle(pi,sigma);
    vector<ZZ_p> reconstructed = reveal(sigma);
    // print_vector(reconstructed);
    vector<ZZ_p> sigma_inv = inversePerm(reconstructed);
    apply_perm_local(rho, sigma_inv);
    unshuffle(pi, rho);
    // return rho;
}
void writeVectorToCSV(const std::vector<NTL::ZZ_p>& data, int pid)
{
    std::string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/output_" + std::to_string(pid) + ".csv";
    std::ofstream file(filename);
    if (file.is_open())
    {
        for (const NTL::ZZ_p& value : data)
        {
            file << value << ",";
        }
        file.close();
        std::cout << "Vector successfully written to CSV file." << std::endl;
    }
    else
    {
        std::cout << "Error opening the file." << std::endl;
    }
}
void mpc::genperm(vector<ZZ_p> k)
{
    vector<ZZ_p> k_i(2*this->n);
    k_i.assign(this->shares.begin(), this->shares.begin() + 2*this->n);
    // vector<ZZ_p> reconstructed_sigma = reveal(k_i);
    // cout <<this->lk <<endl;
    vector<ZZ_p> sigma = genbitperm(k_i);
    vector<ZZ_p> rho(sigma.size());
    for (int i=1; i<this->lk; i++)
    {
        k_i.assign(this->shares.begin() + 2*this->n * i, this->shares.begin() + 2*this->n * (i + 1));
        // vector<ZZ_p> reconstructed_sigma = reveal(k_i);
        // print_vector(reconstructed_sigma);
        vector<ZZ_p> prev_sigma(sigma);
        apply_shared_perm(sigma, k_i);
        rho = genbitperm(k_i);
        
        
        compose(prev_sigma, rho);
        swap(sigma,rho);
        // sigma=rho;    
    }
    vector<ZZ_p> reconstructed = reveal(sigma);
    writeVectorToCSV(reconstructed, this->pid);
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
