#include "mpc.h"
#include "utils.h"
// mpc::mpc(std::atomic<int>& counter) : readyCounter(counter) {}
vector<ZZ_p> convVec(vector<uint32_t> v){
    vector<ZZ_p> converted(v.size());
    for (int i=0; i<v.size(); i++)
    {
        converted[i] = to_ZZ_p(v[i]);
    }
    return converted;
}
vector<uint32_t> convVec(vector<ZZ_p> v){
    vector<uint32_t> converted(v.size());
    for (int i=0; i<v.size(); i++)
    {
        converted[i] = conv<uint32_t>(v[i]);
    }
    return converted;
}

double cal_accuracy(vector<double>& predicted, string& filename, int N, int numCol)
{
   vector<double> actual = getRowFromMatrixFile(filename,N, numCol);
    // vector<double> actual = CSVtoVector(filename);
    
   if (predicted.size() != actual.size()) {
        cout << "predicted size: " + to_string(predicted.size()) << endl;
        cout << "actual size: " + to_string(actual.size()) << endl;
        // throw runtime_error("Vector sizes do not match");
    }

    double sumSquaredDiff = 0.0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        double diff = predicted[i] - actual[i];
        sumSquaredDiff += diff * diff;
    }

    double mse = sumSquaredDiff / predicted.size();
    double rmse = sqrt(mse);
    // cout << string("RMSE: " + to_string(rmse)) << endl;
    auto max = max_element(actual.begin(), actual.end());
    auto min = min_element(actual.begin(), actual.end());
    double normalized = *max - *min;
    return 1.0 - rmse/normalized;

}
bool mpc::initialize(int pid, string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    this->pid = pid;
    globalprng.SetSeed(commonSeed);
    if (!setupChannels(ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2))
    {
        cout << "Network channel failed.\n";
    }
    if (!setupSeeds())
    {
        cout << "PRNG setup failed.\n";
    }
    this->dataowner.recv(this->p);
    ZZ_p::init(to_ZZ(this->p));
    this->p = p;
    return true;
}
bool mpc::setupChannels(string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    Endpoint p_owner(ios, ownerIP, ownerPort, EpMode::Client);
    Endpoint ownersend(ios, ownerIP, toOwnerPort, EpMode::Server);
    Endpoint eprec1(ios, address1, recPort1, EpMode::Client);
    Endpoint eprec2(ios, address2, recPort2, EpMode::Client);
    Endpoint epsend1(ios, address1, sendPort1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendPort2, EpMode::Server);
    this->dataowner = p_owner.addChannel();
    this->toOwner = ownersend.addChannel();
    this->fromPlus = eprec1.addChannel();
    this->fromMinus = eprec2.addChannel();
    this->toPlus = epsend1.addChannel();
    this->toMinus = epsend2.addChannel();
    cout << "Established channels with the computing parties and data owner.\n";
    p_owner.stop();
    ownersend.stop();
    eprec1.stop();
    eprec2.stop();
    epsend1.stop();
    epsend2.stop();
    // ios.stop();
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
        this->localprng = new PRNG(oc::toBlock(this->pid)); //NEEDS CHANGE

        uint32_t plusSeed = this->localprng->get<uint32_t>();
        this->toPlus.send(plusSeed);

        uint32_t minusSeed;
        this->fromMinus.recv(minusSeed);

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

void mpc::ready()
{
    bool empty = areVectorsEmpty();
    if (empty)
        cout << string(to_string(this->pid) + " ready to receive.\n");
    int ready = static_cast<int>(empty);
    this->toOwner.send(ready);
}
void mpc::receiveSecrets()
{
    vector<uint32_t> dest;
    this->dataowner.recv(dest);
    cout << "vector size received.\n";
    // this->p = dest[2];
    // ZZ_p::init(to_ZZ(dest[2])); ///p
    this->n = dest[0]; //number of secrets
    this->lk = dest[1]; //number of bits
    this->inv = dest[3]; //inv 
    vector<uint32_t> result;
    try 
    {
        // receiving shares;
        this->dataowner.recv(result);
        vector<ZZ_p> temp(result.size());
        temp = convVec(result);
        // vector<ZZ_p> temp(result.begin(), result.end());
        swap(this->shares, temp);
        cout <<"secrets received.\n";
        bool prelim = this->identity.empty() && this->zscores.empty();
        int ready = static_cast<int>(prelim);
        this->toOwner.send(ready);
        if(prelim)
        {   
            cout << string(to_string(this->pid) + " ready to receive identity and zscores.\n");
            //receiving identity shares
            vector<uint32_t> identity_result;
            this->dataowner.recv(identity_result);
            vector<ZZ_p> identity_temp(identity_result.begin(), identity_result.end());
            swap(this->identity, identity_temp);
            
            // receiving zscore shares
            double shiftsize;
            this->dataowner.recv(shiftsize);
            this->shiftsize=shiftsize;

            vector<uint32_t> zscore_result;
            this->dataowner.recv(zscore_result);
            vector<ZZ_p> zscore_temp(zscore_result.begin(), zscore_result.end());
            swap(this->zscores, zscore_temp);
        }
        
    }
    catch (const exception &e)
    {
        cout << "except\n";
        cout << e.what() << endl;
    }
    // cout <<this->pid <<endl;
    // print_vector(this->identity);
    cout << string(to_string(this->pid) + " Secrets received.\n");
}

// template class std::vector<ZZ_p>;

void mpc::assertSize(vector<ZZ_p> pi, string tag) {
    for (ZZ_p z : pi) {
        uint32_t value = conv<uint32_t>(z);
        if (value > pi.size() || value < 1) {
            string output = string("[" + to_string(this->pid) + "][" + tag + "] Error! value: " + to_string(value) + " size: " + to_string(pi.size()) + "\n");
            cout << output;
            throw logic_error("["+tag+"]: permutation is bigger than vector size.");
        }
    }    
}
vector<uint32_t> mpc::apply_plaintext_perm(vector<uint32_t> rho, vector<uint32_t> sigma)
{
    //Applying rho on sigma rho o sigma
    if (rho.size() != sigma.size())
        throw invalid_argument("This is regular apply permutation; both pis should be same size.");
    vector<uint32_t> appliedvec(rho.size());
    for (int i=0; i < sigma.size(); i++)
    {
        appliedvec[rho[i]-1] = sigma[i];
    }
    return appliedvec;
}

vector<ZZ_p> mpc::Frand(uint32_t bufferSize) 
{
    if (bufferSize <= 1)
    {
        throw std::invalid_argument("bufferSize must be positive");
    }
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &plus_PRNG = *(plus_it->second);
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    vector<uint32_t> pi(2*bufferSize);
    for (uint32_t i = 0; i < bufferSize; i++)
    {
        pi[i*2+0] = i + 1;
        pi[i*2+1] = i + 1;
    }
    for (uint32_t i = bufferSize - 1; i >= 1; i--)
    {
        // random index chosen to be swapped
        uint32_t j = minus_PRNG.get<uint32_t>() % (i + 1);
        uint32_t temp = pi[i*2+0];
        pi[i*2+0] = pi[j*2+0];
        pi[j*2+0] = temp;

        uint32_t k = plus_PRNG.get<uint32_t>() % (i + 1);
        uint32_t temp2 = pi[i*2+1];
        pi[i*2+1] = pi[k*2+1];
        pi[k*2+1] = temp2;
    } 
    vector<ZZ_p> ZZ_pi(pi.begin(), pi.end());
    assertSize(ZZ_pi, "Frand");
    return ZZ_pi;
}

vector<ZZ_p> mpc::reveal(vector<ZZ_p>& pi, bool isperm)
{
    // cout <<"reveal entered.\n";
    vector<uint32_t> share1;
    vector<uint32_t> share2;
    for (int i=0; i<pi.size(); i+=2)
    {
        uint32_t r = conv<uint32_t>(pi[i]);
        uint32_t q = conv<uint32_t>(pi[i+1]);
        if(r >= this->p) {
            cout << "r too big: " << r << endl;
        }
        if(q >= this->p) {
            cout << "q too big: " << q << endl;
        }
        share1.push_back(r);
        share2.push_back(q);
    }
    this->toPlus.send(share1);
    this->toMinus.send(share2);
    vector<uint32_t> receivedMinus, receivedPlus;
    this->fromMinus.recv(receivedMinus);
    this->fromPlus.recv(receivedPlus);
    
    
    if (receivedMinus != receivedPlus) 
    {
        cout << "error! in " << this->pid << endl;
        print_vector(receivedMinus);
        print_vector(receivedPlus);
        try 
        {
            throw logic_error("received inputs do not match."); // <-- THIS IS SOMETIMES THROWN BUT NOT CAUGHT
        } 
        catch (const std::string& ex) {
            std::cerr << "Exception occurred: " << ex << std::endl;
        }
        close();
    }

    vector<ZZ_p> reconstructed;
    vector<uint32_t> final(share1.size());
    vector<uint32_t> temp(share1.size());
    if (isperm)
    {
        if (this->pid==0)
        {
            temp = apply_plaintext_perm(share2, share1);
            final = apply_plaintext_perm(receivedMinus, temp);
        }
        else if (this->pid==1)
        {
            temp=apply_plaintext_perm(share1, receivedMinus);
            final=apply_plaintext_perm(share2, temp);
        }
        else
        {
            temp=apply_plaintext_perm(receivedMinus,share2);
            final=apply_plaintext_perm(share1, temp);
        }
        reconstructed = convVec(final);
    }

    else 
    {
        for (int i=0; i<receivedMinus.size(); i++)
        {
            if (receivedMinus[i] >= this->p) {
                cout << "too big! Minus: " << receivedMinus[i] << endl;
            }
            reconstructed.push_back(conv<ZZ_p>(receivedMinus[i])+conv<ZZ_p>(share1[i])+conv<ZZ_p>(share2[i]));
            // reconstructed.push_back(conv<ZZ_p>(receivedMinus[i]));
            uint32_t result = conv<uint32_t>(reconstructed[i]);
            if (result >= this->p) {
                cout << "too big! Result: " << result << endl;
            }
        }
    }
    return reconstructed;
    
}

void mpc::reshare(vector<ZZ_p>& shares, int reshareID)
{
    vector<ZZ_p> randoms(3);
    vector<uint32_t> sendvec;
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
            sendvec.push_back(conv<uint32_t>(shares[2*j+1]));
            // this->toPlus.send(conv<uint32_t>(shares[2*j+1]));
        }
        this->toPlus.send(sendvec);
    }
    else if (reshareID == (this->pid + 2) % 3) //share to minus
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
            sendvec.push_back(conv<uint32_t>(shares[2*j]));
        }
        this->toMinus.send(sendvec);
    }
    else 
    {
       vector<ZZ_p> newshare(shares.size());
       vector<uint32_t> minus(shares.size()/2), plus(shares.size()/2);
       this->fromMinus.recv(minus);
       this->fromPlus.recv(plus);
        for (int i=0; i<shares.size()/2; i++)
        {
            newshare[2*i] = to_ZZ_p(minus[i]);
            newshare[2*i+1] = to_ZZ_p(plus[i]);
            // this->fromMinus.recv(newshare[2*i]);
            // this->fromPlus.recv(newshare[2*i+1]);
        }
        // vector<ZZ_p> temp(newshare.begin(), newshare.end());
        shares.swap(newshare); 

    }
}
void mpc::reshareM(vector<vector<ZZ_p>>& shares, int reshareID)
{
    vector<ZZ_p> randoms(3);
    vector<vector<uint32_t>> sendvec(shares.size(), vector<uint32_t>(shares[0].size()));
    if (reshareID == (this->pid + 1) % 3) // share to plus
    {
        auto minus_it = this->seedpair.find((this->pid + 2) % 3);
        PRNG &minus_PRNG = *(minus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(minus_PRNG.get<uint32_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int i=0; i<shares.size(); i++)
        {
            for (int j=0; j<shares[0].size()/2; j++)
            {
                shares[i][2*j]+=randoms[this->pid];
                shares[i][2*j+1]+=randoms[(this->pid +1)%3];
                sendvec[i][j] = conv<uint32_t>(shares[i][2*j+1]);
            }
            this->toPlus.send(sendvec[i]); 
        }
    }
    else if (reshareID == (this->pid + 2) % 3) //share to minus
    {
        auto plus_it = this->seedpair.find((this->pid + 1) % 3);
        PRNG &plus_PRNG = *(plus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(plus_PRNG.get<uint32_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int i=0; i<shares.size(); i++)
        {
            for (int j=0; j<shares[0].size()/2; j++)
            {
                shares[i][2*j]+=randoms[this->pid];
                shares[i][2*j+1]+=randoms[(this->pid +1)%3];
                sendvec[i][j] = conv<uint32_t>(shares[i][2*j]);
            }
            this->toMinus.send(sendvec[i]);
        }
    }
    else 
    {
       vector<vector<ZZ_p>> newshare(shares.size(), vector<ZZ_p>(shares[0].size()));
        for(int i=0; i<shares.size(); i++)
        {
            vector<uint32_t> minus(shares[0].size());
            vector<uint32_t> plus(shares[0].size());
            this->fromMinus.recv(minus);
            this->fromPlus.recv(plus);
            for(int j=0; j<shares[0].size()/2; j++)
            {
                newshare[i][2*j] = to_ZZ_p(minus[j]);
                newshare[i][2*j+1] = to_ZZ_p(plus[j]);
            }
        }
        // vector<vector<ZZ_p>> temp(newshare.begin(), newshare.end());
        shares.swap(newshare); 
    }
}
//is adding alpha beta gamma necessary? 
vector<ZZ_p> mpc::Fmult(vector<ZZ_p> k_i, vector<ZZ_p> s_i)
{
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    PRNG &plus_PRNG = *(plus_it->second);
    vector<ZZ_p> t_i(2);
    ZZ_p ri = k_i[0]*s_i[0] + k_i[1]*s_i[0] + k_i[0]*s_i[1];
    ri+=to_ZZ_p(minus_PRNG.get<uint32_t>()); // result + r1 - r2;
    ri-=to_ZZ_p(plus_PRNG.get<uint32_t>());
    t_i[0]=ri;
    this->toMinus.send(conv<uint32_t>(ri));
    uint32_t received_data;
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

vector<ZZ_p> mpc::inversePerm(vector<ZZ_p> pi)
{
    assertSize(pi, "inversePerm");
    vector<uint32_t> arr2(pi.size());
    vector<ZZ_p> inverse(pi.size());
    for (uint32_t i = 0; i < pi.size(); i++)
        arr2[conv<uint32_t>(pi[i]) - 1] = i + 1;
 
    for (uint32_t i = 0; i < pi.size(); i++)
    {
        inverse[i] = conv<ZZ_p>(arr2[i]);
    }
    return inverse;
}

void mpc::apply_perm_local(bool participate, vector<ZZ_p> &v, vector<ZZ_p> &pi)
{
    if (participate)
    {
        for (ZZ_p z : pi) 
        {
            uint32_t value = conv<uint32_t>(z);
            if (value > pi.size()) 
            {
                cout << "apply_perm_local." << endl;
                cout << "value: " << value << endl;
                cout << "size of vector: " << pi.size() << endl;
                throw logic_error("permutation is bigger than vector size.");
            }
        }
        if (pi.size() != v.size()/2)
            throw std::invalid_argument("local permutation should be half the size of shared vector");
        if (v.size() % 2 != 0)
            throw invalid_argument("v size is not divisible by 2.");
        vector<ZZ_p> v2(v.size());
        for (uint32_t i = 0; i < v.size()/2; i++)
        {
            if (2 * i + 1 >= v.size()) 
            {
                cout << "i too high" << i << "size: " << v.size() << endl;
            }
            try
            {
                uint32_t idx1 = conv<uint32_t>((pi[i])-1)*2;
                uint32_t idx2 = conv<uint32_t>((pi[i])-1)*2+1;
                // v2[(pi_i-1)*2] =v[2*i];
                // v2[(pi_i-1)*2+1]=v[2*i+1];
                v2[idx1] =v[2*i];
                v2[idx2]=v[2*i+1];
            }
            catch (const std::exception &e)
            {
                cerr << "Exception caught: " << e.what() << std::endl;
            }
        }
        // return v2;
        v = std::move(v2);
    }
}
void mpc::apply_perm_localM(bool participate, vector<vector<ZZ_p>> &v, vector<ZZ_p> &pi)
{
    if (participate)
    {
        for (ZZ_p z : pi) 
        {
            uint32_t value = conv<uint32_t>(z);
            if (value > pi.size()) 
            {
                cout << "apply_perm_local." << endl;
                cout << "value: " << value << endl;
                cout << "size of vector: " << pi.size() << endl;
                throw logic_error("permutation is bigger than vector size.");
            }
        }
        if (pi.size() != v.size())
            throw std::invalid_argument("local permutation should be half the size of shared vector");
        // if (v.size() % 2 != 0)
        //     throw invalid_argument("v size is not divisible by 2.");
        vector<vector<ZZ_p>> v2(v.size(), vector<ZZ_p>(v[0].size()));
        for (size_t i=0; i<v.size(); i++)
        {
            uint32_t idx = conv<uint32_t>(pi[i]-1);
            v2[idx] = v[i];
        }
        // for (uint32_t i = 0; i < v.size()/2; i++)
        // {
        //     if (2 * i + 1 >= v.size()) 
        //     {
        //         cout << "i too high" << i << "size: " << v.size() << endl;
        //     }
        //     try
        //     {
        //         uint32_t idx1 = conv<uint32_t>((pi[i])-1)*2;
        //         uint32_t idx2 = conv<uint32_t>((pi[i])-1)*2+1;
        //         // v2[(pi_i-1)*2] =v[2*i];
        //         // v2[(pi_i-1)*2+1]=v[2*i+1];
        //         v2[idx1] =v[2*i];
        //         v2[idx2]=v[2*i+1];
        //     }
        //     catch (const std::exception &e)
        //     {
        //         cerr << "Exception caught: " << e.what() << std::endl;
        //     }
        // }
        // return v2;
        v = std::move(v2);
    }
}
void mpc::shuffle(vector<ZZ_p> &pi, vector<ZZ_p> &a)
{
    assertSize(pi, "Shuffle");
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
        apply_perm_local(true, a, pi_m);
        reshare(a, 1);
        apply_perm_local(true, a, pi_p);
        reshare(a, 2);
        apply_perm_local(false, a, pi);
        reshare(a, 0);
    }
    else if (this->pid == 1)
    {
        apply_perm_local(false, a, pi);
        reshare(a, 1);
        apply_perm_local(true,a, pi_m);
        reshare(a, 2);
        apply_perm_local(true,a, pi_p);
        reshare(a, 0);
    }
    else
    {
        apply_perm_local(true, a, pi_p);
        reshare(a, 1);
        apply_perm_local(false, a, pi);
        reshare(a, 2);
        apply_perm_local(true, a, pi_m);
        reshare(a, 0);
    }
}
void mpc::unshuffle(vector<ZZ_p> &pi, vector<ZZ_p> &b)
{
    assertSize(pi, "unshuffle");
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
        apply_perm_local(false,b,inv_p);
        reshare(b,0);
        apply_perm_local(true,b,inv_p);
        reshare(b,2);
        apply_perm_local(true,b,inv_m);
        reshare(b,1);
    }
    else if (this->pid == 1)
    {
        apply_perm_local(true,b,inv_p);
        reshare(b, 0);
        apply_perm_local(true,b,inv_m);
        reshare(b, 2);
        apply_perm_local(false,b,inv_p);
        reshare(b, 1);
    }
    else
    {
        apply_perm_local(true,b,inv_m);
        reshare(b, 0);
        apply_perm_local(false,b,inv_p);
        reshare(b, 2);
        apply_perm_local(true,b,inv_p);
        reshare(b, 1);
    }
}

void mpc::shuffleM(vector<ZZ_p> &pi, vector<vector<ZZ_p>> &a)
{
    assertSize(pi, "Shuffle");
    // if (pi.size() != a.size())
    //     throw std::invalid_argument("Your shares and pi size are different");
    vector<ZZ_p> pi_m(pi.size()/2), pi_p(pi.size()/2);
    for (int i=0; i<pi.size()/2; i++)
    {
        pi_m[i] = pi[2*i];
        pi_p[i] = pi[2*i+1];
    }
    if (this->pid == 0)
    {   
        apply_perm_localM(true, a, pi_m);
        reshareM(a, 1);
        apply_perm_localM(true, a, pi_p);
        reshareM(a, 2);
        apply_perm_localM(false, a, pi);
        reshareM(a, 0);
    }
    else if (this->pid == 1)
    {
        apply_perm_localM(false, a, pi);
        reshareM(a, 1);
        apply_perm_localM(true,a, pi_m);
        reshareM(a, 2);
        apply_perm_localM(true,a, pi_p);
        reshareM(a, 0);
    }
    else
    {
        apply_perm_localM(true, a, pi_p);
        reshareM(a, 1);
        apply_perm_localM(false, a, pi);
        reshareM(a, 2);
        apply_perm_localM(true, a, pi_m);
        reshareM(a, 0);
    }
}
void mpc::testMatrix(vector<ZZ_p>& pi)
{
    if(this->pid !=0)
    {
        print_vector(this->pheno);
        print_vector(pi);
    }
        
    shuffleM(pi, this->pheno);
    // vector<ZZ_p> revealp = reveal(this->pheno);
    // if(this->pid !=0)
    //     print_vector(this->pheno);
    
}
void mpc::apply_shared_perm(vector<ZZ_p> &rho, vector<ZZ_p> &k)
{
    if (rho.size() != k.size()){
        cout << string("rho: "+ to_string(rho.size()) + ", k: "+ to_string(k.size()));
        throw std::invalid_argument("rho and k size don't match");
    }
        
    vector<ZZ_p> pi=Frand(k.size()/2);
    
    shuffle(pi, rho);
    shuffle(pi, k);
    vector<ZZ_p> reconstructed = reveal(rho, false); 
    assertSize(reconstructed, "apply_shared_perm");
    if (reconstructed.size() == k.size()/2)
    {
        vector<ZZ_p> v2(k.size());
        for (size_t i = 0; i < k.size()/2; i++)
        {
            uint64_t pi_i = conv<uint64_t>(reconstructed[i]);
            v2[(pi_i-1)*2] =k[2*i];
            v2[(pi_i-1)*2+1]=k[2*i+1];
        }
        k = std::move(v2);

    }
    else 
    {
        cout << "reconstructed size: " << reconstructed.size() << endl;
        cout << "vector size:" << k.size() << endl;
        throw logic_error("Apply shared perm: Asking to apply different size permutations");
    }
    // apply_perm_local(k, reconstructed);
    // print_vector(k);
}

void mpc::compose(vector<ZZ_p> &sigma, vector<ZZ_p> &rho)
{
    vector<ZZ_p> pi=Frand(sigma.size()/2);
    shuffle(pi,sigma);
    vector<ZZ_p> reconstructed = reveal(sigma, false);
    assertSize(reconstructed, "compose");
    vector<ZZ_p> sigma_inv = inversePerm(reconstructed);

    if (sigma_inv.size() == rho.size()/2)
    {
        vector<ZZ_p> v2(rho.size());
        for (size_t i = 0; i < rho.size()/2; i++)
        {
            uint64_t pi_i = conv<uint64_t>(sigma_inv[i]);
            v2[(pi_i-1)*2] =rho[2*i];
            v2[(pi_i-1)*2+1]=rho[2*i+1];
        }
        rho = std::move(v2);
    }
    else 
    {
        throw logic_error("Compose: Asking to apply different size permutations");
    }
    unshuffle(pi, rho);
    // return rho;
}
template <typename T>
void writeVectorToCSV(const std::vector<T>& data, int pid, int row, string name)
{
    std::string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name + "_pid" + std::to_string(pid) + "_row" + to_string(row)+ ".csv";
    std::ofstream file(filename);
    if (file.is_open())
    {
        for (const T& value : data)
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

vector<double> mpc::genperm(int row, int numCol)
{
    vector<ZZ_p> k_i(2*this->n);
    k_i.assign(this->shares.begin(), this->shares.begin() + 2*this->n);
    cout <<"assigned" <<endl;
    vector<ZZ_p> inv_rho, inv_cdf;
    inv_rho.assign(this->identity.begin(), this->identity.end());
    inv_cdf.assign(this->zscores.begin(), this->zscores.end());
    try
    {
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
        cout << "loop finished.\n";
        vector<ZZ_p> reconstructed = reveal(sigma, false); 
        apply_shared_perm(sigma, inv_rho);
        apply_shared_perm(inv_rho, inv_cdf);
        vector<ZZ_p> zscores_result = reveal(inv_cdf, false); // check if this changes after running multiple times TODO

        // writeVectorToCSV(zscores_result, this->pid, 0, "unfiltered");
        vector<uint32_t> zscores_uint = convVec(zscores_result);
        vector<double> zscores_unscaled = UnscaleVector_signed(zscores_uint, pow(10,2));
        // writeVectorToCSV(zscores_unscaled, this->pid, row, "Unscaled");
        vector<double> zscores_double = UnshiftVector(zscores_unscaled,this->shiftsize);
        // print_vector(zscores_int);
        // writeVectorToCSV(zscores_int, this->pid, row, "unscaled_unshifted");
        // vector<double> zscores = UnscaleVector_signed(zscores_int, pow(10,5));
        // writeVectorToCSV(zscores, this->pid, row, "unshift_descaled");
        string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/qn.tsv";
        // string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/toyQN.tsv";
        double accuracy = cal_accuracy(zscores_double, zscore_filename, row, numCol);
        if (this->pid == 0)
        {
            cout << string("row: " + to_string(row) + "/ pid" + to_string(this->pid) + ": " + to_string(accuracy)) << endl;
            writeVectorToCSV(zscores_double, this->pid, row, "0726");
        }
            
        return zscores_double;
        // if (pid == 0) {
        //     print_vector(reconstructed);
        //     // print_vector(this->originalData);
        //     vector<int> trimData;
        //     for (int i = 0; i < this->originalData.size(); i++) {
        //         trimData.push_back(floor(this->originalData[i] * 100000));
        //     }
        //     vector<int> sortedData(trimData);
        //     sort(sortedData.begin(), sortedData.end());
        //     vector<int> result;
        //     for (int i = 0; i < trimData.size(); i++) {
        //         int j = 0;
        //         for (; j < sortedData.size() && (sortedData[j] != trimData[i] || find(result.begin(), result.end(), j + 1) != result.end()); j++);
        //         result.push_back(j + 1);
        //     }
        //     print_vector(result);
        //     vector<ZZ_p> resultZZ;
        //     for (int i = 0; i < result.size(); i++) {
        //         resultZZ.push_back(conv<ZZ_p>(result[i]));
        //     }
        //     assert_equal(resultZZ, reconstructed, "genperm result");
        // }
    }
    // catch(const std::logic_error &e){
    //     cerr << "Exception caught: " << e.what() << std::endl;
    // }
    catch (const std::exception &e)
    {
        cerr << "Exception caught: " << e.what() << std::endl;
    }
    
}
void mpc::receiveMatrix()
{
    // uint32_t p;
    
    vector<uint32_t> mat1, mat2;
    this->dataowner.recv(this->shape);
    this->dataowner.recv(mat1);
    this->dataowner.recv(mat2);
    // print_vector(this->shape);
    // cout << mat1.size() <<endl;
    // cout << mat2.size() << endl;
    for (int i=0; i<this->shape[0]; i++)
    {
        vector<ZZ_p> row;
        for (int j=0; j<this->shape[1]; j++)
        {
            // cout << string(to_string(i)+ " " + to_string(j)+"\n");
            // cout << mat1[2*i*shape[1]+2*j] << endl;
            // cout << mat1[2*i*shape[1]+2*j+1] << endl;
            row.push_back(to_ZZ_p(mat1[2*i*shape[1]+2*j]));
            row.push_back(to_ZZ_p(mat1[2*i*shape[1]+2*j+1]));
            // this->geno[i][2*j]=mat1[2*i*shape[1]+2*j];
            // this->geno[i][2*j+1]=mat1[2*i*shape[1]+2*j+1];
        }
        this->geno.push_back(row);
    }
    for (int j=0; j<this->shape[2]; j++)
    {
        vector<ZZ_p> row;
        for (int k=0; k<this->shape[3]; k++)
        {
            row.push_back(to_ZZ_p(mat2[2*j*shape[3]+2*k]));
            row.push_back(to_ZZ_p(mat2[2*j*shape[3]+2*k+1]));
            // this->pheno[j][2*k]=mat2[2*j*shape[3]+2*k];
            // this->pheno[j][2*k+1]=mat2[2*j*shape[3]+2*k+1];
        }
        this->pheno.push_back(row);
    }
    // print_vector(this->pheno);
    // vector<ZZ_p> temp(mat2.begin(), mat2.end());
    // vector<ZZ_p> recoveredp = reveal(temp,false);
    // print_vector(recoveredp);
    // vector<vector<ZZ_p>> multmat = matmult(this->geno, this->pheno, this->shape[0], this->shape[1], this->shape[2], this->shape[3]);
    // for (int i=0; i<multmat.size(); i++)
    // {
    //     vector<ZZ_p> recoveredg = reveal(multmat[i], false);
    //     if (this->pid==0)
    //         print_vector(recoveredg);
    // }
    // vector<ZZ_p> recoveredg = reveal(multmat, false);
    // print_vector(recoveredg);
}
vector<vector<ZZ_p>> transpose(vector<vector<ZZ_p>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size()/2;
    vector<vector<ZZ_p>> transposed(cols, vector<ZZ_p>(rows*2));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][2*i] = matrix[i][2*j];
            transposed[j][2*i+1] = matrix[i][2*j+1];
        }
    }
    return transposed;
}
void mpc::logRatio()
{
    // pseudo-reference is average log counts across all samples, for each gene
    uint32_t gene = this->pheno.size();
    uint32_t sample = this->pheno[0].size();
    // cout << string(to_string(this->pheno.size())+" / "+ to_string(this->pheno[0].size())+"\n");
    vector<uint32_t> sumvec;
    vector<double> pseudoref;
    vector<vector<double>> ratios(this->pheno.size(), vector<double>(this->pheno[0].size()));
    vector<double> sendback;
    for (int i=0; i<gene; i++)
    {
        ZZ_p sum=to_ZZ_p(0);
        for (int j=0; j<sample; j++)
        {
            sum+=this->pheno[i][j];
        }

        // double avg = conv<uint32_t>(sum)/static_cast<double>(sample);
        // cout <<sum<<endl;
        sumvec.push_back(conv<uint32_t>(sum));
    }
    // print_vector(sumvec);
    vector<uint32_t> v1, v2;
    this->toMinus.send(sumvec);
    this->toPlus.send(sumvec);
    this->fromMinus.recv(v1);
    this->fromPlus.recv(v2);
    for (int i=0; i<sumvec.size();i++)
    {
        ZZ_p totalsum=to_ZZ_p(sumvec[i]);
        totalsum+=to_ZZ_p(v1[i]);
        totalsum+=to_ZZ_p(v2[i]);
        // cout <<totalsum;
        pseudoref.push_back(conv<uint32_t>(totalsum)/static_cast<double>(sample));
    }
    // vector<ZZ_p> pseudosum = reveal(pseudoref,false);
    // print_vector(pseudoref);
    // for (int i=0; i<gene; i++)
    // {
    //     for(int j=0; j<sample; j++)
    //     {
    //         ratios[i][j] = conv<uint32_t>(this->pheno[i][j]) - pseudoref[i];
    //         sendback.push_back(ratios[i][j]);
    //     }
    // }
    // print_vector(sendback);
    this->toOwner.send(pseudoref);
}
vector<vector<ZZ_p>> mpc::matmult(vector<vector<ZZ_p>>& mat1, vector<vector<ZZ_p>>& mat2, int row1, int col1, int row2, int col2)
{
    if (col1 != row2)
        throw logic_error("mat1 column size and mat2 row size do not match.");
    vector<vector<ZZ_p>> result(row1, vector<ZZ_p>(col2));
    vector<vector<ZZ_p>> transposed = transpose(mat2);
    for (int i = 0; i < row1; ++i) 
    { // row index
        for (int j = 0; j < col2; ++j) 
        { // col index
            vector<ZZ_p> vec1=mat1[i];
            vector<ZZ_p> vec2=transposed[j];
            // for (int m = 0; m < row2; m++) 
            // {
                
            //     // cout << string("col2:"+to_string(col2)+" m:"+to_string(m)+"/ j:"+to_string(j)+"/ index:"+to_string(2*m*col2+j*col2)+"\n");
            //     vec2[2*m]=this->pheno[2*m*col2+j*2];
            //     vec2[2*m+1]=this->pheno[2*m*col2+j*2+1]; 
            // }
            // vector<ZZ_p> revec2 = reveal(vec2,false);
            // print_vector(vec2);
            ZZ_p sum=to_ZZ_p(0);
            for (int k=0; k<vec2.size()/2; k++)
            {
                sum += vec1[2*k]*vec2[2*k]+vec1[2*k+1]*vec2[2*k]+vec1[2*k]*vec2[2*k+1];
            }
            result[i][j] = sum;
            // cout << string("row"+to_string(i)+" col"+to_string(j)+" :"+to_string(conv<int>(sum)))<<endl;
            // vector<ZZ_p> revec1 = reveal(vec1,false);
            // print_vector(revec1);
            // for (int k = 0; k < col1; ++k) {
            //     result[i * col2 + j] += mat1[i * col1 + k] * mat2[k * col2 + j];
            // }
        }
    }
    vector<vector<ZZ_p>> finalresult(row1, vector<ZZ_p>(col2*2));
    for (int ridx=0; ridx<result.size(); ridx++)
    {
        vector<uint32_t> result1 = convVec(result[ridx]);
        this->toMinus.send(result1);
        vector<uint32_t> result2;
        this->fromPlus.recv(result2);
        if (result1.size() != result2.size())
            throw runtime_error("matmul result shares are not the same size.");
        for (int cidx=0; cidx<result[0].size(); cidx++)
        {
            finalresult[ridx][2*cidx] = to_ZZ_p(result1[cidx]);
            finalresult[ridx][2*cidx+1] = to_ZZ_p(result2[cidx]);
        }
        // vector<ZZ_p> recovered = reveal(finalresult[ridx], false);
        // print_vector(recovered);
    }
    
    // this->toMinus.send(result1);
    // vector<uint32_t> result2;
    // this->fromPlus.recv(result2);
    // if (result1.size() != result2.size())
    //     throw runtime_error("matmul result shares are not the same size.");
    // vector<ZZ_p> finalresult(result.size()*2);
    // for (size_t i=0; i<result1.size(); i++)
    // {
    //     finalresult[i*2] = conv<ZZ_p>(result1[i]);
    //     finalresult[i*2+1] = conv<ZZ_p>(result2[i]);
    // }
    // print_vector(finalresult);
    
    return finalresult;
}
void mpc::close()
{
    // clearVectors();
    cout << "Closing channels.\n";
    // std::lock_guard<std::mutex> lock(mtx);
    this->dataowner.close();
    this->fromPlus.close();
    this->fromMinus.close();
    this->toPlus.close();
    this->toMinus.close();
    // this->ios.stop();
    cout << "channels closed.\n";
}
