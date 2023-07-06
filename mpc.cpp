#include "mpc.h"
#include "utils.h"
// mpc::mpc(std::atomic<int>& counter) : readyCounter(counter) {}
vector<ZZ_p> convVec(vector<uint32_t> v){
    vector<ZZ_p> converted(v.size());
    for (int i=0; i<v.size(); i++)
    {
        converted[i] = conv<ZZ_p>(v[i]);
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

double cal_accuracy(const vector<double>& predicted, const string& filename, int N, int numCol)
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

    // readyCounter = std::ref(counter);
    // ZZ_p::init(to_ZZ(4));
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
    return true;
}
bool mpc::setupChannels(string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    string to_print = "setup Channels: ownerPort: " + to_string(ownerPort);
    to_print.append(" pid " + to_string(this->pid));
    to_print.append(" recPort1 " + to_string(recPort1));
    to_print.append(" recPort2 " + to_string(recPort2));
    to_print.append("\n");
    cout << to_print;

    // IOService ios;
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
        PRNG prng(oc::toBlock(this->pid)); //initializing private prng with pid..?
        uint32_t plusSeed = prng.get<uint32_t>();
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
    this->p = dest[2];
    ZZ_p::init(to_ZZ(dest[2])); ///p
    this->n = dest[0]; //number of secrets
    this->lk = dest[1]; //number of bits
    this->inv = dest[3]; //inv 
    vector<uint32_t> result;
    try 
    {
        // this->dataowner.recv(this->originalData);
        this->dataowner.recv(result);
        // for (int i = 0; i < dest[0]*dest[1]*2; i++) 
        // {
        //     uint32_t received_data;
        //     this->dataowner.recv(received_data);
        //     result.push_back(received_data);
        // }
        vector<ZZ_p> temp(result.begin(), result.end());
        swap(this->shares, temp);
        //receiving identity shares
        vector<uint32_t> identity_result;
        this->dataowner.recv(identity_result);
        vector<ZZ_p> identity_temp(identity_result.begin(), identity_result.end());
        swap(this->identity, identity_temp);
        // for (int i=0; i<dest[0]*2; i++)
        // {
        //     uint32_t share;
        //     this->dataowner.recv(share);
        //     this->identity.push_back(conv<ZZ_p>(share));

        // }
        //receiving zscore shares
        double shiftsize;
        this->dataowner.recv(shiftsize);
        this->shiftsize=shiftsize;

        vector<uint32_t> zscore_result;
        this->dataowner.recv(zscore_result);
        vector<ZZ_p> zscore_temp(zscore_result.begin(), zscore_result.end());
        swap(this->zscores, zscore_temp);
        // for (int i=0; i<dest[0]*2; i++)
        // {
        //     uint32_t share;
        //     this->dataowner.recv(share);
        //     this->zscores.push_back(conv<ZZ_p>(share));
        // }
    }
    catch (const exception &e)
    {
        cout << "except\n";
        cout << e.what() << endl;
    }
    // cout <<this->pid <<endl;
    // print_vector(this->identity);
    cout << "Secrets received.\n";
}

template class std::vector<ZZ_p>;

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
    
    // std::unique_lock<std::mutex> lock(mtx);
    // readyCounter.get().store(readyCounter.get() + 1);
    // int counterValue = readyCounter.get().load();
    // if (counterValue % 3 == 0) {
    //     // cout << "reveal counter: " << readyCounter.get() << " pid: " << this->pid << endl;
    //     // readyCounter.get().store(0); // Reset the counter for the next iteration
    //     // cout << "readyCounter: " << readyCounter.get() << endl;
    //     lock.unlock();
    //     cv.get().notify_all();
    // } else {
    //     int div = readyCounter.get().load() / 3;
    //     // cout << "reveal counter: " << readyCounter.get() << " pid: " << this->pid << endl;
    //     cv.get().wait(lock, [this, div]{ 
    //         // cout << "notified! " << readyCounter.get().load() << " pid: " << this->pid << endl;
    //         return readyCounter.get().load() / 3 > div; 
    //         });
    // }
    // string output = string("reveal exit: ").append(to_string(this->pid)).append("!\n");
    // cout << output;
    
    
    if (receivedMinus != receivedPlus) 
    {
        cout << "error! in " << this->pid << endl;
        // print_vector(receivedMinus);
        // print_vector(receivedPlus);
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
    // mutex mutex;
    // condition_variable cv;

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
            this->toPlus.send(conv<uint32_t>(shares[2*j+1]));
        }
        // cout << this->pid << " finished.\n";

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
            this->toMinus.send(conv<uint32_t>(shares[2*j]));
        }
        // cout << this->pid << " finished.\n";

    }
    else 
    {
       vector<uint32_t> newshare(shares.size());
        for (int i=0; i<shares.size()/2; i++)
        {
            this->fromMinus.recv(newshare[2*i]);
            this->fromPlus.recv(newshare[2*i+1]);
        }
        vector<ZZ_p> temp(newshare.begin(), newshare.end());
        shares.swap(temp); 
        // cout << this->pid << " finished.\n";

    }
        // cout << this->pid << endl;
    // {
    //     std::unique_lock<std::mutex> lock(mtx);
    //     readyCounter.get().store(readyCounter.get() + 1);
    //     int counterValue = readyCounter.get().load();
    //     if (counterValue % 3 == 0) {
    //         // cout << "counter: " << readyCounter.get() << endl;
    //         // readyCounter.get().store(0); // Reset the counter for the next iteration
    //         cv.get().notify_all();
    //     } else {
    //         int div = readyCounter.get().load() / 3;
    //         // cout << "counter: " << readyCounter.get() << endl;
    //         cv.get().wait(lock, [this, div]{ return readyCounter.get().load() / 3 > div; });
    //     }
    // }
    // cout << "reshared." << endl;
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
    // cout << this->pid << " apply_perm_local entered.\n";
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
            // uint64_t pi2 = conv<uint64_t>(pi[2*i+1]);
            // v2[(pi_i-1)*2] =v[2*i];
            // v2[(pi_i-1)*2+1]=v[2*i+1];
        }
        // return v2;
        v = std::move(v2);
    }
    // std::unique_lock<std::mutex> lock(mtx);
    // readyCounter.get().store(readyCounter.get() + 1);
    // int counterValue = readyCounter.get().load();
    // if (counterValue % 3 == 0) 
    // {
    //     // cout << "counter: " << readyCounter.get() << endl;
    //     // readyCounter.get().store(0); // Reset the counter for the next iteration
    //     cv.get().notify_all();
    // } else 
    // {
    //     int div = readyCounter.get().load() / 3;
    //     // cout << "counter: " << readyCounter.get() << endl;
    //     cv.get().wait(lock, [this, div]{ return readyCounter.get().load() / 3 > div; });
    // }
}
// void mpc::apply_perm_local(vector<ZZ_p> &v, vector<ZZ_p> &pi, int partnerID)
// {
//     vector<uint32_t> share1, share2;
//     vector<uint32_t> received;
//     for (int i=0; i<v.size()/2; i++)
//     {
//         uint32_t r = conv<uint32_t>(v[2*i]);
//         uint32_t q = conv<uint32_t>(v[2*i+1]);
//         share1.push_back(r);
//         share2.push_back(q);
//     }
//     auto partner_it = this->seedpair.find(partnerID);
//     PRNG &sharedPRNG = *(partner_it->second);
//     if (partnerID==(this->pid+1)%3) //partner is plus
//     {
//         this->toPlus.send(share1);
//         this->fromPlus.recv(received);
//     }
//     if (partnerID==(this->pid+2)%3) //partner is minus
//     {
//         this->toMinus.send(share2);
//         this->fromMinus.recv(received);
//     }
    
//     vector<ZZ_p> reconstructed, newvec(received.size()), to_share(3);
//     for (int i=0; i<received.size(); i++)
//     {
//         ZZ_p result_p = conv<ZZ_p>(received[i])+conv<ZZ_p>(share1[i])+conv<ZZ_p>(share2[i]);
//         reconstructed.push_back(result_p);
//         uint32_t idx = conv<uint32_t>(pi[i]);
//         newvec[idx-1] = result_p;
//     }
//     for (int j=0; j<newvec.size(); j++)
//     {
//         to_share[0] = to_ZZ_p(sharedPRNG.get<uint32_t>());
//         to_share[1] = to_ZZ_p(sharedPRNG.get<uint32_t>());
//         to_share[2] = newvec[j] - to_share[0] - to_share[1];
//         v[2*j] = to_share[this->pid];
//         v[2*j+1] = to_share[(this->pid+1)%3];
//     }
//     // return newvec;
// }

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
    // std::unique_lock<std::mutex> lock(mtx);
    // readyCounter.get().store(readyCounter.get() + 1);
    // int counterValue = readyCounter.get().load();
    // if (counterValue % 3 == 0) 
    // {
    //     // cout << "Close counter: " << readyCounter.get() << endl;
    //     // readyCounter.get().store(0); // Reset the counter for the next iteration
    //     cv.get().notify_all();
    // } else 
    // {
    //     int div = readyCounter.get().load() / 3;
    //     // cout << "Close counter: " << readyCounter.get() << endl;
    //     cv.get().wait(lock, [this, div]{ return readyCounter.get().load() == 0  > div; });
    // }
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
        apply_shared_perm(sigma, this->identity);
        apply_shared_perm(this->identity, this->zscores);
        vector<ZZ_p> zscores_result = reveal(this->zscores, false); // check if this changes after running multiple times TODO

        // writeVectorToCSV(zscores_result, this->pid, 0, "unfiltered");
        vector<uint32_t> zscores_uint = convVec(zscores_result);
        vector<double> zscores_unscaled = UnscaleVector_signed(zscores_uint, pow(10,2));
        // writeVectorToCSV(zscores_unscaled, this->pid, row, "Unscaled");
        vector<double> zscores_double = UnshiftVector(zscores_unscaled,this->shiftsize);
        // print_vector(zscores_int);
        // writeVectorToCSV(zscores_int, this->pid, row, "unscaled_unshifted");
        // vector<double> zscores = UnscaleVector_signed(zscores_int, pow(10,5));
        // writeVectorToCSV(zscores, this->pid, row, "unshift_descaled");
        string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/accuracy_combination/correctQN.tsv";
        // string zscore_filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/toy_QN/toyQN.tsv";
        double accuracy = cal_accuracy(zscores_double, zscore_filename, row, numCol);
        if (this->pid == 0)
        {
            cout << string("row: " + to_string(row) + "/ pid" + to_string(this->pid) + ": " + to_string(accuracy)) << endl;
            // writeVectorToCSV(zscores_double, this->pid, row, "0628");
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
    catch(const std::logic_error &e){
        cerr << "Exception caught: " << e.what() << std::endl;
    }
    catch (const std::exception &e)
    {
        cerr << "Exception caught: " << e.what() << std::endl;
    }
    
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
