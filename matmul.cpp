#include "mpc.h"
#include "utils.h"
#include <cmath>
bool MPCinitialize(int pid, string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    // this->pid = pid;
    // globalprng.SetSeed(commonSeed);
    if (!MPCsetupChannels(ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2))
    {
        cout << "Network channel failed.\n";
    }
    // if (!setupSeeds())
    // {
    //     cout << "PRNG setup failed.\n";
    // }
    return true;
}
bool MPCsetupChannels(string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    Channel dataowner;
    Channel toOwner;
    Channel fromPlus;
    Channel fromMinus;
    Channel toPlus;
    Channel toMinus;  
    Endpoint p_owner(ios, ownerIP, ownerPort, EpMode::Client);
    Endpoint ownersend(ios, ownerIP, toOwnerPort, EpMode::Server);
    Endpoint eprec1(ios, address1, recPort1, EpMode::Client);
    Endpoint eprec2(ios, address2, recPort2, EpMode::Client);
    Endpoint epsend1(ios, address1, sendPort1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendPort2, EpMode::Server);
    dataowner = p_owner.addChannel();
    toOwner = ownersend.addChannel();
    fromPlus = eprec1.addChannel();
    fromMinus = eprec2.addChannel();
    toPlus = epsend1.addChannel();
    toMinus = epsend2.addChannel();
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
bool isPortOpen(int port) {
    int sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) {
        perror("socket");
        return false;
    }

    sockaddr_in addr{};
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    addr.sin_addr.s_addr = htonl(INADDR_ANY);

    int result = bind(sockfd, (struct sockaddr*)&addr, sizeof(addr));
    if (result < 0) {
        close(sockfd);
        return false;
    }

    close(sockfd);
    return true;
}
int main()
{
    
    int startingPort = 12300;
    int assignedPorts=0;
    vector<int> openPorts;
    while (assignedPorts < 48)
    {
        int port = startingPort;
        if(isPortOpen(port))
        {
            openPorts.push_back(port);
            assignedPorts++;
        }
        startingPort++;
    }
    // print_vector(openPorts);
    vector<string> address{"localhost","localhost","localhost","localhost"};

    auto start = chrono::high_resolution_clock::now();
    int numCores = thread::hardware_concurrency();
    // mutex mtx;
    // condition_variable cv;
    // std::atomic<int> readyCounter(0);
    vector<thread> threads;
    vector<vector<double>> resultVec1,resultVec2,resultVec3,resultVec4;
    thread thread1([&]() {
        startMPCset(openPorts, 0, address, 0, 1, samplecount, zscorefile, 0, resultVec1);
    });
    thread thread2([&]() {
        startMPCset(openPorts, 12, address, 1, 2,samplecount,zscorefile, 0, resultVec2);
    });
    thread thread3([&]() {
        startMPCset(openPorts, 24, address, 2, 3,samplecount,zscorefile, 0, resultVec3);
    });
    thread thread4([&]() {
        startMPCset(openPorts, 36, address, 3, 4,samplecount,zscorefile, 0, resultVec4);
    });
    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));
    threads.emplace_back(move(thread3));
    threads.emplace_back(move(thread4));
    for (auto& thread : threads) {
        thread.join();
    }
    cout << resultVec1.size() << endl;
    // cout << resultVec2[0].size() << endl;
    resultVec1.insert(resultVec1.end(), resultVec2.begin(), resultVec2.end());
    cout << resultVec1.size() << endl;

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    cout << "Execution time: " << durationInminutes << " minutes" << endl;

    return 0;
}