#include "accessible_volume.h"

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <chrono>

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <input_dump> <mode> <output_csv> <type:sigma>...\n"
                  << " e.g. dump.demixing_equil.100000 default result_ex3_1.csv 1:0.3550 2:0.3000\n";

        return 1;
    }

    std::string inFile = argv[1];
    // std::cout << "[Debug] inFile = " << inFile << "\n";
    std::string mode = argv[2];
    // std::cout << "[Debug] mode = " << mode << "\n";
    std::string outFile = argv[3];
    // std::cout << "[Debug] outFile = " << outFile << "\n";

    int modeNum;
    if (mode == "graphite")
    {
        modeNum = 0;
    }
    else if (mode == "cnt")
    {
        modeNum = 1;
    }
    else if (mode == "default")
    {
        modeNum = 2;
    }
    else
    {
        std::cerr << "모드 입력 오류\n";

        return 1;
    }

    std::unordered_map<int,double> sigmaMap;

    for(int i = 4; i < argc; ++i)
    {
        std::string s = argv[i];
        auto pos = s.find(':');
        if (pos == std::string::npos)
        {
            continue;
        }
        int t = std::atoi(s.substr(0, pos).c_str());
        double v = std::atof(s.substr(pos+1).c_str());
        sigmaMap[t] = v;

        // std::cout << "[Debug] sigmaMap[t] = " << sigmaMap[t] << "\n";
    }

    Box simulationBox;
    std::vector<Atom> atoms;

    auto tStart = std::chrono::high_resolution_clock::now(); // 시간 측정

    try
    {
        // std::cout << "[Debug] try\n";

        Parse(inFile, simulationBox, atoms, modeNum, sigmaMap);
    }
    catch (const std::exception &e)
    {
        std::cerr << "파싱 오류: " << e.what() << "\n";

        return 1;
    }

    // 지름
    auto diameters = ComputeDiameters(simulationBox, atoms, modeNum, sigmaMap);

    // Vacc
    auto distribution = ComputeAccessibleVolume(diameters, simulationBox, modeNum);

    auto tEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> computeTime = tEnd - tStart;

    // CSV 출력
    try
    {
        WriteCSV(distribution, outFile);
    }
    catch (const std::exception &e)
    {
        std::cerr << "출력 오류: " << e.what() << "\n";

        return 1;
    }

    std::cout << "완료. 결과를 '" << outFile << "'에 저장했습니다.\n";
    std::cout << "계산 시간: " << computeTime.count() << "초\n";
    
    return 0;
}