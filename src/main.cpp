#include "accessible_volume.h"

#include <iostream>
#include <algorithm>

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <input_dump> <mode: graphite|cnt> [output_csv]\n";

        return 1;
    }

    std::string inFile = argv[1];
    std::string mode = argv[2];
    std::string outFile = (argc >= 4 ? argv[3] : inFile + ".csv");

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

    Box simulationBox;
    std::vector<Atom> atoms;

    try
    {
        Parse(inFile, simulationBox, atoms, modeNum);
    }
    catch (const std::exception &e)
    {
        std::cerr << "파싱 오류: " << e.what() << "\n";

        return 1;
    }

    // 지름
    auto diameters = ComputeDiameters(simulationBox, atoms, modeNum);

    // Vacc
    auto distribution = ComputeAccessibleVolume(diameters, simulationBox, modeNum);

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
    
    return 0;
}