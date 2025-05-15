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

    Box simulationBox;
    std::vector<Atom> atoms;

    try
    {
        Parse(inFile, simulationBox, atoms);
    }
    catch (const std::exception &e)
    {
        std::cerr << "파싱 오류: " << e.what() << "\n";

        return 1;
    }

    // 박스 부피 Vbox [nm^3]
    double Vbox = std::fabs(simulationBox.xhi - simulationBox.xlo)
        * std::fabs(simulationBox.yhi - simulationBox.ylo)
        * std::fabs(simulationBox.zhi - simulationBox.zlo);

    // 지름
    auto diameters = ComputeDiameters(simulationBox, atoms);

    // Vacc
    auto distribution = ComputeAccessibleVolume(diameters, Vbox);

    // // CSV 출력
    // try
    // {
    //     WriteCSV(distribution, outFile);
    // }
    // catch (const std::exception &e)
    // {
    //     std::cerr << "출력 오류: " << e.what() << "\n";

    //     return 1;
    // }

    // std::cout << "완료. 결과를 '" << outFile << "'에 저장했습니다.\n";
    
    return 0;
}










// #define _USE_MATH_DEFINES

// #include "Geometry.h"
// #include "Potential.h"
// #include "SphereFinder.h"
// #include "Histogram.h"

// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <random>
// #include <cstdlib>
// #include <filesystem>

// // Example 1 (이중 슬릿)
// void generateDoubleSlitAtoms(double boxX, double boxY, double boxZ,
//     double slitW1, double slitW2,
//     double neighborDist,
//     std::vector<Atom> &atoms)
// {
//     atoms.clear();

//     // 육각형 격자
//     const double nd = neighborDist; // 인접 원자 거리
//     const double dy = nd * std::sqrt(3.0) / 2.0; // 육각형 배치했을 때 바로 위(y) 원자와의 y 좌표 차이

//     // 두 슬릿은 각각 z = 0 ~ slitW1, z = slitW1 ~ slitW1+slitW2
//     std::vector<double> wallZ = { 0.0, slitW1, slitW1 + slitW2 };

//     for (double z : wallZ)
//     {
//         // xy 평면에 육각형 격자
//         int nx = static_cast<int>(boxX / nd) + 2; // x 개수
//         int ny = static_cast<int>(boxY / dy) + 2; // y 개수

//         for (int j = 0; j < ny; ++j)
//         {
//             for (int i = 0; i < nx; ++i)
//             {
//                 double x = i * nd + (j % 2 ? nd/2 : 0);
//                 double y = j * dy;

//                 if (x >= 0 && x <= boxX && y >= 0 && y <= boxY)
//                 {
//                     atoms.push_back({ "C", {x, y, z} });
//                 }
//             }
//         }
//     }
// }

// // Example 2 (탄소나노튜브)
// void generateTubeAtoms(double radius, double height,
//     double neighborDist,
//     std::vector<Atom> &atoms)
// {
//     atoms.clear();

//     int nCirc = static_cast<int>(2 * M_PI * radius / neighborDist); // 일단 이런 식으로 int 형으로 나오도록 하자
//     if (nCirc < 6) nCirc = 6; // 원주에 배치되는 원자를 6개 이상으로는 만들자
//     double dz = neighborDist * std::sqrt(3.0) / 2.0; // 육각형 배치했을 때 바로 위(z) 원자와의 z 좌표 차이
//     int nz = static_cast<int>(height / dz) + 2; // z 개수

//     for (int k = 0; k < nCirc; ++k)
//     {
//         double phi = 2*M_PI * k / nCirc;
//         double cx = radius * std::cos(phi); // 원주에서 x 좌표
//         double cy = radius * std::sin(phi); // 원주에서 y 좌표

//         for (int m = 0; m < nz; ++m)
//         {
//             double cz = m * dz; // 원통에서 z 좌표

//             if (cz >= 0 && cz <= height)
//             {
//                 atoms.push_back({ "C", {cx, cy, cz} });
//             }
//         }
//     }
// }



// // cfg 파일에 파라미터들 입력하면 쉽게 사용 가능
// std::vector<std::string> loadConfig(const std::string &path)
// {
//     std::ifstream ifs(path);

//     if (!ifs)
//     {
//         std::cerr << "Error: Could not open config file " << path << "\n";
//         std::exit(1);
//     }

//     std::vector<std::string> tokens;
//     std::string line, tok;

//     while (std::getline(ifs, line))
//     {
//         std::istringstream iss(line);

//         while (iss >> tok)
//         {
//             // 주석 무시
//             if (tok[0] == '#')
//             {
//                 break;
//             }

//             tokens.push_back(tok);
//         }
//     }

//     return tokens;
// }



// int main(int argc, char** argv)
// {
//     int exampleType;
//     std::cout << "Select example number (1: double slit, 2: nanotube): ";
//     std::cin >> exampleType;

//     double probeSigma, probeEpsilon;
//     std::cout << "probe sigma [nm]: ";    std::cin >> probeSigma;
//     std::cout << "probe epsilon [K]: ";   std::cin >> probeEpsilon;

//     double carbonSigma, carbonEpsilon, neighborDist;
//     std::cout << "carbon sigma [nm]: ";   std::cin >> carbonSigma;
//     std::cout << "carbon epsilon [K]: ";  std::cin >> carbonEpsilon;
//     std::cout << "neighbor distance [nm]: ";
//     std::cin >> neighborDist;

//     double boxX, boxY, boxZ, slitW1, slitW2;
//     double tubeR, tubeH;

//     if (exampleType == 1)
//     {
//         std::cout<<"box X Y Z [nm]: "; std::cin>>boxX>>boxY>>boxZ;
//         std::cout<<"slit widths [nm] (w1 w2): "; std::cin>>slitW1>>slitW2;
//     }
//     else
//     {
//         std::cout<<"tube radius [nm]: "; std::cin>>tubeR;
//         std::cout<<"tube height [nm]: "; std::cin>>tubeH;
//     }

//     int numSamples;
//     std::cout << "Number of Monte Carlo samples: ";
//     std::cin >> numSamples;

//     // Atom 생성
//     std::vector<Atom> atoms;

//     if (exampleType == 1)
//     {
//         generateDoubleSlitAtoms(boxX, boxY, boxZ,
//             slitW1, slitW2,
//             neighborDist,
//             atoms);
//     }
//     else
//     {
//         generateTubeAtoms(tubeR, tubeH,
//             neighborDist,
//             atoms);
//     }

//     // 시뮬레이션 박스 계산
//     Box box;

//     if (exampleType == 1)
//     {
//         // 앞 세 개는 최소, 뒤 세 개는 최대
//         box = {0,0,0, boxX,boxY,boxZ};
//     }
//     else
//     {
//         // 앞 세 개는 최소, 뒤 세 개는 최대
//         box = {-tubeR, -tubeR, 0, tubeR, tubeR, tubeH};
//     }

//     double volBox = (box.x_max - box.x_min) * (box.y_max - box.y_min) * (box.z_max - box.z_min);

//     // LJ 파라미터 초기화
//     LJPotential::probeParam = {probeSigma, probeEpsilon};
//     LJPotential::baseParams = {{ "C", {carbonSigma, carbonEpsilon} }};
//     LJPotential::init(atoms);

//     // 몬테카를로 & 히스토그램
//     std::mt19937 gen(std::random_device{}());
//     std::uniform_real_distribution<double> dx(box.x_min, box.x_max);
//     std::uniform_real_distribution<double> dy(box.y_min, box.y_max);
//     std::uniform_real_distribution<double> dz(box.z_min, box.z_max);
//     Histogram hist(0.05, 2.0);

//     int success=0;
//     const double minDiameter = neighborDist * 0.5;

//     for (int i=0; i<numSamples; ++i)
//     {
//         Vec3 point{ dx(gen), dy(gen), dz(gen) };

//         if (LJPotential::computePotential(point, atoms) <= 0.0)
//         {
//             double d = SphereFinder::findLargestSphereDiameter(point, atoms);

//             if (d > minDiameter)
//             {
//                 ++success;

//                 hist.addVolume(d, volBox/numSamples);
//             }
//         }
//     }

//     double f = double(success)/numSamples;
//     double Vacc = f*volBox;
//     std::cout<<"Accessible volume: "<<Vacc<<" nm^3 ("<<f*100<<"%)\n";

//     std::filesystem::create_directories("results");
//     std::string csv = "results/psd_data.csv";
//     hist.saveCSV(csv);

//     std::cout<<"CSV saved to "<<csv<<"\n";
//     std::cout<<"Run the Python plotting script to generate the graph.\n";

//     // 콘솔창이 바로 닫히는 것을 방지
//     std::cout << "\nPress Enter to exit...";
//     std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//     std::cin.get();

//     return 0;
// }