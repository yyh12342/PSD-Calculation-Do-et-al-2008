#include "Geometry.h"
// #include "PDBParser.h"
#include "Potential.h"
#include "SphereFinder.h"
#include "Histogram.h"

#include <iostream>
#include <random>
#include <cstdlib>

// Geometry
void generateDoubleSlitAtoms(double boxX, double boxY, double boxZ,
    double slitW1, double slitW2,
    double neighborDist,
    std::vector<Atom> &atoms)
{
    atoms.clear();

    // 육각형 격자
    const double nd = neighborDist; // 인접 원자 거리
    const double dy = nd * std::sqrt(3.0) / 2.0; // 육각형 배치했을 때 바로 위(y) 원자와의 y 좌표 차이

    // 두 슬릿은 각각 z = 0 ~ slitW1, z = slitW1 ~ slitW1+slitW2
    // TODO 이거 wallZ 왜 2개지?
    std::vector<double> wallZ = { 0.0, slitW1, slitW1 + slitW2, boxZ };
    wallZ = { 0.0, slitW1, slitW1 + slitW2 };

    for (double z : wallZ)
    {
        // xy 평면에 육각형 격자
        int nx = static_cast<int>(boxX / nd) + 2; // x 개수
        int ny = static_cast<int>(boxY / dy) + 2; // y 개수

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                double x = i * nd + (j % 2 ? nd/2 : 0);
                double y = j * dy;

                if (x >= 0 && x <= boxX && y >= 0 && y <= boxY)
                {
                    atoms.push_back({ "C", {x, y, z} });
                }
            }
        }
    }
}

// Geometry
void generateTubeAtoms(double radius, double height,
    double neighborDist,
    std::vector<Atom> &atoms)
{
    atoms.clear();

    // TODO 탄소나노튜브 육각형 격자 다시 확인
    int nCirc = static_cast<int>(2 * M_PI * radius / neighborDist); // 일단 이런 식으로 int 형으로 나오도록 하자
    if (nCirc < 6) nCirc = 6; // 원주에 배치되는 원자를 6개 이상으로는 만들자
    double dz = neighborDist * std::sqrt(3.0) / 2.0; // 육각형 배치했을 때 바로 위(z) 원자와의 z 좌표 차이
    int nz = static_cast<int>(height / dz) + 2; // z 개수

    for (int k = 0; k < nCirc; ++k)
    {
        double rad = 2*M_PI * k / nCirc;
        double cx = radius * std::cos(rad); // 원주에서 x 좌표
        double cy = radius * std::sin(rad); // 원주에서 y 좌표

        for (int m = 0; m < nz; ++m)
        {
            double cz = m * dz; // 원통에서 z 좌표

            if (cz >= 0 && cz <= height)
            {
                atoms.push_back({ "C", {cx, cy, cz} });
            }
        }
    }
}



int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage:\n"
                  << "  Example1 (double slit):\n"
                  << "    psdcalc 1 probeSigma probeEpsilon "
                     "boxX boxY boxZ slitW1 slitW2 carbonSigma carbonEpsilon neighborDist [numSamples]\n"
                  << "  Example2 (tube):\n"
                  << "    psdcalc 2 probeSigma probeEpsilon "
                     "tubeRadius tubeHeight carbonSigma carbonEpsilon neighborDist [numSamples]\n";
        
        return true;
    }

    std::string pdbFile = argv[1]; // argv[1]에 pdb 파일 경로 설정
    // 논문에 따라 기본 샘플 수를 100000으로 정함
    int numSamples = 100000;

    if (argc >= 3)
    {
        numSamples = std::atoi(argv[2]); // 아스키에서 int로, argv[2]는 그냥 샘플 개수 넘기는 거

        if (numSamples <= 0)
        {
            numSamples = 100000;
        }
    }

    std::string outputFile;

    if (argc >= 4)
    {
        outputFile = argv[3]; // argv[3]에 csv 파일 경로 설정
    }

    // PDB 파싱
    std::vector<Atom> atoms;
    Box box;
    PDBParser parser;

    if (!parser.parseFile(pdbFile, atoms, box))
    {
        return 1;
    }

    double Lx = box.x_max - box.x_min;
    double Ly = box.y_max - box.y_min;
    double Lz = box.z_max - box.z_min;
    double boxVolume = Lx * Ly * Lz;

    if (boxVolume <= 0)
    {
        std::cerr << "Error: Simulation box volume is zero or negative.\n";

        return 1;
    }

    // LJ 파라미터 초기화
    LJPotential::init(atoms);

    // 랜덤 생성 설정
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distX(box.x_min, box.x_max);
    std::uniform_real_distribution<double> distY(box.y_min, box.y_max);
    std::uniform_real_distribution<double> distZ(box.z_min, box.z_max);

    // 히스토그램
    // 임의로 0.01nm 간격으로 최대 2nm 설정해놓음
    Histogram histogram(0.01, 2.0);
    long long successCount = 0;

    // 몬테카를로
    for (int i = 0; i < numSamples; ++i)
    {
        Vec3 point = { distX(gen), distY(gen), distZ(gen) };
        double U = LJPotential::computePotential(point, atoms);

        if (U <= 0.0)
        {
            successCount++;
            double diameter = SphereFinder::findLargestSphereDiameter(point, atoms);
            double dV = boxVolume / numSamples;
            histogram.addVolume(diameter, dV);
        }
    }

    // csv 또는 콘솔창 출력
    double fraction = (double) successCount / numSamples;
    double accessibleVol = fraction * boxVolume;
    std::cout << "Total accessible volume = " << accessibleVol 
              << " nm^3 (" << (fraction * 100.0) << "% of box volume)\n";

    if (!outputFile.empty())
    {
        histogram.saveCSV(outputFile);
        std::cout << "PSD data saved to " << outputFile << std::endl;
    }
    else
    {
        histogram.print();
    }

    return 0;
}