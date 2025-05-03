#include "PDBParser.h"
#include "Potential.h"
#include "SphereFinder.h"
#include "Histogram.h"

#include <iostream>
#include <random>
#include <cstdlib>

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] // arguments vector, argv[0]는 실행 경로
                  << " <structure.pdb> [num_samples] [output.csv]\n";

        return 1;
    }

    std::string pdbFile = argv[1]; // argv[1]에 pdb 파일 경로 설정
    // 임의로 기본 샘플 수를 100000으로 정함
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