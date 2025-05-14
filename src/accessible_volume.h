#ifndef ACCESSIBLE_VOLUME_H
#define ACCESSIBLE_VOLUME_H

#include <vector>
#include <string>
#include <utility>

// 원자 좌표
struct Atom
{
    double x;
    double y;
    double z;
};

// 시뮬레이션 박스 정보
struct Box
{
    double xlo; // low
    double xhi; // high
    double ylo;
    double yhi;
    double zlo;
    double zhi;
    bool px; // 각 축의 주기적 경계 flag
    bool py;
    bool pz;
};

// 입력 파일 파싱
void Parse(const std::string &filename,
    Box &simulationBox,
    std::vector<Atom> &atoms);

// 몬테카를로로 지름 계산
std::vector<double> ComputeDiameters(const Box &simulationBox,
    const std::vector<Atom> &atoms);

double ComputeDmax(const Box &box,
    const std::vector<Atom> &atoms);

// histogram 구해서 (diameter, V_acc,j) 쌍 리스트 생성
std::vector<std::pair<double,double>> ComputeAccessibleVolume(const std::vector<double> &diameters,
    double Vbox);

// CSV 출력
void WriteCSV(const std::vector<std::pair<double,double>> &distribution,
    const std::string &outputFilename);

#endif