#ifndef ACCESSIBLE_VOLUME_H
#define ACCESSIBLE_VOLUME_H

#include <vector>
#include <string>
#include <utility>
#include <unordered_map>

// 원자 좌표
struct Atom
{
    int type;
    double x;
    double y;
    double z;
    double sigmaMix;
    double sigmaMixSquare;
};

// 시뮬레이션 박스 정보
struct Box
{
    double xlo; // low
    double xhi; // high
    double ylo; // low
    double yhi; // high
    double zlo; // low
    double zhi; // high
    bool px; // x축의 주기적 경계 flag
    bool py; // y축의 주기적 경계 flag
    bool pz; // z축의 주기적 경계 flag
};



static inline double PotentialAtPoint(
    double px, double py, double pz,
    const std::vector<Atom> &atoms,
    const std::unordered_map<int,double> &sigmaMap);

// 입력 파일 파싱
void Parse(
    const std::string &filename,
    Box &simulationBox,
    std::vector<Atom> &atoms,
    int modeNum,
    const std::unordered_map<int,double> &sigmaMap);

// 몬테카를로로 지름 계산
std::vector<double> ComputeDiameters(
    const Box &simulationBox,
    const std::vector<Atom> &atoms,
    int modeNum,
    const std::unordered_map<int,double> &sigmaMap);

// histogram 구해서 (diameter, V_acc,j) 쌍 리스트 생성
std::vector<std::pair<double,double>> ComputeAccessibleVolume(
    const std::vector<double> &diameters,
    const Box &simulationBox,
    int modeNum);

// CSV 출력
void WriteCSV(
    const std::vector<std::pair<double,double>> &distribution,
    const std::string &outputFilename);

#endif