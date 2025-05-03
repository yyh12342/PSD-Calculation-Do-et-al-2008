#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Geometry.h"

#include <string>
#include <unordered_map>
#include <vector>

struct LJParam
{
    double sigma;
    double epsilon;
};

class LJPotential
{
public:
    static LJParam probeParam; // 아르곤 파라미터
    static std::unordered_map<std::string, LJParam> baseParams; // 다른 원소 파라미터
    static void init(const std::vector<Atom> &atoms);
    static double computePotential(const Vec3 &point, const std::vector<Atom> &atoms);
    static double getCombinedSigma(const std::string &elem);
private:
    static std::unordered_map<std::string, LJParam> combinedParams; // 결합 파라미터 (인데 Example 1&2에서는 안 씀)
};

#endif