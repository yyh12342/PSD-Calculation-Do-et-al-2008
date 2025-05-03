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
    static void init(const std::vector<Atom> &atoms);
    static double computePotential(const Vec3 &point, const std::vector<Atom> &atoms);
    static double getCombinedSigma(const std::string &elem);
private:
    static LJParam probeParam;
    static std::unordered_map<std::string, LJParam> baseParams;
    static std::unordered_map<std::string, LJParam> combinedParams;
};

#endif