#ifndef SPHEREFINDER_H
#define SPHEREFINDER_H

#include "Geometry.h"
#include "Potential.h"

#include <vector>

class SphereFinder
{
public:
    static double findLargestSphereDiameter(const Vec3 &point, const std::vector<Atom> &atoms);
};

#endif