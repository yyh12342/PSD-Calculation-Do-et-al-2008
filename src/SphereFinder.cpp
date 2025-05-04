#include "Geometry.h"
#include "SphereFinder.h"

#include <cmath>
#include <limits>

// 세 점으로 구를 찾아서 지름 계산 [nm]
double SphereFinder::findLargestSphereDiameter(const Vec3 &rA, const std::vector<Atom> &atoms)
{
    int idx1 = -1;
    int idx2 = -1;
    int idx3 = -1;
    double d1 = std::numeric_limits<double>::max();
    double d2 = std::numeric_limits<double>::max();
    double d3 = std::numeric_limits<double>::max();

    // 1) 논문의 계산 방법에 따라 가장 가까운 세 원자 C1, C2, C3 찾기
    for (int i = 0; i < atoms.size(); ++i)
    {
        double dist2 = distanceSquared(rA, atoms[i].position);

        if (dist2 < d1)
        {
            d3 = d2;
            idx3 = idx2;
            d2 = d1;
            idx2 = idx1;
            d1 = dist2;
            idx1 = i;
        }
        else if (dist2 < d2)
        {
            d3 = d2;
            idx3 = idx2;
            d2 = dist2;
            idx2 = i;
        }
        else if (dist2 < d3)
        {
            d3 = dist2;
            idx3 = i;
        }
    }

    Vec3 C1 = atoms[idx1].position;
    Vec3 C2 = atoms[idx2].position;
    Vec3 C3 = atoms[idx3].position;

    // 2) B 계산
    Vec3 v1 = rA - C1;
    Vec3 v2 = rA - C2;
    double c1 = dot(v1, v1);
    double c2 = dot(v2, v2);
    double denom1 = 2*dot(v1, v2);
    double lambda1 = (c2 - c1) / (denom1 != 0 ? denom1 : 1e-12);
    Vec3 rB = rA + v1 * lambda1;

    // 3) D 계산
    Vec3 v3 = (rB - C1) + (rB - C2);
    double b1 = distanceSquared(rB, C3) - distanceSquared(rB, C1);
    double denom2 = 2*dot(v3, C1 - C3);
    double lambda2 = b1 / (denom2 != 0 ? denom2 : 1e-12);
    Vec3 rD = rB + v3 * lambda2;

    // 4) 각 C 방향에 대해 퍼텐셜이 0이 되는 점 탐색
    static const double R_FACTOR = std::pow(2.0, 1.0/6.0);

    auto findRmax = [&](const Vec3 &Ci)
    {
        Vec3 dir = normalize(rD - Ci);
        double lo = 0.0, hi = 5.0;

        auto U = [&](double lam)
        {
            Vec3 p = rD + dir * lam;

            return LJPotential::computePotential(p, atoms);
        };
        
        // hi는 potential(hi) < 0 여야 하는데 그냥 안전하게 고정
        for(int iter=0; iter<40; ++iter)
        {
            double mid = 0.5*(lo + hi);

            if(U(mid) > 0)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        return hi;
    };

    double r1 = findRmax(C1);
    double r2 = findRmax(C2);
    double r3 = findRmax(C3);

    double rMax = std::max({r1, r2, r3});
    if(rMax < 0)
    {
        rMax = 0;
    }

    return 2.0 * rMax;
}