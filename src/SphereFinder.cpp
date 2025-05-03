#include "Geometry.h"
#include "SphereFinder.h"

#include <cmath>
#include <limits>

// 세 점으로 구를 찾아서 지름 계산 [nm]
double SphereFinder::findLargestSphereDiameter(const Vec3 &point, const std::vector<Atom> &atoms)
{
    int idx1 = -1;
    int idx2 = -1;
    int idx3 = -1;
    double d1 = std::numeric_limits<double>::max();
    double d2 = std::numeric_limits<double>::max();
    double d3 = std::numeric_limits<double>::max();

    // 논문의 계산 방법에 따라 가장 가까운 세 점 좌표 계산
    for (int i = 0; i < atoms.size(); ++i)
    {
        double dist2 = distanceSquared(point, atoms[i].position);

        if (dist2 < d1)
        {
            d3 = d2; idx3 = idx2;
            d2 = d1; idx2 = idx1;
            d1 = dist2; idx1 = i;
        }
        else if (dist2 < d2)
        {
            d3 = d2; idx3 = idx2;
            d2 = dist2; idx2 = i;
        }
        else if (dist2 < d3)
        {
            d3 = dist2; idx3 = i;
        }
    }

    if (idx3 < 0)
    {
        return 0.0;
    }

    const Vec3 &P1 = atoms[idx1].position;
    const Vec3 &P2 = atoms[idx2].position;
    const Vec3 &P3 = atoms[idx3].position;
    
    // 세 점 가지고 구의 외심 계산
    Vec3 v12 = {P2.x - P1.x, P2.y - P1.y, P2.z - P1.z};
    Vec3 v13 = {P3.x - P1.x, P3.y - P1.y, P3.z - P1.z};
    Vec3 norm = cross(v12, v13);
    double normLen2 = dot(norm, norm);
    Vec3 D;

    if (normLen2 < 1e-12) // 세 점이 거의 일직선 상에 위치
    {
        D.x = 0.5 * (P1.x + P2.x); // 가장 가까운 두 점 중간 좌표 반환
        D.y = 0.5 * (P1.y + P2.y);
        D.z = 0.5 * (P1.z + P2.z);
    }
    else
    {
        double d12 = dot(v12, v12);
        double d13 = dot(v13, v13);
        Vec3 cx1 = cross(v13, norm);
        Vec3 cx2 = cross(norm, v12);
        Vec3 numerator = { cx1.x * d12 + cx2.x * d13,
                           cx1.y * d12 + cx2.y * d13,
                           cx1.z * d12 + cx2.z * d13 };
        D.x = P1.x + numerator.x / (2 * normLen2);
        D.y = P1.y + numerator.y / (2 * normLen2);
        D.z = P1.z + numerator.z / (2 * normLen2);
    }
    
    // 최대 반경 계산
    static const double R_FACTOR = std::pow(2.0, 1.0/6.0);
    double maxRadius = std::numeric_limits<double>::max();

    for (const Atom &atom : atoms)
    {
        double dist = distance(D, atom.position); // D와 atom 좌표 사이의 거리
        double sigma = LJPotential::getCombinedSigma(atom.element);
        double r0 = R_FACTOR * sigma;
        double allowable = dist - r0; // 거기서 LJ zero-potential 반경을 뺀 값

        if (allowable < maxRadius) // 거기서 최솟값 구하기
        {
            maxRadius = allowable;
        }

    }

    if (maxRadius < 0)
    {
        maxRadius = 0;
    }

    return 2.0 * maxRadius; // 지름으로 반환
}