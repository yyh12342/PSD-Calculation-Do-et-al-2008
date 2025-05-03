#include "Geometry.h"
#include "Potential.h"

#include <cmath>
#include <iostream>

// 아르곤 파라미터
// 시그마: 입자 간 상호작용(반발력)하는 길이 [nm]
// 입실론: 퍼텐셜 E 곡선의 깊이
LJParam LJPotential::probeParam = {0.3405, 119.8};
// 기본 원소 파라미터
std::unordered_map<std::string, LJParam> LJPotential::baseParams;
// 결합 파라미터
std::unordered_map<std::string, LJParam> LJPotential::combinedParams;

// LJ 파라미터 초기화
void LJPotential::init(const std::vector<Atom> &atoms)
{
    if (baseParams.empty())
    {
        baseParams["C"] = {0.340, 28.0};
        baseParams["O"] = {0.296, 60.0};
        baseParams["N"] = {0.325, 36.0};
        baseParams["H"] = {0.250, 15.0};
        baseParams["Si"] = {0.420, 34.0};
        baseParams["Zn"] = {0.255, 62.0};
        baseParams["Li"] = {0.182, 24.0};
    }
    
    combinedParams.clear();

    // 각 원소들과 아르곤의 결합 파라미터 계산
    for (const Atom &atom : atoms)
    {
        const std::string &elem = atom.element;

        if (combinedParams.find(elem) != combinedParams.end()) // 계산했던 원소면
        {
            continue;
        }

        LJParam base = {0.340, 50.0}; // 혹시나 baseParams에 해당하는 원소가 없을 경우의 디폴트 값

        auto iter = baseParams.find(elem);
        if (iter != baseParams.end())
        {
            base = iter->second;
        }
        else
        {
            static bool warned = false;

            if (!warned)
            {
                std::cerr << "Warning: Unknown element \"" << elem 
                          << "\" using default LJ parameters." << std::endl;
                warned = true;
            }
        }
        
        // Lorentz–Berthelot 결합
        // 이종 간 상호작용에서 파라미터를 섞는 표준 방법
        LJParam comb;
        comb.sigma = 0.5 * (base.sigma + probeParam.sigma);
        comb.epsilon = std::sqrt(base.epsilon * probeParam.epsilon);
        combinedParams[elem] = comb;
    }
}

// 퍼텐셜 계산
double LJPotential::computePotential(const Vec3 &point, const std::vector<Atom> &atoms)
{
    double totalU = 0.0; // 0보다 낮은 값 찾아야 함

    for (const Atom &atom : atoms)
    {
        const LJParam &lj = combinedParams[atom.element];

        // 거리 제곱 계싼
        double dx = atom.position.x - point.x;
        double dy = atom.position.y - point.y;
        double dz = atom.position.z - point.z;
        double r2 = dx*dx + dy*dy + dz*dz;

        if (r2 == 0)
        {
            return 1e9; // 동일한 위치일 때 반발력이 큰 것을 임의로 1e9라 함
        }
        
        double inv_r2 = (lj.sigma * lj.sigma) / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;
        double u_pair = 4 * lj.epsilon * (inv_r12 - inv_r6);
        totalU += u_pair;
        
        if (totalU > 1e3) // 반발이 충분히 큰 것을 임의로 1e3이라 함
        {
            return totalU; // 그대로 조기 종료
        }
    }

    return totalU;
}

// 결합 시그마 반환
double LJPotential::getCombinedSigma(const std::string &elem)
{
    auto iter = combinedParams.find(elem);

    if (iter != combinedParams.end())
    {
        return iter->second.sigma;
    }
    
    return probeParam.sigma; // 없으면 아르곤 시그마 반환
}