#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>
#include <algorithm>
#include <vector>
#include <string>

// vector3
struct Vec3
{
    double x;
    double y;
    double z;

    // 벡터 덧셈
    Vec3 operator+(const Vec3 &o) const
    {
        return { x + o.x, y + o.y, z + o.z };
    }
    // 벡터 뺄셈
    Vec3 operator-(const Vec3 &o) const
    {
        return { x - o.x, y - o.y, z - o.z };
    }
    // 스칼라 곱
    Vec3 operator*(double s) const
    {
        return { x * s, y * s, z * s };
    }
};

// 원자 정보
struct Atom
{
    std::string element; // 원소 기호
    Vec3 position; // xyz 좌표 [nm]
};

// 시뮬레이션 박스
struct Box
{
    double x_min, y_min, z_min;
    double x_max, y_max, z_max;
};

// 벡터 내적
inline double dot(const Vec3 &a, const Vec3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// 벡터 외적
inline Vec3 cross(const Vec3 &a, const Vec3 &b)
{
    return
    {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

// 거리 제곱
inline double distanceSquared(const Vec3 &a, const Vec3 &b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;

    return dx*dx + dy*dy + dz*dz;
}

// 거리
inline double distance(const Vec3 &a, const Vec3 &b)
{
    return ::sqrt(distanceSquared(a, b));
}

#endif