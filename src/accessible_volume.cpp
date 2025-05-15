#include "accessible_volume.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <random>
#include <algorithm>

static const double sigmaAr = 0.3405; // 아르곤 충돌 지름 σ [nm]
// static const double epsilonAr = 0.9970; // 아르곤 depth ε
static const double sigmaC = 0.3400; // 탄소 충돌 지름 σ [nm]
// static const double epsilonC = 0.2300; // 탄소 depth ε

// static const double cutoffFactor = 2.5; // mix σ 기준 컷오프

// static const int M = 1000000; // Monte Carlo 시도 횟수
static const int M = 10; // 디버그용
static const double Dmax = 2.0; // 최대 직경 [nm]
static const int bin = 50; // bin 개수

void Parse(const std::string &filename, // 구조 입력 파일
    Box &simulationBox, // 시뮬레이션 박스
    std::vector<Atom> &atoms) // 각 원자
{
    std::ifstream infile(filename); // 구조 입력 파일
    if (!infile.is_open())
    {
        throw std::runtime_error("입력 파일을 열 수 없습니다.");
    }
    
    std::string line; // 구조 입력 파일에서 각 atom에 대한 ㅈ어보
    while (std::getline(infile, line))
    {
        // find() 반환값은 size_t 타입
        if (line.find("ITEM: NUMBER OF ATOMS") == 0) // 반환된 위치가 0인지를 묻는 비교
        {
            break;
        }
    }

    size_t N = 0; // 원자 개수 N
    infile >> N;
    while (std::getline(infile, line))
    {
        if (line.find("ITEM: BOX BOUNDS") == 0)
        {
            break;
        }
    }

    // 경계 flag 추출
    {
        // istringstream : string을 입력 받아 다른 형식으로 변환
        std::istringstream iss(line); // line을 iss로 저장
        std::vector<std::string> strVec;
        std::string str;
        while (iss >> str) // string인 str로 파싱, 변수 형식에 맞게 변환
        {
            // push_back : 임시 객체에 값 복사 후, vector 마지막에 그걸 삽입
            strVec.push_back(str); // iss에 있는 애들 다 strVec에 삽입
        }

        // // pp면 periodic이므로 주기적 보정을 해야 함
        // simulationBox.px = (strVec.size()>3 && strVec[3]=="pp");
        // simulationBox.py = (strVec.size()>4 && strVec[4]=="pp");
        // simulationBox.pz = (strVec.size()>5 && strVec[5]=="pp");
    }

    // 경계
    double tmp; // 버리는 값
    infile >> simulationBox.xlo >> simulationBox.xhi >> tmp;
    infile >> simulationBox.ylo >> simulationBox.yhi >> tmp;
    infile >> simulationBox.zlo >> simulationBox.zhi >> tmp;

    // 옹스트롬 변환
    simulationBox.xlo *= 0.1;
    simulationBox.xhi *= 0.1;
    simulationBox.ylo *= 0.1;
    simulationBox.yhi *= 0.1;
    simulationBox.zlo *= 0.1;
    simulationBox.zhi *= 0.1;

    // std::cout << "[DEBUG] simulationBox.xlo="<<simulationBox.xlo<<"\n";
    // std::cout << "[DEBUG] simulationBox.xhi="<<simulationBox.xhi<<"\n";
    // std::cout << "[DEBUG] simulationBox.ylo="<<simulationBox.ylo<<"\n";
    // std::cout << "[DEBUG] simulationBox.yhi="<<simulationBox.yhi<<"\n";
    // std::cout << "[DEBUG] simulationBox.zlo="<<simulationBox.zlo<<"\n";
    // std::cout << "[DEBUG] simulationBox.zhi="<<simulationBox.zhi<<"\n";

    while (std::getline(infile, line))
    {
        if (line.rfind("ITEM: ATOMS", 0) == 0)
        {
            break;
        }
    }

    // scaled 판단
    std::istringstream headerIss(line);
    std::vector<std::string> headerCols;
    std::string col;
    while (headerIss >> col)
    {
        headerCols.push_back(col);
    }

    bool scaled = false;

    for (auto &c : headerCols)
    {
        if (c == "xs" || c == "ys" || c == "zs")
        {
            scaled = true;

            break;
        }
    }

    atoms.clear();
    // reserve : N의 용량을 미리 확보
    atoms.reserve(N);

    // 각 줄 마지막 세 열을 좌표로 하기
    for (size_t i = 0; i < N; ++i)
    {
        // getline : string 형을 받아옴
        std::getline(infile, line);
        if (line.empty())
        {
            --i; // i 다시 줄이고 continue

            continue;
        }

        std::istringstream iss(line);
        std::vector<std::string> col;
        std::string str;
        while (iss >> str)
        {
            col.push_back(str);
        }

        if (col.size() < 3)
        {
            throw std::runtime_error("atom line 형식 오류");
        }

        double xs = std::stod(col[col.size()-3]);
        double ys = std::stod(col[col.size()-2]);
        double zs = std::stod(col[col.size()-1]);

        // 입력 파일에서 xs ys zs면 scaled, x y z면 real 좌표임
        if (scaled)
        {
            xs = simulationBox.xlo + xs*(simulationBox.xhi - simulationBox.xlo);
            ys = simulationBox.ylo + ys*(simulationBox.yhi - simulationBox.ylo);
            zs = simulationBox.zlo + zs*(simulationBox.zhi - simulationBox.zlo);
        }

        atoms.push_back({ xs,ys,zs });
    }

    infile.close();
}

std::vector<double> ComputeDiameters(const Box &simulationBox,
    const std::vector<Atom> &atoms)
{
    std::vector<double> diameters;
    diameters.reserve(M);

    // double Lx = simulationBox.xhi - simulationBox.xlo; // 시뮬레이션 박스 크기 x 길이
    // double Ly = simulationBox.yhi - simulationBox.ylo; // 시뮬레이션 박스 크기 y 길이
    // double Lz = simulationBox.zhi - simulationBox.zlo; // 시뮬레이션 박스 크기 z 길이

    // TODO 시뮬레이션 박스 크기를 다공성 재료 내부로 한정해야 함
    // 지금은 전체 범위임
    std::mt19937_64 gen{std::random_device{}()}; // 몬테카를로
    std::uniform_real_distribution<double> randX(simulationBox.xlo, simulationBox.xhi),
        randY(simulationBox.ylo, simulationBox.yhi),
        randZ(simulationBox.zlo, simulationBox.zhi);

    // 혼합 파라미터
    double sigma_mix = 0.5*(sigmaAr + sigmaC);
    // double eps_mix = std::sqrt(epsilonAr * epsilonC);
    // double cutoff = cutoffFactor * sigma_mix;

    // double sigmaSquare = sigma*sigma; // 충돌 거리 제곱

    for (int i = 0; i < M; ++i)
    {
        // A 랜덤
        double Ax = randX(gen); // A의 x 성분
        double Ay = randY(gen); // A의 y 성분
        double Az = randZ(gen); // A의 z 성분
        // double Ax = 5;
        // double Ay = 5;
        // double Az = 0.4;

        // C1 C2 C3 찾기 및 φ(r_A) 체크
        double min1 = 1e9; // AC1 제곱이 될 값
        double min2 = 1e9; // AC2 제곱이 될 값
        double min3 = 1e9; // AC3 제곱이 될 값

        size_t numC1 = 0; // C1의 번호
        size_t numC2 = 0; // C2의 번호
        size_t numC3 = 0; // C3의 번호

        // TODO 혹시 시간을 더 줄이고 싶으면 tree 같은 거 써서 검색 시스템 만들기
        for (size_t j = 0; j < atoms.size(); ++j)
        {
            // A부터 C까지의 거리
            // C는 임의의 원자
            double dACx = Ax - atoms[j].x; // AC의 x 길이
            double dACy = Ay - atoms[j].y; // AC의 y 길이
            double dACz = Az - atoms[j].z; // AC의 z 길이

            // // 주기적 경계면 특정 상황에서 반대쪽 거를 쓰는 게 더 가까움
            // if (simulationBox.px)
            // {
            //     if (dACx > 0.5 * Lx)
            //     {
            //         dACx -= Lx;
            //     }

            //     if (dACx < -0.5 * Lx)
            //     {
            //         dACx += Lx;
            //     }
            // }
            // if (simulationBox.py)
            // {
            //     if (dACy > 0.5 * Ly)
            //     {
            //         dACy -= Ly;
            //     }

            //     if (dACy < -0.5 * Ly)
            //     {
            //         dACy += Ly;
            //     }
            // }
            // if (simulationBox.pz)
            // {
            //     if (dACz > 0.5 * Lz)
            //     {
            //         dACz -= Lz;
            //     }

            //     if (dACz < -0.5 * Lz)
            //     {
            //         dACz += Lz;
            //     }
            // }

            double dACsquare = dACx*dACx + dACy*dACy + dACz*dACz; // AC 길이 제곱

            // 거리에 따라 C1 C2 C3 설정
            if (dACsquare < min1)
            {
                min3 = min2;
                numC3 = numC2;
                min2 = min1;
                numC2 = numC1;
                min1 = dACsquare;
                numC1 = j;
            }
            else if (dACsquare < min2)
            {
                min3 = min2;
                numC3 = numC2;
                min2 = dACsquare;
                numC2 = j;
            }
            else if (dACsquare < min3)
            {
                min3 = dACsquare;
                numC3 = j;
            }
        }

        // std::cout << "[DEBUG] min1="<<min1<<"\nmin2="<<min2<<"\nmin3="<<min3<<"\n";

        // min1을 통해 성공 여부 결정
        if (min1 < sigma_mix * sigma_mix) // 퍼텐셜이 0 이상이면 실패
        {
            // std::cout << "[DEBUG] min1="<<min1<<" < sigma_mix * sigma_mix"<<sigma_mix * sigma_mix<<"\n";
            continue;
        }

        // numC1 numC2 numC3번째 atoms를 각각 C1 C2 C3로 설정
        const Atom &C1 = atoms[numC1];
        const Atom &C2 = atoms[numC2];
        const Atom &C3 = atoms[numC3];

        double dAC1x = Ax - C1.x; // AC1의 x 길이
        double dAC1y = Ay - C1.y; // AC1의 ㅛ 길이
        double dAC1z = Az - C1.z; // AC1의 z 길이
        double dAC2x = Ax - C2.x; // AC2의 x 길이
        double dAC2y = Ay - C2.y; // AC2의 ㅛ 길이
        double dAC2z = Az - C2.z; // AC2의 z 길이

        double dot12 = dAC1x*dAC2x + dAC1y*dAC2y + dAC1z*dAC2z;
        double lam1 = 0.0;
        double denom1 = (min1 - dot12);
        if (fabs(denom1) > 1e-12)
        {
            lam1 = (min2 - min1) / (2.0 * denom1);
        }

        // B 계산
        double Bx = Ax + lam1 * dAC1x; // B의 x 성분
        double By = Ay + lam1 * dAC1y; // B의 y 성분
        double Bz = Az + lam1 * dAC1z; // B의 z 성분

        // rCB 벡터
        double dBC1x = Bx - C1.x; // BC1의 x 길이
        double dBC1y = By - C1.y; // BC1의 y 길이
        double dBC1z = Bz - C1.z; // BC1의 z 길이
        double dBC2x = Bx - C2.x; // BC2의 x 길이
        double dBC2y = By - C2.y; // BC2의 y 길이
        double dBC2z = Bz - C2.z; // BC2의 z 길이
        double dBC3x = Bx - C3.x; // BC3의 x 길이
        double dBC3y = By - C3.y; // BC3의 y 길이
        double dBC3z = Bz - C3.z; // BC3의 z 길이

        double dBC1pBC2x = dBC1x + dBC2x; // BC1 x 길이 + BC2 x 길이
        double dBC1pBC2y = dBC1y + dBC2y; // BC1 y 길이 + BC2 y 길이
        double dBC1pBC2z = dBC1z + dBC2z; // BC1 z 길이 + BC2 z 길이

        double dBC1Square = dBC1x * dBC1x + dBC1y * dBC1y + dBC1z * dBC1z; // BC1 길이 제곱
        double dBC3Square = dBC3x * dBC3x + dBC3y * dBC3y + dBC3z * dBC3z; // BC3 길이 제곱

        double dot1p = dBC1x * dBC1pBC2x + dBC1y * dBC1pBC2y + dBC1z * dBC1pBC2z;
        double dot3p = dBC3x * dBC1pBC2x + dBC3y * dBC1pBC2y + dBC3z * dBC1pBC2z;
        double lam2 = 0.0;
        double denom2 = (dot1p - dot3p);
        if (fabs(denom2) > 1e-12)
        {
            lam2 = (dBC3Square - dBC1Square) / (2.0 * denom2);
        }
        
        // D 계산
        double Dx = Bx + lam2 * dBC1pBC2x; // D의 x 성분
        double Dy = By + lam2 * dBC1pBC2y; // D의 y 성분
        double Dz = Bz + lam2 * dBC1pBC2z; // D의 z 성분

        auto dist = [&](const Atom &C){double vx = Dx - C.x;
            double vy = Dy - C.y;
            double vz = Dz - C.z;
            return std::sqrt(vx*vx + vy*vy + vz*vz);
        };

        double dC1 = dist(C1);
        double dC2 = dist(C2);
        double dC3 = dist(C3);

        // std::cout << "[DEBUG] dC1="<<dC1<<"\ndC2="<<dC2<<"\ndC3="<<dC3<<"\n";

        double radius = dC2 - sigma_mix; // 반경 radius

        if (radius <= 0)
        {
            // std::cout << "[DEBUG] radius="<<radius<<" <= 0\n";
            continue;
        }

        double diameter = 2.0 * radius; // 지름 diameter
        
        diameters.push_back(diameter);

        // std::cout<<"[DEBUG] D="<<diameter<<"\n\n";
        std::cout<<"[DEBUG] Ax="<<Ax<<"\n"<<"Ay="<<Ay<<"\n"<<"Az="<<Az<<"\n"
            <<"min1="<<min1<<"\n"
            <<"Dx="<<Dx<<"\n"<<"Dy="<<Dy<<"\n"<<"Dz="<<Dz<<"\n"
            <<"dC1="<<dC1<<"\n"<<"dC2="<<dC2<<"\n"<<"dC3="<<dC3<<"\n"
            <<"D="<<diameter<<"\n";
    }

    return diameters;
}

std::vector<std::pair<double,double>> ComputeAccessibleVolume(const std::vector<double> &diameters,
    double Vbox)
{
    std::vector<int> Mj(bin, 0);

    for (double D : diameters)
    {
        if (D <= 0 || D > Dmax)
        {
            continue;
        }

        int idx = std::min((int)(D / Dmax * bin), bin - 1);

        Mj[idx]++;
    }

    std::vector<std::pair<double,double>> dist;
    dist.reserve(bin);

    double dV = 1 / (double)M;

    for (int j = 0; j < bin; ++j)
    {
        double diamMid = (j + 0.5) * (Dmax / bin);
        double Vacc_j = Mj[j] * dV;

        // emplace_back : 임시 객체 없이, vector 마지막에 바로 삽입
        dist.emplace_back(diamMid, Vacc_j);
    }

    return dist;
}

void WriteCSV(const std::vector<std::pair<double,double>> &distribution,
    const std::string &outputFilename)
{
    std::ofstream ofs(outputFilename);
    if (!ofs)
    {
        throw std::runtime_error("출력 파일을 열 수 없습니다.");
    }

    ofs << "diameter(nm), accessible volume(nm3)\n";
    ofs << std::fixed;
    ofs.precision(6);

    for (auto &dv : distribution)
    {
        ofs << dv.first << ", " << dv.second << "\n";
    }

    ofs.close();
}