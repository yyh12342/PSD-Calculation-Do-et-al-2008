#define _USE_MATH_DEFINES

#include "accessible_volume.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <random>
#include <algorithm>

static const double sigmaAr = 0.3405; // 아르곤 충돌 지름 σ [nm]
static const double sigmaC = 0.3400; // 탄소 충돌 지름 σ [nm]

static const int M = 1000000; // Monte Carlo 시도 횟수
// static const int M = 10; // 디버깅용 Monte Carlo 시도 횟수
static const double Dmax = 2.0; // 지름 분포 최대 직경 [nm]
static const int bin = 50; // bin 개수

// cnt 모드일 때 사용
int hAxis = 0; // cnt 높이 축
int rAxis1 = 0; // cnt 지름 축 1
int rAxis2 = 0; // cnt 지름 축 2
double center[3]; // cnt 중심
double rSquare = 0.0; // cnt 반경 제곱

void Parse(
    const std::string &filename, // 구조 입력 파일
    Box &simulationBox, // 시뮬레이션 박스
    std::vector<Atom> &atoms,
    int modeNum) // 각 원자
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
    }

    // BOX BOUNDS
    // 범용적인 구조가 graphite일 때와 같은 BOX BOUNDS 형식을 가진다고 가정
    // TODO 만약 그렇지 않다면 수정 필요
    if (modeNum == 1)
    {
        infile >> simulationBox.xlo >> simulationBox.xhi;
        infile >> simulationBox.ylo >> simulationBox.yhi;
        infile >> simulationBox.zlo >> simulationBox.zhi;
    }
    else // graphite 모드, 범용적인 구조
    {
        double tmp; // 버리는 값
        infile >> simulationBox.xlo >> simulationBox.xhi >> tmp;
        infile >> simulationBox.ylo >> simulationBox.yhi >> tmp;
        infile >> simulationBox.zlo >> simulationBox.zhi >> tmp;
    }

    // 옹스트롬 변환
    simulationBox.xlo *= 0.1;
    simulationBox.xhi *= 0.1;
    simulationBox.ylo *= 0.1;
    simulationBox.yhi *= 0.1;
    simulationBox.zlo *= 0.1;
    simulationBox.zhi *= 0.1;

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
        // TODO 혹시나 세 축 중 일부만 scaled면, 각각 설정해줘야 함
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

std::vector<double> ComputeDiameters(
    const Box &simulationBox,
    const std::vector<Atom> &atoms,
    int modeNum)
{
    std::vector<double> diameters;
    diameters.reserve(M);

    if (modeNum == 1) // cnt일 때 실린더의 축 결정
    {
        // 원자가 분포되어 있는 길이가 가장 긴 축이 높이 축이라고 가정
        // TODO 만약 그렇지 않은 구조라면 수정해야 함
        double minX = 1e9;
        double maxX = -1e9;
        double minY = 1e9;
        double maxY = -1e9;
        double minZ = 1e9;
        double maxZ = -1e9;

        for (auto &A : atoms)
        {
            minX = std::min(minX, A.x);
            maxX = std::max(maxX, A.x);
            minY = std::min(minY, A.y);
            maxY = std::max(maxY, A.y);
            minZ = std::min(minZ, A.z);
            maxZ = std::max(maxZ, A.z);
        }

        double spanX = maxX - minX;
        double spanY = maxY - minY;
        double spanZ = maxZ - minZ;

        if (spanX > spanY && spanX > spanZ)
        {
            hAxis = 0;
        }
        else if (spanY > spanX && spanY > spanZ)
        {
            hAxis = 1;
        }
        else
        {
            hAxis = 2;
        }

        rAxis1 = (hAxis + 1) % 3;
        rAxis2 = (hAxis + 2) % 3;

        double min1 = 1e9;
        double max1 = -1e9;
        double min2 = 1e9;
        double max2 = -1e9;

        for (auto &A : atoms)
        {
            double c1 = (&A.x)[rAxis1];
            double c2 = (&A.x)[rAxis2];

            min1 = std::min(min1, c1);
            max1 = std::max(max1, c1);
            min2 = std::min(min2, c2);
            max2 = std::max(max2, c2);
        }

        center[rAxis1] = 0.5 * (min1 + max1);
        center[rAxis2] = 0.5 * (min2 + max2);

        for (auto &A : atoms)
        {
            double dr1 = (&A.x)[rAxis1] - center[rAxis1];
            double dr2 = (&A.x)[rAxis2] - center[rAxis2];

            rSquare = std::max(rSquare, dr1*dr1 + dr2*dr2);
        }
    }

    std::mt19937_64 gen{std::random_device{}()}; // 몬테카를로
    std::uniform_real_distribution<double> randX(simulationBox.xlo, simulationBox.xhi);
    std::uniform_real_distribution<double> randY(simulationBox.ylo, simulationBox.yhi);
    std::uniform_real_distribution<double> randZ(simulationBox.zlo, simulationBox.zhi);

    double sigma_mix = 0.5*(sigmaAr + sigmaC); // 혼합 파라미터

    const double ztol = 0.1; // graphite 모드일 때 그래핀 층 구분
    int accepted = 0; // A가 다공성 재료 내부에 삽입된 횟수

    while (accepted < M)
    {
        // A 랜덤
        double Ax = randX(gen); // A의 x 성분
        double Ay = randY(gen); // A의 y 성분
        double Az = randZ(gen); // A의 z 성분

        // cnt 내부 삽입 여부 체크
        double coords[3] = { Ax, Ay, Az };

        if (modeNum == 1)
        {
            double dr1 = coords[rAxis1] - center[rAxis1];
            double dr2 = coords[rAxis2] - center[rAxis2];

            if (dr1*dr1 + dr2*dr2 > rSquare)
            {
                continue; // cnt 외부에 있으면 continue
            }

            // std::cout<<"[DEBUG] dr1*dr1 + dr2*dr2="<<dr1*dr1 + dr2*dr2<<"\n"<<
            //     "rSquare="<<rSquare<<"\n\n";
        }

        // C1 C2 C3 찾기 및 퍼텐셜 체크
        double min1 = 1e9; // AC1 제곱이 될 값
        double min2 = 1e9; // AC2 제곱이 될 값
        double min3 = 1e9; // AC3 제곱이 될 값

        size_t numC1 = 0; // C1의 번호
        size_t numC2 = SIZE_MAX; // C2의 번호
        size_t numC3 = SIZE_MAX; // C3의 번호

        for (size_t j = 0; j < atoms.size(); ++j)
        {
            double dCAx = Ax - atoms[j].x; // AC의 x 길이
            double dCAy = Ay - atoms[j].y; // AC의 y 길이
            double dCAz = Az - atoms[j].z; // AC의 z 길이

            double dCASquare = dCAx*dCAx + dCAy*dCAy + dCAz*dCAz; // AC 길이 제곱

            if (dCASquare < min1)
            {
                min1 = dCASquare;
                numC1 = j;
            }
        }

        // min1을 통해 성공 여부 결정
        if (min1 < sigma_mix * sigma_mix) // 퍼텐셜이 0 이상이면 실패
        {
            ++accepted;

            // std::cout<<"[DEBUG] min1="<<min1<<"\n"<<
            //     "sigma_mix * sigma_mix="<<sigma_mix * sigma_mix<<"\n\n";

            continue;
        }

        if (modeNum == 0) // graphite 모드
        {
            for(size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == numC1)
                {
                    continue;
                }

                if (std::fabs(atoms[j].z - atoms[numC1].z) < ztol)
                {
                    continue;
                }

                double dCAx = Ax - atoms[j].x; // AC의 x 길이
                double dCAy = Ay - atoms[j].y; // AC의 y 길이
                double dCAz = Az - atoms[j].z; // AC의 z 길이

                double dCASquare = dCAx*dCAx + dCAy*dCAy + dCAz*dCAz; // AC 길이 제곱

                if (dCASquare < min2)
                {
                    min2 = dCASquare;
                    numC2 = j;
                }
            }

            if (numC2 == SIZE_MAX)
            {
                for(size_t j = 0; j < atoms.size(); ++j)
                {
                    if (j == numC1)
                    {
                        continue;
                    }

                    double dCAx = Ax - atoms[j].x; // AC의 x 길이
                    double dCAy = Ay - atoms[j].y; // AC의 y 길이
                    double dCAz = Az - atoms[j].z; // AC의 z 길이

                    double dCASquare = dCAx*dCAx + dCAy*dCAy + dCAz*dCAz; // AC 길이 제곱

                    if (dCASquare < min2)
                    {
                        min2 = dCASquare;
                        numC2 = j;
                    }
                }
            }

            for(size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == numC1 || j == numC2)
                {
                    continue;
                }

                double dCAx = Ax - atoms[j].x; // AC의 x 길이
                double dCAy = Ay - atoms[j].y; // AC의 y 길이
                double dCAz = Az - atoms[j].z; // AC의 z 길이

                double dCASquare = dCAx*dCAx + dCAy*dCAy + dCAz*dCAz; // AC 길이 제곱

                if (dCASquare < min3)
                {
                    min3 = dCASquare;
                    numC3 = j;
                }
            }
        }
        else // cnt 모드, 범용적인 구조
        {
            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == numC1)
                {
                    continue;
                }

                double dCAx = Ax - atoms[j].x; // AC의 x 길이
                double dCAy = Ay - atoms[j].y; // AC의 y 길이
                double dCAz = Az - atoms[j].z; // AC의 z 길이

                double dCASquare = dCAx*dCAx + dCAy*dCAy + dCAz*dCAz; // AC 길이 제곱

                if (dCASquare < min2)
                {
                    min3 = min2;
                    numC3 = numC2;
                    min2 = dCASquare;
                    numC2 = j;
                }
                else if (dCASquare < min3)
                {
                    min3 = dCASquare;
                    numC3 = j;
                }
            }
        }

        // numC1 numC2 numC3번째 atoms를 각각 C1 C2 C3로 설정
        const Atom &C1 = atoms[numC1];
        const Atom &C2 = atoms[numC2];
        const Atom &C3 = atoms[numC3];

        double dC1Ax = Ax - C1.x; // AC1의 x 길이
        double dC1Ay = Ay - C1.y; // AC1의 ㅛ 길이
        double dC1Az = Az - C1.z; // AC1의 z 길이
        double dC2Ax = Ax - C2.x; // AC2의 x 길이
        double dC2Ay = Ay - C2.y; // AC2의 ㅛ 길이
        double dC2Az = Az - C2.z; // AC2의 z 길이
        double dC1C2x = C2.x - C1.x;
        double dC1C2y = C2.y - C1.y;
        double dC1C2z = C2.z - C1.z;

        double dC1ASquare = dC1Ax*dC1Ax + dC1Ay*dC1Ay + dC1Az*dC1Az;
        double dC2ASquare = dC2Ax*dC2Ax + dC2Ay*dC2Ay + dC2Az*dC2Az;

        double denom1 = dC1C2x*dC1Ax + dC1C2y*dC1Ay + dC1C2z*dC1Az;
        if (fabs(denom1) < 1e-12) // A가 수직 위에 있는 경우 점 B가 존재하지 못함
        {
            continue;
        }

        double lam1 = (dC2ASquare - dC1ASquare) / (2 * denom1);

        // B 계산
        double Bx = Ax + lam1 * dC1Ax; // B의 x 성분
        double By = Ay + lam1 * dC1Ay; // B의 y 성분
        double Bz = Az + lam1 * dC1Az; // B의 z 성분

        // rCB 벡터
        double dC1Bx = Bx - C1.x; // BC1의 x 길이
        double dC1By = By - C1.y; // BC1의 y 길이
        double dC1Bz = Bz - C1.z; // BC1의 z 길이
        double dC2Bx = Bx - C2.x; // BC2의 x 길이
        double dC2By = By - C2.y; // BC2의 y 길이
        double dC2Bz = Bz - C2.z; // BC2의 z 길이
        double dC3Bx = Bx - C3.x; // BC3의 x 길이
        double dC3By = By - C3.y; // BC3의 y 길이
        double dC3Bz = Bz - C3.z; // BC3의 z 길이
        double dC1C3x = C3.x - C1.x;
        double dC1C3y = C3.y - C1.y;
        double dC1C3z = C3.z - C1.z;

        double dC1BSquare = dC1Bx * dC1Bx + dC1By * dC1By + dC1Bz * dC1Bz; // BC1 길이 제곱
        double dC3BSquare = dC3Bx * dC3Bx + dC3By * dC3By + dC3Bz * dC3Bz; // BC3 길이 제곱

        double denom2 = dC1C3x*(dC1Bx+dC2Bx) + dC1C3y*(dC1By+dC2By) + dC1C3z*(dC1Bz+dC2Bz);
        if (fabs(denom2) < 1e-12)
        {
            continue;
        }

        double lam2 = (dC3BSquare - dC1BSquare) / (2 * denom2);
        
        // D 계산
        double Dx = Bx + lam2 * (dC1Bx+dC2Bx); // D의 x 성분
        double Dy = By + lam2 * (dC1By+dC2By); // D의 y 성분
        double Dz = Bz + lam2 * (dC1Bz+dC2Bz); // D의 z 성분

        auto dist = [&](const Atom &C){
            double vx = Dx - C.x;
            double vy = Dy - C.y;
            double vz = Dz - C.z;
            return std::sqrt(vx*vx + vy*vy + vz*vz);
        };

        double dC1 = dist(C1);
        double dC2 = dist(C2);
        double dC3 = dist(C3);

        double radius = std::min({dC1, dC2, dC3}) - sigma_mix; // 반경 radius

        if (radius > 0)
        {
            double diameter = 2.0 * radius; // 지름 diameter

            diameters.push_back(diameter);

            // std::cout<<"[DEBUG] Ax="<<Ax<<"\n"<<"Ay="<<Ay<<"\n"<<"Az="<<Az<<"\n"
            //     <<"min1="<<min1<<"\n"
            //     <<"Dx="<<Dx<<"\n"<<"Dy="<<Dy<<"\n"<<"Dz="<<Dz<<"\n"
            //     <<"dC1="<<dC1<<"\n"<<"dC2="<<dC2<<"\n"<<"dC3="<<dC3<<"\n"
            //     <<"D="<<diameter<<"\n\n";
        }

        ++accepted;
    }

    // std::cout<<"[DEBUG] accepted="<<accepted<<"\n\n";

    return diameters;
}

std::vector<std::pair<double,double>> ComputeAccessibleVolume(
    const std::vector<double> &diameters,
    const Box &simulationBox,
    int modeNum)
{
    double Vbox;

    if (modeNum == 1) // cnt
    {
        double dx = simulationBox.xhi - simulationBox.xlo;
        double dy = simulationBox.yhi - simulationBox.ylo;
        double dz = simulationBox.zhi - simulationBox.zlo;

        double radius = std::sqrt(rSquare);
        double height = (hAxis == 0 ? dx : (hAxis == 1 ? dy : dz));

        Vbox = M_PI * radius * radius * height;
    }
    else // graphite와 범용적인 구조
    {
        Vbox = std::fabs(simulationBox.xhi - simulationBox.xlo)
            * std::fabs(simulationBox.yhi - simulationBox.ylo)
            * std::fabs(simulationBox.zhi - simulationBox.zlo);
    }
    
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

void WriteCSV(
    const std::vector<std::pair<double,double>> &distribution,
    const std::string &outputFilename)
{
    std::ofstream ofs(outputFilename);
    if (!ofs)
    {
        throw std::runtime_error("출력 파일을 열 수 없습니다.");
    }

    ofs << "diameter(nm), reduced accessible volume\n";
    ofs << std::fixed;
    ofs.precision(6);

    for (auto &dv : distribution)
    {
        ofs << dv.first << ", " << dv.second << "\n";
    }

    ofs.close();
}