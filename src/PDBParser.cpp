// #include "PDBParser.h"

// #include <fstream>
// #include <sstream>
// #include <cctype>
// #include <algorithm>
// #include <iostream>

// // PDB 분석
// bool PDBParser::parseFile(const std::string &filename, // pdb 파일 경로
//     std::vector<Atom> &atoms, // 원자 리스트를 반환할 참조
//     Box &box) // 시뮬레이션한 박스 좌표 반환할 참조
// {
//     std::ifstream file(filename); // 파일 입력

//     if (!file)
//     {
//         std::cerr << "Error: Could not open PDB file " << filename << std::endl;
        
//         return false;
//     }

//     atoms.clear();

//     // 박스 좌표 초기화
//     // TODO 임의로 volume 크기를 2nm로 세팅했는데, 나중에 바꿔야 함
//     box.x_min = box.y_min = box.z_min = 1e9;
//     box.x_max = box.y_max = box.z_max = -1e9;
//     std::string line;
//     // CRYST1: 단위 셀 정보 (unit cube)
//     double box_a=0, box_b=0, box_c=0; // 위치
//     double box_alpha=90, box_beta=90, box_gamma=90; // 각도
//     bool gotCryst = false; // CRYST1를 찾을 수 있는지 확인

//     while (std::getline(file, line)) // 파일의 모든 줄을 한 줄씩 반복
//     {
//         // 파일 읽기
//         if (line.rfind("CRYST1", 0) == 0) // CRYST1로 시작하면
//         {
//             // PDB 포맷은 고정폭 텍스트이므로, xyz 좌표열도 고정되어 있음
//             // TODO PDB 사용하지 않을 경우 수정해야 함
//             if (line.size() >= 54) // 파일이 충분한 길이인지 확인
//             {
//                 try
//                 {
//                     std::string a_str = line.substr(6, 9); // 문자열 추출 (6번쨰부터 9개를 추출)
//                     std::string b_str = line.substr(15, 9);
//                     std::string c_str = line.substr(24, 9);
//                     std::string alpha_str = line.substr(33, 7);
//                     std::string beta_str = line.substr(40, 7);
//                     std::string gamma_str = line.substr(47, 7);

//                     box_a = std::stod(a_str); // string to double
//                     box_b = std::stod(b_str);
//                     box_c = std::stod(c_str);
//                     if (!alpha_str.empty()) box_alpha = std::stod(alpha_str);
//                     if (!beta_str.empty()) box_beta = std::stod(beta_str);
//                     if (!gamma_str.empty()) box_gamma = std::stod(gamma_str);

//                     gotCryst = true;
//                 }
//                 catch(...) // 나머지 모든 예외
//                 {
//                     gotCryst = false;
//                 }
//             }
//         }

//         // PDB에는 일반 원자와 비표준 원자로, 두 가지 원자 레코드가 있음
//         // 즉 모든 원자를 탐색한다는 뜻
//         if (line.rfind("ATOM", 0) == 0 ||
//             line.rfind("HETATM", 0) == 0) // 헤테로 원자
//         {
//             // 원자 좌표 추출
//             std::string xs, ys, zs; // 원자 좌표

//             if (line.size() >= 54) // 파일이 충분한 길이일 때, 고정 폭으로 문자열 추출
//             {
//                 xs = line.substr(30, 8);
//                 ys = line.substr(38, 8);
//                 zs = line.substr(46, 8);
//             }
//             else
//             {
//                 std::istringstream iss(line); // 각 타입 별로 추출
//                 std::string token;
//                 int col = 0;

//                 while (iss >> token) // string 타입인 token을 전부 추출
//                 {
//                     col++;

//                     // 6, 7, 8번째 문자열이 각각 원자 좌표 x, y, z 나타내도록 했는데, 나중에 바꿔야 함
//                     if (col == 6) xs = token;
//                     if (col == 7) ys = token;
//                     if (col == 8) zs = token;
//                 }
//             }

//             if (xs.empty() || ys.empty() || zs.empty()) // 파싱 실패하면
//             {
//                 continue;
//             }

//             double x = std::stod(xs);
//             double y = std::stod(ys);
//             double z = std::stod(zs);

//             // TODO 단위가 옹스트롬임을 상정하고 nm로 단위 변환했는데, 나중에 필요 시 수정
//             x *= 0.1;
//             y *= 0.1;
//             z *= 0.1;

//             // 원소 기호 추출
//             std::string elem; // 원소 기호 저장할 문자열

//             // element 필드가 77, 78열
//             if (line.size() >= 78)
//             {
//                 elem = line.substr(76, 2);
//             }

//             elem.erase(std::remove_if(elem.begin(), elem.end(), ::isspace), elem.end()); // 앞뒤 공백 제거

//             if (elem.empty()) // element 열이 없으면
//             {
//                 if (line.size() >= 16) // 원자 이름은 13~16 열
//                 {
//                     std::string name = line.substr(12, 4); // atom name 열에서 유추
//                     name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end());
//                     name.erase(std::remove_if(name.begin(), name.end(), ::isdigit), name.end());
                    
//                     if (!name.empty())
//                     {
//                         if (name.size() >= 2 && std::isupper(name[0]) && std::isupper(name[1])) // isupper: 대문자인지 확인
//                         {
//                             elem = name.substr(0, 2);
//                         }
//                         else
//                         {
//                             elem = name.substr(0, 1);
//                         }
//                     }
//                 }
//             }

//             // 원소 기호 첫번째 문자는 대문자, 두번째 문자는 소문자로
//             if (elem.size() == 1)
//             {
//                 elem[0] = std::toupper(elem[0]);
//             }
//             else if (elem.size() == 2)
//             {
//                 elem[0] = std::toupper(elem[0]);
//                 elem[1] = std::tolower(elem[1]);
//             }
            
//             // Atom 객체 생성
//             Atom atom;
//             atom.element = elem;
//             atom.position = {x, y, z};
//             atoms.push_back(atom); // 임시 객체에 값 복사 후, vector에 삽입

//             // 박스 좌표 업데이트
//             if (x < box.x_min)
//             {
//                 box.x_min = x;
//             }
//             if (x > box.x_max)
//             {
//                 box.x_max = x;
//             }
//             if (y < box.y_min)
//             {
//                 box.y_min = y;
//             }
//             if (y > box.y_max)
//             {
//                 box.y_max = y;
//             }
//             if (z < box.z_min)
//             {
//                 box.z_min = z;
//             }
//             if (z > box.z_max)
//             {
//                 box.z_max = z;
//             }
//         }
//     }

//     file.close();

//     if (atoms.empty())
//     {
//         std::cerr << "Error: No atoms found in PDB file." << std::endl;
        
//         return false;
//     }

//     if (gotCryst)
//     {
//         if (std::fabs(box_alpha - 90.0) < 1e-2 &&
//             std::fabs(box_beta - 90.0) < 1e-2 &&
//             std::fabs(box_gamma - 90.0) < 1e-2)
//         {
//             box.x_min = 0.0;
//             box.y_min = 0.0;
//             box.z_min = 0.0;
//             box.x_max = box_a * 0.1;
//             box.y_max = box_b * 0.1;
//             box.z_max = box_c * 0.1;
//         }
//         else
//         {
//             std::cerr << "Warning: Non-orthogonal unit cell angles detected; using atomic bounding box instead." << std::endl;
//         }
//     }

//     return true;
// }