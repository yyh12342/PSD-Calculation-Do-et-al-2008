# PSD-Calculation-Do-et-al-2008  

분자 정보를 기반으로 기공 크기와 부피 분포를 탐색  

실제 기공은 불균일한 모양을 가지므로, accessible volume으로 계산  
accessible volume: 구 형태의 최대 volume  
perfect slit pore일 때 accessible volume이 0이 될 수 있음  
accessible volume의 합계를 총 기공 부피로 함  

probe 입자로 아르곤 사용하며 accessible volume에서 0 이하의 퍼텐셜을 가짐  
시뮬레이션 박스의 임의의 위치에 아르곤 입자를 m번 삽입하여 accessible volume을 결정  
기공 크기 상한은 고체에 따른 최대 기공 크기에 달려있음  

단순히 계산하면 accessible volume은 삽입 성공 확률 f에 box volume 곱한 값  

세 고체 원자 C의 좌표를 통해 accessible volume 구해야 함  
C1: 퍼텐셜이 0 이하인 무작위 점 A에서 가장 가까운 고체 원자  
C2: 점 A에서 두번째로 가까운 고체 원자  
벡터 C1A 상에서 C1까지의 거리와 C2 까지의 거리가 같도록 점 B 설정  
C3: 점 A에서 세번째로 가까운 고체 원자  
벡터 C1B + 벡터 C2B 상에서 C1, C2, C3 까지의 거리가 서로 같도록 점 D 설정  
벡터 C1D나 벡터 C2D 또는는 벡터 C3D에서 퍼텐셜이 0인 점 탐색  
세 점 중 하나를 지나며 점 D를 중심으로 하는 가장 큰 구를 선택하면, 그 내부가 0 이하의 퍼텐셜 (점 A도 내부에 포함됨)  

개별 accessible volume의 합이 총 accessible volume  
이 accessible volume은 말 그대로 접근 가능한 볼륨이며 실질적인 물리적 직경보다 작음  

퍼텐셜의 계산은 LJ potential로  



# 알고리듬  

1. 시뮬레이션 박스 선택  
2. 기공 크기 범위를 선택하고 bin으로 나눔  
3. 점 A 선택 후 퍼텐셜이 0 이하인지 확인  
4. C1 선택  
5. C2 선택하고 점 B 선택  
6. C3 선택하고 점 D 선택  
7. 벡터 C1D, C2D, C3D에서 각각 퍼텐셜이 0인 점 선택  
8. 세 점 중에서 D까지의 거리가 가장 먼 점 선택  
9. 이 점을 지나고 점 D를 중심으로 하는 구가 accessible volume  



# Example  

1. 이중 슬릿 기공  
2. 탄소나노튜브  



# 코드 구조  

- PDBParser  
.pdb 파일을 읽고, 원소와 좌표를 반환  
처음엔 범용적으로 사용할 수 있도록 코드를 설계했는데, 생각해보니 Example 1과 2의 환경이 평행한 이중 슬릿과 탄소나노튜브로 상이함  
따라서 특정 상황에서 적용할 수 있도록 코드를 단순화하였고, 이에 따라 PDBParser는 당분간 사용하지 않음  

- Potnetial  
파라미터, 퍼텐셜 계산  

- SphereFinder  
세 좌표를 통해 기공 반경 계산  

- Histogram  
각 지름 bin 별로 volume 저장  



# 빌드 방법  

mkdir build  
cd build  
cmake .. -G "Visual Studio 17 2022" -A x64  

cmake --build . --config Release  



# 실행 방법  

- PSDCalculator  
cd C:\Users\yooyh\Project\PSD-Calculation-Do-et-al-2008  
build\Release\PSDCalculator.exe  

- plot  
py plot_psd.py  
이 코드에서 csv와 png 파일명 수정 가능  

파이썬 버전은 3.12.7  



# 빌드 시 오류 해결  

- 한글 깨질 때  
vscode 실행 후, 오른쪽 하단 인코딩에서 UTF-8 with BOM  



# 입력값  

PSDCalculator에 넣을 테스트 입력값  

1  
0.3405  
119.8  
0.340  
28.0  
0.142  
5 5 4  
1.0 2.0  
100000  

2  
0.3405  
119.8  
0.340  
28.0  
0.142  