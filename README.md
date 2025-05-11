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



# 알고리즘  

1. 시뮬레이션 박스 설정  
텍스트(아마 LAMMPS)로 입력  
2. 지름 분포를 bin j로 나눔  
처음에는 모든 bin에서 Vacc가 0  
성공 횟수 Mj를 저장  
3. 랜덤 위치 rA로 퍼텐셜 phi(rA) 계산  
퍼텐셜이 0 이하면 성공  
(a) 제일 가까운 원자 C1 찾기  
(b) 두번째로 가까운 원자 C2 찾기  
C1, C2까지의 거리가 같은 B 설정  
세번쨰로 가까운 원자 C3 찾기  
C1, C2, C3까지의 거리가 같은 D 설정  
(C) D에서 C1, C2, C3까지 가다가 퍼텐셜이 0인 점 찾기  
이 중에서 가운데 값인 C2를 사용하자  
지름 계산  
4. M번 반복 (약 100만)  
5. Vacc 계산  
dV = Vbox / M  
Vacc_j = Mj * dV (각 bin에 대해)  



# Example  

1. graphite 이중 슬릿 기공  
2. cnt 탄소나노튜브  



# 코드 구조  

아래 코드들 이제 안 쓰고 accessible_volume만 씀  

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

cd build  
cmake .. -G "Visual Studio 17 2022" -A x64  
cmake --build . --config Release  



# 실행 방법  

cd ..  
build\Release\PSDCalculator.exe inputs/dump.graphite1 graphite results/result_ex1_1.csv  

(PSDCalculator.exe {input 파일 이름} {파일 형식} {csv 파일 이름})  

- 파일명  
Ex1 : inputs/dump.graphite1  
Ex2 : inputs/dump.cnt.armchair.r7.largebox.rlx  

- 파일 형식  
Ex1 : graphite  
Ex2 : cnt  



# plotting 방법  

py plot_psd.py  

이 코드에서 입력할 csv와 png 파일명 수정 가능  

파이썬 버전은 3.12.7  



# 빌드 시 오류 해결  

- 한글 깨질 때  
vscode 실행 후, 오른쪽 하단 인코딩에서 UTF-8 with BOM  