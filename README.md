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
지름 계산  
4. M번 반복 (약 100만)  
5. Vacc 계산  
dV = Vbox / M  
Vacc_j = Mj * dV (각 bin에 대해)  



# Example  

1. graphite 이중 슬릿 기공  
2. cnt 탄소나노튜브  



# 코드 구조  

- accessible_volume  
- main  
- plot_psd : 플롯 생성  
- dump_to_numpy : openPNM 검증을 위해 dump 파일을 변환  



# 빌드 방법  

cd build  
cmake .. -G "Visual Studio 17 2022" -A x64  
cmake --build . --config Release  



# 실행 방법  

cd ..  
build\Release\PSDCalculator.exe inputs/dump.graphite1 graphite results/result_ex1_1.csv  

(PSDCalculator.exe {input 파일명명} {파일 형식} {csv 파일명})  

- 파일명  
Ex1 : inputs/dump.graphite1  
Ex2 : inputs/dump.cnt.armchair.r7.largebox.rlx  

- 파일 형식  
Ex1 : graphite (graphite뿐만 아니라 원자들이 한 평면에 존재하는 구조면 graphite 형식으로 적용 가능)  
Ex2 : cnt (범용적인 구조와 계산 과정은 같지만, 점 A를 cnt 내부에 위치시키는 제약 조건 존재)  
default : 범용적인 구조  



# plotting 방법  

py plot_psd.py results/result_ex1_1.csv results/plots/plot_ex1_1.png  

(plot_psd.py {csv 파일명} {설정할 plot 파일명})  

파이썬 버전은 3.12.7  



# 빌드 시 오류 해결  

- 한글 깨질 때  
vscode 실행 후, 오른쪽 하단 인코딩에서 UTF-8 with BOM  



# 참고  

- openPNM  
- zeo++  