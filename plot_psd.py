import pandas as pd
import matplotlib.pyplot as plt
import os

# results에서 csv 파일 읽고 그래프 생성하기 위한 코드

df = pd.read_csv("results/psd_data.csv") # 여기서 csv 파일명 변경
total = df['Volume(nm^3)'].sum()
df['ReducedAccessibleVolume'] = df['Volume(nm^3)'] / total

# 디렉토리 생성
os.makedirs("results/plots", exist_ok=True)

# 그래프
plt.figure()
plt.plot(df['Diameter(nm)'], df['ReducedAccessibleVolume'], '-o')
plt.xlabel('Diameter (nm)')
plt.ylabel('Reduced Accessible Volume')
plt.title('Pore Size Distribution')
plt.grid(True)
plt.tight_layout()
plt.savefig("results/psd_plot.png") # 여기서 png 파일명 변경

print("Plot saved")