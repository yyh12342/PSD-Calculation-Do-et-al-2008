import numpy as np
import porespy as ps
import openpnm as op
import matplotlib.pyplot as plt

# 1) 3D boolean 이미지 로드 (True=기공, False=고체)
grid = np.load('dump_image.npy')
pore_space = ~grid

# 2) 해상도 및 다운샘플
voxel_size = 0.2  # [nm]

downsample_factor = 2
if downsample_factor > 1:
    sx, sy, sz = pore_space.shape
    nx, ny, nz = sx // downsample_factor, sy // downsample_factor, sz // downsample_factor
    reduced = np.zeros((nx, ny, nz), dtype=bool)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                block = pore_space[
                    i*downsample_factor:(i+1)*downsample_factor,
                    j*downsample_factor:(j+1)*downsample_factor,
                    k*downsample_factor:(k+1)*downsample_factor
                ]
                reduced[i,j,k] = np.any(block)
    pore_space = reduced

print("다운샘플된 pore_space shape:", pore_space.shape)
print("voxel_size (nm):", voxel_size)

# 3) snow2() 호출 시 chunk_size, num_workers를 개별 인자로 넘김
snow_output = ps.networks.snow2(
    phases=pore_space,            # True=기공, False=고체
    voxel_size=voxel_size,        # 한 voxel 크기 [nm]
    r_max=5,                      # peak 검출 시 사용할 최대 반경 (픽셀 단위)
    chunk_size=[128, 128, 128],   # Dask 청크 크기
    num_workers=1                 # 워커 개수 1로 제한
)

# 4) OpenPNM 네트워크로 변환
net = op.io.network_from_porespy(snow_output.network)

# 5) pore.diameter 추출 후 히스토그램 그리기
diameters = net['pore.diameter']   # 단위: nm (voxel_size 기준)
plt.figure(figsize=(6,4))
plt.hist(diameters, bins=30, edgecolor='black', alpha=0.7)
plt.xlabel('Pore Diameter [nm]')
plt.ylabel('Number of Pores')
plt.title('Pore Size Distribution (SNOW2, 옵션 수정)')
plt.tight_layout()
plt.show()


# import numpy as np
# import porespy as ps
# import matplotlib.pyplot as plt

# # True = pore, False = solid
# pore_space = ~np.load('dump_image.npy')
# voxel_size = 0.1  # 또는 0.2, 다운샘플에 맞춰서

# bins, counts = ps.metrics.pore_size_distribution(
#     pore_space,
#     voxel_size=voxel_size
# )
# plt.bar(bins, counts, width=bins[1]-bins[0], edgecolor='k')
# plt.xlabel('Pore Diameter [nm]')
# plt.ylabel('Frequency')
# plt.title('PSD via pore_size_distribution')
# plt.show()