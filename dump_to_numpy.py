import numpy as np
import sys
import os

def parse_dump(filename):

    with open(filename, 'r') as f:
        lines = f.readlines()

    # 1) TIMESTEP 건너뛰기
    i = 0
    if lines[i].startswith("ITEM: TIMESTEP"):
        i += 2

    # 2) NUMBER OF ATOMS
    assert lines[i].startswith("ITEM: NUMBER OF ATOMS")
    i += 1
    N = int(lines[i].strip())
    i += 1

    # 3) BOX BOUNDS
    assert lines[i].startswith("ITEM: BOX BOUNDS")
    box_info = lines[i].strip().split()[3:]
    i += 1
    xlo, xhi, _ = map(float, lines[i].split())
    i += 1
    ylo, yhi, _ = map(float, lines[i].split())
    i += 1
    zlo, zhi, _ = map(float, lines[i].split())
    i += 1

    # 4) ATOMS 헤더
    assert lines[i].startswith("ITEM: ATOMS")
    header = lines[i].split()[2:]
    i += 1
    is_scaled = False
    if any(col in ("xs","ys","zs") for col in header):
        is_scaled = True

    # 5) 원자 좌표 읽기
    positions = np.zeros((N,3), dtype=float)
    xcol = header.index("x") if "x" in header else (header.index("xs") if "xs" in header else -1)
    ycol = header.index("y") if "y" in header else (header.index("ys") if "ys" in header else -1)
    zcol = header.index("z") if "z" in header else (header.index("zs") if "zs" in header else -1)
    if xcol < 0 or ycol < 0 or zcol < 0:
        raise RuntimeError("ATOM 헤더에 x,y,z 혹은 xs,ys,zs 컬럼이 없습니다.")

    for idx in range(N):
        tokens = lines[i+idx].split()
        xs = float(tokens[xcol])
        ys = float(tokens[ycol])
        zs = float(tokens[zcol])
        if is_scaled:
            xs = xlo + xs*(xhi - xlo)
            ys = ylo + ys*(yhi - ylo)
            zs = zlo + zs*(zhi - zlo)
        positions[idx,0] = xs
        positions[idx,1] = ys
        positions[idx,2] = zs

    return positions, np.array([xlo,xhi,ylo,yhi,zlo,zhi])



def make_voxel_image(positions, box, voxel_size=0.1, atom_radius=0.34):

    xlo,xhi,ylo,yhi,zlo,zhi = box
    nx = int(np.ceil((xhi-xlo)/voxel_size))
    ny = int(np.ceil((yhi-ylo)/voxel_size))
    nz = int(np.ceil((zhi-zlo)/voxel_size))

    # 좌표별 격자 index 대응 함수
    def coord_to_index(x,y,z):
        ix = int((x - xlo) / voxel_size)
        iy = int((y - ylo) / voxel_size)
        iz = int((z - zlo) / voxel_size)
        return ix, iy, iz

    # 3D 볼륨 초기화: False=빈 공간, True=고체
    grid = np.zeros((nx, ny, nz), dtype=bool)

    # 원자 하나씩 복셀화: 반경(atom_radius) 이내는 모두 True로 설정
    r_pix = int(np.ceil(atom_radius / voxel_size))
    for (x,y,z) in positions:
        ix,iy,iz = coord_to_index(x,y,z)
        x0 = max(ix-r_pix, 0);   x1 = min(ix+r_pix+1, nx)
        y0 = max(iy-r_pix, 0);   y1 = min(iy+r_pix+1, ny)
        z0 = max(iz-r_pix, 0);   z1 = min(iz+r_pix+1, nz)
        for ixi in range(x0, x1):
            dx = (ixi + 0.5)*voxel_size + xlo - x  # 복셀 중심 좌표 vs 원자 중심
            dx2 = dx*dx
            for iyi in range(y0, y1):
                dy = (iyi + 0.5)*voxel_size + ylo - y
                dy2 = dy*dy
                for izi in range(z0, z1):
                    dz = (izi + 0.5)*voxel_size + zlo - z
                    if dx2 + dy2 + dz*dz <= atom_radius*atom_radius:
                        grid[ixi, iyi, izi] = True

    # “True” 영역(고체) 외에는 자동으로 “기공” 취급(False)
    return grid, (nx,ny,nz), voxel_size


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_image.py <dump 파일명> [voxel_size_nm]")
        sys.exit(1)

    dumpfile = sys.argv[1]
    vs = float(sys.argv[2]) if len(sys.argv)>2 else 0.1  # 기본 0.1 nm 복셀
    print("Dump 파일 읽는 중 →", dumpfile)
    pos, box = parse_dump(dumpfile)
    print("원자 개수:", pos.shape[0], "박스 크기(NM):", box)

    print("격자 생성 (voxel_size=%.3f nm) 중..." % vs)
    grid, dims, vsz = make_voxel_image(pos, box, voxel_size=vs, atom_radius=0.17)
    nx,ny,nz = dims
    print(f"격자 크기: {nx}×{ny}×{nz}")

    # numpy 배열 저장
    outname = os.path.splitext(os.path.basename(dumpfile))[0] + "_image.npy"
    print("결과 이미지를 저장:", outname)
    np.save(outname, grid)
    print("완료.")