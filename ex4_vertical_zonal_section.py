import numpy as np
import matplotlib.pyplot as plt

# set parameters (ex4_vertical_section.py と同じ構成)

expid = 'ex4-case1'

nbgn = 60
nend = 60
nskp = 1

im = 100
jm = 40
km = 12
dx = 100.e3
dy = 100.e3
dz = 50.e0

# set coordinates

xc = (np.arange(im+2) - 0.5) * dx / 1000.  # [km]
yc = (np.arange(jm+2) - 0.5) * dy / 1000.  # [km]
zc = (np.arange(km+2) - 0.5 - km) * dz     # [m] (負が深さ)

XC, YC = np.meshgrid(xc, yc, indexing='ij')

dtyp = np.dtype([
    ('time','<d'),
    ('u','<' + str((im+2)*(jm+2)*(km+2)) + 'd'),
    ('v','<' + str((im+2)*(jm+2)*(km+2)) + 'd'),
    ('e','<' + str((im+2)*(jm+2))       + 'd'),
    ('w','<' + str((im+2)*(jm+2)*(km+2)) + 'd'),
    ('p','<' + str((im+2)*(jm+2)*(km+2)) + 'd'),
    ('t','<' + str((im+2)*(jm+2)*(km+2)) + 'd'),
    ('s','<' + str((im+2)*(jm+2)*(km+2)) + 'd'),
    ('r','<' + str((im+2)*(jm+2)*(km+2)) + 'd'),
])

# ユーザー入力：東西-鉛直断面の j（緯度方向index）
while True:
    try:
        j_zonal = int(input(f'東西-鉛直断面を見る j のindexを入力してください（1~{jm}）: '))
        if 1 <= j_zonal <= jm:
            break
        else:
            print(f'1から{jm}の範囲で入力してください。')
    except ValueError:
        print('整数を入力してください。')

# 見やすさのため：w は桁が小さいことが多いので矢印だけ拡大（物理値は変えない）
w_arrow_scale = 7000.0  # 必要なら 50, 1000 などに調整

for n in range(nbgn, nend+1, nskp):

    # read data
    fname = expid + f'.n{n:06}'
    with open(fname, 'rb') as fp:
        chunk = np.fromfile(fp, dtype=dtyp)

    time = chunk['time'][0]
    u = chunk['u'][0].reshape((im+2, jm+2, km+2), order="F")
    v = chunk['v'][0].reshape((im+2, jm+2, km+2), order="F")
    e = chunk['e'][0].reshape((im+2, jm+2),       order="F")
    w = chunk['w'][0].reshape((im+2, jm+2, km+2), order="F")
    p = chunk['p'][0].reshape((im+2, jm+2, km+2), order="F")
    t = chunk['t'][0].reshape((im+2, jm+2, km+2), order="F")
    s = chunk['s'][0].reshape((im+2, jm+2, km+2), order="F")
    r = chunk['r'][0].reshape((im+2, jm+2, km+2), order="F")

    # shift staggered variables to the cell center（ex4_vertical_section.py と同じ）
    uc = np.zeros([im+2, jm+2, km+2])
    vc = np.zeros([im+2, jm+2, km+2])
    wc = np.zeros([im+2, jm+2, km+2])

    uc[1:im, :, :] = 0.5 * (u[0:im-1, :, :] + u[1:im, :, :])
    uc[0, :, :]    = -uc[1, :, :]
    uc[im+1, :, :] = -uc[im, :, :]

    vc[:, 1:jm, :] = 0.5 * (v[:, 0:jm-1, :] + v[:, 1:jm, :])
    vc[:, 0, :]    = -vc[:, 1, :]
    vc[:, jm+1, :] = -vc[:, jm, :]

    wc[:, :, 1:km] = 0.5 * (w[:, :, 0:km-1] + w[:, :, 1:km])
    wc[:, :, 0]    = -wc[:, :, 1]
    wc[:, :, km+1] = -wc[:, :, km]

    # ===== 東西-鉛直断面（j固定）：(x, z) =====
    j_idx = j_zonal

    # 断面データ（内部格子のみ）
    u_sec = uc[1:im+1, j_idx, 1:km+1]   # (im, km)
    w_sec = wc[1:im+1, j_idx, 1:km+1]
    t_sec =  t[1:im+1, j_idx, 1:km+1]
    r_sec =  r[1:im+1, j_idx, 1:km+1]

    # メッシュ
    Xc_sec = xc[1:im+1]      # [km]
    Zc_sec = zc[1:km+1]      # [m]
    X_mesh, Z_mesh = np.meshgrid(Xc_sec, Zc_sec, indexing='ij')  # (im, km)

    # 速度の大きさ（色付け用：物理値のまま）
    uw_abs = np.sqrt(u_sec**2 + w_sec**2)

    # quiver を間引き（必要なら調整）
    sx = 2  # x方向 2点に1本
    sz = 1  # z方向は全て
    Xq = X_mesh[::sx, ::sz]
    Zq = Z_mesh[::sx, ::sz]
    Uq = u_sec[::sx, ::sz]
    Wq = (w_sec * w_arrow_scale)[::sx, ::sz]  # 矢印だけ拡大

    # ===== plot =====
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    # 1) (u, w) ベクトル（ex4.py風：矢印を色で流速表示）
    ax = axes[0]
    Q = ax.quiver(
        Xq, Zq, Uq, Wq, uw_abs[::sx, ::sz],
        cmap='jet', pivot='mid'
    )
    plt.colorbar(Q, ax=ax, label='Speed sqrt(u^2+w^2) [m/s]', shrink=0.8)
    ax.set_ylabel('Depth [m]')
    ax.set_title(f'Velocity (u,w) on Zonal-Vertical Section (y={yc[j_idx]:.0f} km, j={j_idx})'
                 f'\n(w arrow scaled ×{w_arrow_scale:g})')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([zc[1], zc[km]])

    # 2) 温度
    ax = axes[1]
    cf = ax.contourf(X_mesh, Z_mesh, t_sec, levels=20, cmap='jet')
    ax.contour(X_mesh, Z_mesh, t_sec, levels=10, colors='k', linewidths=0.5, alpha=0.5)
    plt.colorbar(cf, ax=ax, label='Temperature [°C]')
    ax.set_ylabel('Depth [m]')
    ax.set_title(f'Temperature on Zonal-Vertical Section (y={yc[j_idx]:.0f} km, j={j_idx})')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([zc[1], zc[km]])

    # 3) 密度
    ax = axes[2]
    cf = ax.contourf(X_mesh, Z_mesh, r_sec, levels=20, cmap='viridis')
    ax.contour(X_mesh, Z_mesh, r_sec, levels=10, colors='k', linewidths=0.5, alpha=0.5)
    plt.colorbar(cf, ax=ax, label='Density [kg/m³]')
    ax.set_xlabel('X [km]')
    ax.set_ylabel('Depth [m]')
    ax.set_title(f'Density on Zonal-Vertical Section (y={yc[j_idx]:.0f} km, j={j_idx})')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([zc[1], zc[km]])

    fig.tight_layout()
    plt.show()

# 最後の図を保存（ex4_vertical_section.py に合わせた書き方）
fig.savefig(f'ex4_zonal_vertical_j{j_zonal:03d}.png')
plt.close()
