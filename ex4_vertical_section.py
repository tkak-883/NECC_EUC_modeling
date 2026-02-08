import numpy as np
import matplotlib.pyplot as plt

# set parameters

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

xc = (np.arange(im+2)-0.5 ) * dx / 1000.
yc = (np.arange(jm+2)-0.5 ) * dy / 1000.
zc = (np.arange(km+2)-0.5 - km) * dz 

print(zc)

XC,YC = np.meshgrid(xc, yc, indexing = 'ij') 

dtyp = np.dtype( [('time','<d'),\
                  ('u','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('v','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('e','<'+str((im+2)*(jm+2))+'d'),\
                  ('w','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('p','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('t','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('s','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('r','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ] )

# ユーザー入力：子午面のiインデックス
while True:
    try:
        i_meridional = int(input(f'子午面を見るiのindexを入力してください（1~{im}）: '))
        if 1 <= i_meridional <= im:
            break
        else:
            print(f'1から{im}の範囲で入力してください。')
    except ValueError:
        print('整数を入力してください。')

for n in range(nbgn,nend+1,nskp):

    # read data
    
    fname = expid + f'.n{n:06}'
    fp=open(fname,'rb')
    chunk = np.fromfile(fp, dtype=dtyp)
    fp.close()

    time = chunk['time'][0]
    u = chunk['u'][0].reshape((im+2,jm+2,km+2),order="F")
    v = chunk['v'][0].reshape((im+2,jm+2,km+2),order="F")
    e = chunk['e'][0].reshape((im+2,jm+2),order="F")
    w = chunk['w'][0].reshape((im+2,jm+2,km+2),order="F")
    p = chunk['p'][0].reshape((im+2,jm+2,km+2),order="F")
    t = chunk['t'][0].reshape((im+2,jm+2,km+2),order="F")
    s = chunk['s'][0].reshape((im+2,jm+2,km+2),order="F")
    r = chunk['r'][0].reshape((im+2,jm+2,km+2),order="F")

    # shift staggerd variables to the cell center
    
    uc=np.zeros([im+2,jm+2,km+2])
    vc=np.zeros([im+2,jm+2,km+2])
    wc=np.zeros([im+2,jm+2,km+2])
    uc[1:im,:,:]=0.5*(u[0:im-1,:,:]+u[1:im,:,:])
    uc[0,:,:]=-uc[1,:,:]
    uc[im+1,:,:]=-uc[im,:,:]
    vc[:,1:jm,:]=0.5*(v[:,0:jm-1,:]+v[:,1:jm,:])
    vc[:,0,:]=-vc[:,1,:]
    vc[:,jm+1,:]=-vc[:,jm,:]
    wc[:,:,1:km]=0.5*(w[:,:,0:km-1]+w[:,:,1:km])
    wc[:,:,0]=-wc[:,:,1]
    wc[:,:,km+1]=-wc[:,:,km]

    ### Vertical section at the meridional plane
    # 子午面（特定の経度）での断面
    i_idx = i_meridional  # ユーザーが指定したインデックス
    
    # 子午面での東西流速の鉛直断面 (y, z)
    u_section = uc[i_idx, 1:jm+1, 1:km+1]
    
    # 子午面での温度の鉛直断面
    t_section = t[i_idx, 1:jm+1, 1:km+1]
    
    # 子午面での密度の鉛直断面
    r_section = r[i_idx, 1:jm+1, 1:km+1]

    # メッシュグリッドの作成（Y軸とZ軸）
    YC_section = yc[1:jm+1]
    ZC_section = zc[1:km+1]
    Y_mesh, Z_mesh = np.meshgrid(YC_section, ZC_section, indexing='ij')
    
    # プロット作成
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # 1. 東西流速の鉛直断面
    ax = axes[0]
    levels = np.linspace(-0.2, 0.4, 21)
    cf = ax.contourf(Y_mesh, Z_mesh, u_section, levels=levels, cmap='RdBu_r', extend='both')
    ax.contour(Y_mesh, Z_mesh, u_section, levels=[0], colors='k', linewidths=2)
    plt.colorbar(cf, ax=ax, label='Zonal Velocity U [m/s]')
    ax.set_ylabel('Depth [m]')
    ax.set_title(f'Zonal Velocity at Meridional Plane (x={xc[i_idx]:.0f} km, i={i_idx})')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([zc[1], zc[km]]) 
    
    # 2. 温度の鉛直断面
    ax = axes[1]
    cf = ax.contourf(Y_mesh, Z_mesh, t_section, levels=20, cmap='jet')
    ax.contour(Y_mesh, Z_mesh, t_section, levels=10, colors='k', linewidths=0.5, alpha=0.5)
    plt.colorbar(cf, ax=ax, label='Temperature [°C]')
    ax.set_ylabel('Depth [m]')
    ax.set_title(f'Temperature at Meridional Plane (x={xc[i_idx]:.0f} km, i={i_idx})')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([zc[1], zc[km]]) 
    
    # 3. 密度の鉛直断面
    ax = axes[2]
    cf = ax.contourf(Y_mesh, Z_mesh, r_section, levels=20, cmap='viridis')
    ax.contour(Y_mesh, Z_mesh, r_section, levels=10, colors='k', linewidths=0.5, alpha=0.5)
    plt.colorbar(cf, ax=ax, label='Density [kg/m³]')
    ax.set_xlabel('Y [km]')
    ax.set_ylabel('Depth [m]')
    ax.set_title(f'Density at Meridional Plane (x={xc[i_idx]:.0f} km, i={i_idx})')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([zc[1], zc[km]])
    
    fig.tight_layout()
    plt.show()

fig.savefig(f'ex4_vertical_5year_para_i{i_meridional:03d}.png')
plt.close()