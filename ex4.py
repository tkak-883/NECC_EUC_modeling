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

    ### Velocity: ベクトルの長さは同じにして、色で流速を示す
    uc = uc[:,:,9]
    vc = vc[:,:,9]
    u_abs=np.sqrt( pow(uc,2) + pow(vc,2) )
    #uc = uc / u_abs
    #vc = vc / u_abs

    fig, ax = plt.subplots(figsize=(8,6))

    ax.contour(XC,YC,e[:,:])
    ax.grid()
    ax.set_title('Sea Surface Height')
    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')

    Q = ax.quiver(XC[1:-1,1:-1], YC[1:-1,1:-1], \
                  uc[1:-1,1:-1], vc[1:-1,1:-1], u_abs[1:-1,1:-1], cmap='jet', pivot='mid' )
    plt.colorbar(Q, label='Velocity [m/s]', shrink=0.6)

    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    plt.show()

fig.savefig('report_5year_k_9_para.png')
plt.close()