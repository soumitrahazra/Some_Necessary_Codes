import os
import numpy as np
import matplotlib.pyplot as plt
import ZDIUtils as zdi
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
import sys

if sys.version_info[0] == 3:
    xrange = range

def createZDImap(alm,mapname,path="./"):
    zdimap=path+mapname
    lmax = int((-3+np.sqrt(9+8*len(str(alm))))/2.)
    if not os.path.exists(zdimap):
        f=open(zdimap,'w')
        f.write("Br spherical harmonices coeffs\n")
        f.write("{:d} {:d} {:d}\n".format(len(alm),0,0))
        k=0
        for i in xrange(1,lmax+1):
            for j in range(i+1):
                f.write("{:1d} {:1d} {:+.05e} {:+.05e}\n".format(i,j,alm[k].real,alm[k].imag))
                k=k+1
        f.write("\n")
        for i in xrange(1,lmax+1):
            for j in range(i+1):
                f.write("{:d} {:d} {:+.05e} {:+.05e}\n".format(i,j,0.,0.))
        f.write("\n")
        for i in range(1,lmax+1):
            for j in range(i+1):
                f.write("{:d} {:d} {:+.05e} {:+.05e}\n".format(i,j,0.,0.))
        f.write("\n")
        f.close()
    
#nso_fits=fits.open('./adapt40311_02a012_201703191600_i00005600n0.fts')
nso_fits=fits.open('adapt40311_02a012_201703191600_i00005600n0.fts')
#hdul=fits.open(nso_fits)

mag=nso_fits[0].data

lmax=1
# Projection on lmax harmonics
ngrid=(mag[0].shape)[0]
theta,phi=np.mgrid[0:np.pi:(ngrid)*1j, 0:2*np.pi:(2*ngrid)*1j]
spharm=zdi.mysph(lmax,theta,phi)
#Mlm,alm,all=spharm.spherical_harmonics_decomposition(mag[0][::-1,0],theta[:,0],lmax=lmax,sym="axis")
Mlm,alm,all=spharm.spherical_harmonics_decomposition(mag[0][::-1,:],theta,lmax=lmax)

#createZDImap(alm,"adapt201703191600_lmax{}".format(lmax))
createZDImap(alm,"adapt20100610_lmax{}".format(lmax))
br,bt,bp=spharm.multipolarExpansion(alm,rb=1.0,rsph=1.0)

print(zdi.cmpMagFlux(theta,phi,mag[0],1,abs=False),zdi.cmpMagFlux(theta,phi,br,1,abs=False))

fig,(ax1,ax2)=plt.subplots(2,1)

cmap="binary"
im1=ax1.pcolormesh(phi[0,:],theta[:,0],mag[0],cmap=cmap)
divider=make_axes_locatable(ax1)
cax=divider.append_axes("right",size="10%",pad=0.1)
cbar=plt.colorbar(im1,cax=cax)

im2=ax2.pcolormesh(phi[0,:],theta[:,0],br,cmap=cmap)
divider=make_axes_locatable(ax2)
cax=divider.append_axes("right",size="10%",pad=0.1)
cbar=plt.colorbar(im2,cax=cax)


plt.show()


