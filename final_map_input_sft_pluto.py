# %load final_map_input.py
import numpy as np
from astropy.io import fits
import scipy.special as scisp
from scipy import interpolate
import matplotlib.pyplot as plt
import sys, os
import glob
from scipy.io import readsav

dirname= '/home/soumitra/Solar_Wind_Study_France/zdimap'
input_file= dirname + '/synoptic_br_sft_2021.12.04_00_00_00.fits'
map_name='SFT_20211204D'
maps_name='SFT_20211204D_modified'
mapd_name='SFT_20211204D_rescaled'
map = 0
l_max = 40
nb_modes_tot = (l_max+1)*(l_max+2)//2 - 1 # computing total number of modes
check = 'yes'
plot='yes'

if sys.version_info[0] == 3:
    xrange = range


input_data = fits.getdata(input_file, ext=0)
input_data=np.nan_to_num(input_data)
wilcox2gauss=0.01
sft2gauss=1.0
input_data=sft2gauss*input_data
shape = np.shape(input_data)
nb_th = shape[0]
nb_phi = shape[1]
Br11 = input_data
theta = np.linspace(0.,np.pi,nb_th)
phi = np.linspace(0.,2.0*np.pi,nb_phi)
Theta = np.tile(theta, (nb_phi,1)).T
Phi = np.tile(phi, (nb_th,1))
bmax = 30
#plt.figure(figsize=[11,5])
#ax1 = plt.subplot(111)
#ax1.pcolormesh(np.rad2deg(phi),np.rad2deg(theta),Br[::-1,:],cmap='bwr',shading='auto',vmax=bmax,vmin=-bmax)
#ax1.set_title('SFT (2017.08.21)')
#plt.show()

Br11 = Br11[:,:]

# shifting to the right in phi by 180 degrees
br1=np.zeros((nb_th,nb_phi))
rot=nb_th #501
for i in range(nb_th):
	for j in range(nb_phi-rot):
			j1=j + rot
			br1[i,j1]=Br11[i,j]

for i in range(nb_th):
		for j in range(nb_phi-rot,nb_phi):
			j1=j-rot+1
			br1[i,j1]=Br11[i,j]
Br = br1

#################################################################
# Projection of Br
#################################################################
print('Beginning of projection')
dtheta = np.tile(np.concatenate([np.diff(theta),[theta[-1]-theta[-2]]]),(nb_phi,1)).T
dphi = np.tile(np.concatenate([np.diff(phi),[phi[1]-phi[0]]]),(nb_th,1))
coefbr = np.zeros(nb_modes_tot, dtype=complex)
ylm = np.zeros((nb_modes_tot, nb_th, nb_phi), dtype=complex)
coefs = np.zeros((l_max+1, l_max+2), dtype=complex)
coefg = np.zeros((l_max+1, l_max+2), dtype=complex)
mod = 0
sumg=0
for l in xrange(1, l_max+1):
  for m in xrange(0, l+1):
    print('{} {} {}'.format(l, m, mod))
    ylm[mod] = scisp.sph_harm(m, l, Phi, Theta)
    ylm_c = np.conj(ylm[mod])
    integrand_a = Br*ylm_c
    integrand_a = integrand_a*np.sin(Theta)*dphi*dtheta 
    coefbr[mod] = 1.0*np.sum(integrand_a)
    coefs[l,m]= 1.0*np.sum(integrand_a)
    sumg= sumg+ (np.real(coefs[l,m]))**2 + (np.imag(coefs[l,m]))**2
# Scaled by 2.0
    coefg[l,m]= 2.0*np.sum(integrand_a)
    #if (m!=0):
    #   coefbr[mod]= 0.0
    mod = mod+1
  #print('%2d%2d'%l,m)
print('End of projection')


coefd = np.zeros((l_max+1, l_max+2), dtype=complex)
for l in xrange(1, l_max+1):
  sums=0.0
  for m in xrange(0, l+1):
    sums=sums + (np.real(coefs[l,m]))**2 + (np.imag(coefs[l,m]))**2
    coefd[l,0]= np.sign(np.real(coefs[l,0]))*np.sqrt(sums)
    if (m!=0):
         coefd[l,m]= 0.0



coeft = np.zeros(l_max+1)
coefm=np.zeros(l_max+1)
coefp=np.zeros(l_max+1)
modt = np.zeros(l_max+1)
for l in xrange(1, l_max+1):
    modt[l]=l
    coeft[l]= (np.real(coefd[l,0]))**2 + (np.imag(coefd[l,0]))**2
    coefp[l]= (np.real(coefs[l,0]))**2 + (np.imag(coefs[l,0]))**2
    coefm[l]= (np.real(coefg[l,0]))**2 + (np.imag(coefg[l,0]))**2


#################################################################
# Writing map
#################################################################
print('Writing coefficients')
F = open(map_name,"w") # opening file to write coefficients
F.write('ZDI map \n') # writing a line of text
F.write('{} 3 -3\n'.format(nb_modes_tot)) # writing total number of modes in the file
# Alpha
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, np.real(coefbr[mod]), np.imag(coefbr[mod]))) # writing
    mod = mod+1
F.write('\n') # line skip
# Beta
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, np.real(coefbr[mod]), np.imag(coefbr[mod]))) # writing
    mod = mod+1
F.write('\n') # line skip
# Gamma
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, np.real(coefbr[mod]), np.imag(coefbr[mod]))) # writing
    mod = mod+1
F.write('\n') # line skip
F.close() # closing file
print('File closed. \n')



#################################################################
# Writing maps
#################################################################
print('Writing coefficients')
F = open(maps_name,"w") # opening file to write coefficients
F.write('ZDIS map \n') # writing a line of text
F.write('{} 3 -3\n'.format(nb_modes_tot)) # writing total number of modes in the file
# Alpha
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, np.real(coefd[l,m]), np.imag(coefd[l,m]))) # writing
    mod = mod+1
F.write('\n') # line skip
# Beta
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, np.real(coefd[l,m]), np.imag(coefd[l, m]))) # writing
    mod = mod+1
F.write('\n') # line skip
# Gamma
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, np.real(coefd[l, m]), np.imag(coefd[l, m]))) # writing
    mod = mod+1
F.write('\n') # line skip
F.close() # closing file
print('File closed. \n')

#################################################################
# Writing maps
#################################################################
print('Writing coefficients')
F = open(mapd_name,"w") # opening file to write coefficients
F.write('ZDI map \n') # writing a line of text
F.write('{} 3 -3\n'.format(nb_modes_tot)) # writing total number of modes in the file
# Alpha
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, (np.real(coefd[l,m])+np.imag(coefd[l,m])), 0.0)) # writing
    mod = mod+1
F.write('\n') # line skip
# Beta
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, (np.real(coefd[l,m])+np.imag(coefd[l,m])), 0.0)) # writing
    mod = mod+1
F.write('\n') # line skip
# Gamma
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, (np.real(coefd[l,m])+np.imag(coefd[l,m])), 0.0)) # writing
    mod = mod+1
F.write('\n') # line skip
F.close() # closing file
print('File closed. \n')

#################################################################
# Check reconstruction
#################################################################
if (check == 'yes'):
  Br_reconst = np.zeros((nb_th, nb_phi)) 
  mod = 0
  for l in xrange(1, l_max+1):
    for m in xrange(0, l+1):
      Br_reconst = Br_reconst + np.real(coefbr[mod]*ylm[mod])
      mod = mod+1
  fig = plt.figure(figsize=[11,9])
  ax1 = plt.subplot(211)
  ax1.pcolormesh(np.rad2deg(phi),np.rad2deg(theta),Br,cmap='bwr',shading='auto',vmax=bmax,vmin=-bmax)
  ax1.set_title('SFT (20211204)')
  ax2 = plt.subplot(212)
  ax2.pcolormesh(np.rad2deg(phi),np.rad2deg(theta),Br_reconst,cmap='bwr',shading='auto',vmax=bmax,vmin=-bmax)
  ax2.set_title('Reconstructed (20211204)')
  #plt.plot(theta, Br[:,0], marker='o', label='SFT')
  #plt.plot(theta, Br_reconst[:,0], marker='o', label='pluto')
  #plt.legend()
  plt.subplots_adjust(hspace=0.15)
  fig.savefig('pfss_reconstructed_final_new.png',dpi=200)

if (plot=='yes'):
   fig=plt.figure()
   plt.plot(modt, coefp, '-k', label='axisymmetric energy')
   plt.plot(modt, coeft, '-b', label='modified axisymmetric energy')
   plt.plot(modt, coefm, '-r', label='scaled axisymmetric energy')
   plt.ylabel('alpha^2_l0')
   plt.xlabel('L')
   plt.legend()
   fig.show()





