from astropy.io import fits
import numpy as np
import scipy.special as scisp
import matplotlib.pyplot as plt
import sys, os

#input_file = '/home/bperri/Python_Pluto/pipelines/adapt/adapt40311_02a012_201703191600_i00005600n0.fts'
input_file = '/home/soumitra/Final_Study/Map_Preparation_adapt/adapt41311_03k012_201207122000_i00015200n1.fts'
#map_name = 'adapt_{}'.format(input_file.split('/')[-1][18:26])
map_name='adapt_20120712_1200f_lmax15'
map = 0
l_max = 15
nb_modes_tot = (l_max+1)*(l_max+2)//2 - 1 # computing total number of modes
check = 'no'

if sys.version_info[0] == 3:
    xrange = range

#################################################################
# Opening fits
#################################################################
input_data = fits.getdata(input_file, ext=0)
shape = np.shape(input_data)
nb_maps = shape[0]
nb_th = shape[1]
nb_phi = shape[2]
Br = input_data[map,:,:]
theta = np.linspace(0.,np.pi,nb_th)
phi = np.linspace(0.,2.0*np.pi,nb_phi)
Theta = np.tile(theta, (nb_phi,1)).T
Phi = np.tile(phi, (nb_th,1))

#################################################################
# Projection of Br
#################################################################
print('Beginning of projection')
dtheta = np.tile(np.concatenate([np.diff(theta),[theta[-1]-theta[-2]]]),(nb_phi,1)).T
dphi = np.tile(np.concatenate([np.diff(phi),[phi[1]-phi[0]]]),(nb_th,1))
coefbr = np.zeros(nb_modes_tot, dtype=np.complex)
ylm = np.zeros((nb_modes_tot, nb_th, nb_phi), dtype=np.complex)
mod = 0
for l in xrange(1, l_max+1):
  for m in xrange(0, l+1):
    #print('{} {} {}'.format(l, m, mod))
    ylm[mod] = scisp.sph_harm(m, l, Phi, Theta)
    ylm_c = np.conj(ylm[mod])
    integrand_a = Br*ylm_c
    integrand_a = integrand_a*np.sin(Theta)*dphi*dtheta
    coefbr[mod] = np.sum(integrand_a)
   # if (m!=0):
   #    coefbr[mod]= 0.0
    mod = mod+1
print('End of projection')

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
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, 0, 0)) # writing
    mod = mod+1
F.write('\n') # line skip
# Gamma
mod = 0
for l in xrange(1,l_max+1): # loop on l
  for m in xrange(0,l+1): # loop on m
    F.write('{} {} {: .6e} {: .6e} \n'.format(l, m, 0, 0)) # writing
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
  fig = plt.figure()
  plt.plot(theta, Br[:,0], marker='o', label='adapt')
  plt.plot(theta, Br_reconst[:,0], marker='o', label='pluto')
  plt.legend()
  fig.show()
