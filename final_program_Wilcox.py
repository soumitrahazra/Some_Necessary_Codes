import numpy as np
from astropy.io import fits
import scipy.special as scisp
import matplotlib.pyplot as plt
import sys, os
import glob
from scipy.io import readsav

# -*- coding: utf-8 -*-

def getListOfFiles(dirName):
    # create a list of file and sub directories
    # names in the given directory
    listOfFile = sorted(os.listdir(dirName))
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(fullPath):
            allFiles = allFiles + sorted(glob.glob(fullPath+ '/*.fits'))
        else:listoffiles=sorted(glob.glob(dirname+'/*.F.fits'))
            allFiles.append(fullPath)
    
    return allFiles

#dirname= '/home/soumitra/Final_Study/data'
dirname= '/home/soumitra/Final_Study/harmonic-analysis/Wilcox_Magnetogram'
path= '/home/soumitra/Final_Study/Figure/'
#listoffiles=sorted(getListOfFiles(dirname))
#listoffiles=sorted(glob.glob(dirname+'/*/*.fits'))

#listoffiles= sorted(glob.glob("/home/soumitra/Final_Study/data/data2/*.fits"))
ns=len(listoffiles)-1
print(ns)
l_max = 10
nb_modes_tot = (l_max+1)*(l_max+2)//2 - 1 # computing total number of modes

if sys.version_info[0] == 3:
    xrange = range

# defining the necessary arrays
axial_dipole_coef= np.zeros(ns)
equatorial_dipole_coef= np.zeros(ns)
equatorial_axial = np.zeros(ns)
coeff20= np.zeros(ns)
coeff21= np.zeros(ns)
coeff22= np.zeros(ns)
coeff30= np.zeros(ns)
coeff31= np.zeros(ns)
coeff32= np.zeros(ns)
coeff33= np.zeros(ns)
dipole_energy= np.zeros(ns)
quadrapolar_energy= np.zeros(ns)
quad_dip_ratio= np.zeros(ns)
time = np.zeros(ns)

#Final loop for spherical harmonic expansion and corresponding calculation
for i in xrange(0, ns):
    input_file= listoffiles[i]
    # opening the data fits file
    #fits_image_filename = fits.util.get_testdata_filepath(input_file)
    #hdulist= fits.open(input_file, ignore_missing_end=True)
    #input_data=hdulist[0].data
    #hdr=hdulist[0].header
    #hdulist.close()
    input_data = fits.getdata(input_file, ext=0)
    wilcox2gauss=0.01
    input_data=wilcox2gauss*input_data

    shape=np.shape(input_data)
    nb_th = shape[0]
    nb_phi = shape[1]
    Br = input_data
    Br= np.nan_to_num(Br)
    theta = np.linspace(0.,np.pi,nb_th)
    phi = np.linspace(0.,2.0*np.pi,nb_phi)
    Theta = np.tile(theta, (nb_phi,1)).T
    Phi = np.tile(phi, (nb_th,1))
    # Projection of Br
    dtheta = np.tile(np.concatenate([np.diff(theta),[theta[-1]-theta[-2]]]),(nb_phi,1)).T
    dphi = np.tile(np.concatenate([np.diff(phi),[phi[1]-phi[0]]]),(nb_th,1))
    coefbr = np.zeros(nb_modes_tot, dtype=np.complex)
    ylm = np.zeros((nb_modes_tot, nb_th, nb_phi), dtype=np.complex)
    subb= np.zeros(((l_max +1)*(l_max +1), (l_max +1)*(l_max +1)), dtype=np.complex)
    mod = 0
    for l in xrange(1, l_max+1):
        for m in xrange(0, l+1):
          #print('{} {} {}'.format(l, m, mod))
          ylm[mod] = scisp.sph_harm(m, l, Phi, Theta)
          ylm_c = np.conj(ylm[mod])
          integrand_a = Br*ylm_c
          integrand_a = integrand_a*np.sin(Theta)*dphi*dtheta
          coefbr[mod] = np.sum(integrand_a)
          subb[l, m] = coefbr[mod]
          mod = mod+1

    axial_dipole_coef[i]= np.real(subb[1,0])
    equatorial_dipole_coef[i] = np.sqrt(( np.real(subb[1,1])*np.real(subb[1,1])) + (np.imag(subb[1,1])*np.imag(subb[1,1])))
    equatorial_axial[i]= equatorial_dipole_coef[i]/ axial_dipole_coef[i]

    coeff20[i] = np.real(subb[2,0])
    coeff21[i]= np.real(subb[2,1])
    coeff22[i]= np.real(subb[2,2])

    coeff30[i]= np.real(subb[3,0])
    coeff31[i] = np.real(subb[3,1])
    coeff32[i]= np.real(subb[3,2])
    coeff33[i]= np.real(subb[3,3])

    bxd= np.real(subb[1,1])
    byd= np.imag(subb[1,1])
    bzd= np.real(subb[1,0])
    dipole_energy[i]= (bxd*bxd + byd*byd + bzd*bzd)

    bx1q= np.real(subb[2,0])
    bx2q= np.real(subb[2,1])
    bx3q= np.imag(subb[2,1])
    bx4q= np.real(subb[2,2])
    bx5q= np.imag(subb[2,2])
    quadrapolar_energy[i] = (bx1q*bx1q + bx2q*bx2q + bx3q*bx3q + bx4q*bx4q + bx5q*bx5q)

    quad_dip_ratio[i] = quadrapolar_energy[i]/dipole_energy[i]
    time[i] = 1642.0 +i


t2s=[1640, 1687, 1754, 1821, 1888, 1955, 2022, 2090, 2157]
l2s=[1977, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015]

xtimes=np.loadtxt("sunspot_data.txt")[:,0]
yss= np.loadtxt("sunspot_data.txt")[:,1]
#print(axial_dipole_coef)
# figure 1
#plt.subplot(3,1,1)
#plt.plot(time, axial_dipole_coef)
#plt.subplot(3,1,2)
#plt.plot(time, equatorial_dipole_coef)
#plt.subplot(3,1,3)
#plt.plot(time, equatorial_axial)
#plt.show()

s= readsav('harmonic_power_wilcox.sav')


#Figure 1
f, axarr = plt.subplots(3, sharex=True)
axarr[0].plot(time, abs(axial_dipole_coef))
axarr[0].set_ylabel('Axial Dipole Coefficient (G)')
axarr[1].semilogy(time, abs(equatorial_dipole_coef))
axarr[1].set_ylabel('Equatorial Dipole Magnitude (G)')
axarr[2].semilogy(time, abs(equatorial_axial))
axarr[2].set_ylabel('Equatorial to Axial Ratio')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[2].set_xticks(t2s)
axarr[2].set_xticklabels(l2s)
axarr[2].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw1.png',dpi=500)
plt.show()

#Figure 2
f, axarr = plt.subplots(3, sharex=True)
axarr[0].semilogy(time, abs(dipole_energy))
axarr[0].set_ylabel('Energy in Dipole (erg)')
axarr[1].semilogy(time, abs(quadrapolar_energy))
axarr[1].set_ylabel('Energy in Quadrapole (erg)')
axarr[2].semilogy(time, abs(quad_dip_ratio))
axarr[2].set_ylabel('Quadrapolar vs Dipolar ratio')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[2].set_xticks(t2s)
axarr[2].set_xticklabels(l2s)
axarr[2].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw2.png',dpi=500)
plt.show()

#Figure 3
f, axarr = plt.subplots(3, sharex=True)
axarr[0].semilogy(time, abs(coeff20))
axarr[0].set_ylabel('L=2 m=|0| magnitude (G)')
axarr[1].semilogy(time, abs(coeff21))
axarr[1].set_ylabel('L=2, m=|1| magnitude (G)')
axarr[2].semilogy(time, abs(coeff22))
axarr[2].set_ylabel('L=2, m=|2| magnitude (G)')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[2].set_xticks(t2s)
axarr[2].set_xticklabels(l2s)
axarr[2].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw3.png',dpi=500)
plt.show()


#Figure 4
f, axarr = plt.subplots(4, sharex=True)
axarr[0].semilogy(time, abs(coeff30))
axarr[0].set_ylabel('L=3 m=|0| magnitude (G)')
axarr[1].semilogy(time, abs(coeff31))
axarr[1].set_ylabel('L=3, m=|1| magnitude (G)')
axarr[2].semilogy(time, abs(coeff32))
axarr[2].set_ylabel('L=3, m=|2| magnitude (G)')
axarr[3].semilogy(time, abs(coeff33))
axarr[3].set_ylabel('L=3, m=|3| magnitude (G)')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[3].set_xticks(t2s)
axarr[3].set_xticklabels(l2s)
axarr[3].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw4.png',dpi=500)
plt.show()

#Figure 5
f, axarr = plt.subplots(3, sharex=True)
axarr[0].plot(time, -axial_dipole_coef)
axarr[0].plot(time, -np.real(s.brcoeff[:,0,1]), '-r')
axarr[0].legend(('Our Python code', 'IDL code (Mark)'))
axarr[0].axhline(y=0)
axarr[0].set_ylabel('Axial Dipole Coefficient (G)')
axarr[1].plot(time, equatorial_dipole_coef)
axarr[1].plot(time, np.sqrt(s.br2energy[:,1,1]), '-r')
axarr[1].legend(('Our Python code','IDL code (Mark)'))
axarr[1].set_ylabel('Equatorial Dipole Magnitude (G)')
axarr[2].semilogy(time, abs(equatorial_axial))
axarr[2].semilogy(time, np.sqrt(s.br2energy[:,1,1]/s.br2energy[:,0,1]), '-r')
axarr[2].legend(('Our Python code', 'IDL code (Mark)'))
axarr[2].set_ylabel('Equatorial to Axial Ratio')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[2].set_xticks(t2s)
axarr[2].set_xticklabels(l2s)
axarr[2].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw5.png',dpi=500)
plt.show()

#Figure 6
f, axarr = plt.subplots(3, sharex=True)
axarr[0].plot(time, abs(dipole_energy))
axarr[0].plot(time, s.dipole_energy, '-r')
axarr[0].legend(('Our Python code', 'IDL code (Mark)'))
axarr[0].set_ylabel('Energy in Dipole (erg)')
axarr[1].plot(time, abs(quadrapolar_energy))
axarr[1].plot(time, s.quadrapole_energy, '-r')
axarr[1].legend(('Our Python code', 'IDL code (Mark)'))
axarr[1].set_ylabel('Energy in Quadrapole (erg)')
axarr[2].semilogy(time, abs(quad_dip_ratio))
axarr[2].semilogy(time, s.quadrapole_energy/s.dipole_energy, '-r')
axarr[2].legend(('Our Python code', 'IDL code (Mark)'))
axarr[2].set_ylabel('Quadrapolar vs Dipolar ratio')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[2].set_xticks(t2s)
axarr[2].set_xticklabels(l2s)
axarr[2].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw6.png',dpi=500)
plt.show()

#Figure 7
f, axarr = plt.subplots(3, sharex=True)
axarr[0].plot(time, coeff20)
axarr[0].plot(time, coeff20)
axarr[0].set_ylabel('L=2 m=|0| magnitude (G)')
axarr[1].plot(time, coeff21)
axarr[1].set_ylabel('L=2, m=|1| magnitude (G)')
axarr[2].plot(time, coeff22)
axarr[2].set_ylabel('L=2, m=|2| magnitude (G)')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[2].set_xticks(t2s)
axarr[2].set_xticklabels(l2s)
axarr[2].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw7.png',dpi=500)
plt.show()


#Figure 8
f, axarr = plt.subplots(4, sharex=True)
axarr[0].plot(time, coeff30)
axarr[0].set_ylabel('L=3 m=|0| magnitude (G)')
axarr[1].plot(time, coeff31)
axarr[1].set_ylabel('L=3, m=|1| magnitude (G)')
axarr[2].plot(time, coeff32)
axarr[2].set_ylabel('L=3, m=|2| magnitude (G)')
axarr[3].plot(time, coeff33)
axarr[3].set_ylabel('L=3, m=|3| magnitude (G)')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
axarr[3].set_xticks(t2s)
axarr[3].set_xticklabels(l2s)
axarr[3].set_xlabel('Time (Yrs)')
plt.savefig(path+'figw8.png',dpi=500)
plt.show()
