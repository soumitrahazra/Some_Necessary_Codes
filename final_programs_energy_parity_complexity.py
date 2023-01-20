import numpy as np
from astropy.io import fits
import scipy.special as scisp
from scipy import interpolate
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
        else:
            allFiles.append(fullPath)
    
    return allFiles

#dirname= '/home/soumitra/Final_Study/data'
#dirname= '/home/soumitra/Final_Study/MW_data'
#dirname= '/home/soumitra/Final_Study/MDI_Data/Final_MDI'
dirname= '/home/soumitra/Final_Study/harmonic-analysis/Wilcox_Magnetogram'
#path= '/home/soumitra/Final_Study/figure/MW_figure/'
#path= '/home/soumitra/Final_Study/figure/MDI/'
path= '/home/soumitra/Final_Study/Final_Plot_Sacha/figures/'
#listoffiles=sorted(getListOfFiles(dirname))
#listoffiles=sorted(glob.glob(dirname+'/*/*.fits'))
listoffiles=sorted(glob.glob(dirname+'/*.F.fits'))
#listoffiles=sorted(glob.glob(dirname+'/*.fits'))
#listoffiles= sorted(glob.glob("/home/soumitra/Final_Study/data/data2/*.fits"))
ns=len(listoffiles)-1
print(ns)
l_max = 40
nb_modes_tot = (l_max+1)*(l_max+2)//2 - 1 # computing total number of modes

if sys.version_info[0] == 3:
    xrange = range

# defining the necessary arrays
total_surface_energy= np.zeros(ns)
total_axial_surface_energy= np.zeros(ns)
total_nonaxial_surface_energy=np.zeros(ns)
dipolar_coeff=np.zeros(ns)
fdipole= np.zeros(ns)
dipolar_coeffss=np.zeros(ns)
fdipoles= np.zeros(ns)
faxial=np.zeros(ns)
fnonaxial=np.zeros(ns)
time = np.zeros(ns)
times = np.zeros(ns)
parity=np.zeros(ns)
symenergy=np.zeros(ns)
antisymenergy=np.zeros(ns)
complexity= np.zeros(ns)

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
    input_data=np.nan_to_num(input_data)
    wilcox2gauss=0.01
    input_data=wilcox2gauss*input_data
    # interpolation is necessary only for MDI data
    #x= np.arange(0,input_data.shape[1], 1)
    #y= np.arange(0,input_data.shape[0], 1)
    #f= interpolate.interp2d(x,y, input_data)
    #xnew= np.arange(0,input_data.shape[1]/6, 1)
    #ynew= np.arange(0,input_data.shape[0]/6, 1)
    #input_data= f(xnew, ynew)
    #input_data=np.nan_to_num(input_data)

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
    coefbrs = np.zeros(nb_modes_tot, dtype=np.complex)
    ylm = np.zeros((nb_modes_tot, nb_th, nb_phi), dtype=np.complex)
    subb= np.zeros(((l_max +1)*(l_max +1), (l_max +1)*(l_max +1)), dtype=np.complex)
    subp= np.zeros(((l_max +1)*(l_max +1), (l_max +1)*(l_max +1)), dtype=np.complex)
    mod = 0
    subbs=0
    subps=0
    subes=0
    subos=0
    subss=0
    for l in xrange(1, l_max+1):
        for m in xrange(0, l+1):
          #print('{} {} {}'.format(l, m, mod))
          ylm[mod] = scisp.sph_harm(m, l, Phi, Theta)
          ylm_c = np.conj(ylm[mod])
          integrand_a = Br*ylm_c
          integrand_a = integrand_a*np.sin(Theta)*dphi*dtheta
          coefbr[mod] = np.sum(integrand_a)
          coefbrs[mod] = np.sum(integrand_a)
          subb[l, m] = coefbr[mod]
          subp[l,m]= coefbrs[mod]
          subbs= subbs + (np.real(subb[l,m]))**2 + (np.imag(subb[l,m]))**2
          subss=subss + l*((np.real(subb[l,m]))**2 + (np.imag(subb[l,m]))**2)
          numd= l +m
          if (numd % 2) == 0:
             subes=subes + (np.real(subb[l,m]))**2 + (np.imag(subb[l,m]))**2
          else:
             subos=subos + (np.real(subb[l,m]))**2 + (np.imag(subb[l,m]))**2
          if (m!=0):
             coefbrs[mod]= 0.0 
             subp[l,m]=coefbrs[mod]
          subps=subps + (np.real(subp[l,m]))**2 + (np.imag(subp[l,m]))**2
          mod = mod+1
    total_surface_energy[i]= subbs
    complexity[i]=subss/subbs
    total_axial_surface_energy[i]= subps
    total_nonaxial_surface_energy[i]= total_surface_energy[i]-total_axial_surface_energy[i]
    dipolar_coeff[i]=(np.real(subb[1,0]))**2 + (np.imag(subb[1,0]))**2 + (np.real(subb[1,1]))**2 + (np.imag(subb[1,1]))**2
    fdipole[i]= dipolar_coeff[i]/ total_surface_energy[i]
    dipolar_coeffss[i]=(np.real(subb[1,0]))**2 + (np.imag(subb[1,0]))**2 
    fdipoles[i]= dipolar_coeffss[i]/ total_axial_surface_energy[i]
    faxial[i]= total_axial_surface_energy[i]/total_surface_energy[i]
    fnonaxial[i]=total_nonaxial_surface_energy[i]/total_surface_energy[i]
    parity[i]= (subes-subos)/ (subes +subos)
    symenergy[i]= subes
    antisymenergy[i]=subos
    
    time[i] = 1642.0 +i
    times[i]= 1976.403 + i*(27.27/365.25)
    #time[i] = 1910.0 +i
   


t2s=[1624, 1691, 1758, 1825, 1892, 1958, 2025, 2092, 2159]
l2s=[1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015]
#t2s=[1910, 1950, 2000, 2050, 2100]
#l2s=[1996, 1999, 2003, 2006, 2010]

xtimes=np.loadtxt("sunspot_data.txt")[:,0]
yss= np.loadtxt("sunspot_data.txt")[:,1]

#Figure

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(time, total_surface_energy)
ax1.set_ylabel('Total Surface energy')
ax1.set_xlabel('Time (years)')
ax2 = ax1.twinx()
ax2.plot(time, fdipole, 'r-')
ax2.set_ylabel('$f_{dip}$', color='r', fontsize=16)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_xticks(t2s)
ax2.set_xticklabels(l2s)

#plt.savefig('surf_fdipoleaasy.png', dpi=500)
plt.show()

#Figure

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(time, total_axial_surface_energy)
ax1.set_ylabel('Total Axial Surface Energy')
ax1.set_xlabel('Time (years)')
ax2 = ax1.twinx()
ax2.plot(time, fdipoles, 'r-')
ax2.set_ylabel('$f_{axial-dip}$', color='r', fontsize=16)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_xticks(t2s)
ax2.set_xticklabels(l2s)

#plt.savefig('axial_surf_fdipoleabbssy.png', dpi=500)
plt.show()

#Figure

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(time, total_surface_energy)
ax1.set_ylabel('Total Surface Energy')

ax2 = ax1.twinx()
ax2.plot(time, total_axial_surface_energy, 'r-')
ax2.set_ylabel('Total Axial Surface Energy', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_xticks(t2s)
ax2.set_xticklabels(l2s)

#plt.savefig('axial_surf_all.png', dpi=500)
plt.show()

#Figure

fig=plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(time, total_surface_energy, '-k')
ax1.set_ylabel('Total Surface energy', fontsize=14)
#ax1.set_xlabel('Time (years)', fontsize=14)
ax2 = ax1.twinx()
ax2.plot(time, fdipole, 'r-')
ax2.set_ylabel('$f_{dip}$', color='r', fontsize=16)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_xticks(t2s)
ax2.set_xticklabels(l2s)
ax3 = fig.add_subplot(212)
ax3.plot(time, total_axial_surface_energy, '-k')
ax3.set_ylabel('Total Axial Surface Energy', fontsize=14)
ax3.set_xlabel('Time (years)', fontsize=14)
ax4 = ax3.twinx()
ax4.plot(time, faxial, 'r-')
ax4.set_ylabel('$f_{axial}$', color='r', fontsize=16)
for tl in ax4.get_yticklabels():
    tl.set_color('r')
ax4.set_xticks(t2s)
ax4.set_xticklabels(l2s)

plt.show()


#Figure

fig=plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(time, total_surface_energy, '-k')
ax1.set_ylabel('Total Surface energy', fontsize=14)
#ax1.set_xlabel('Time (years)', fontsize=14)
ax2 = ax1.twinx()
ax2.plot(time, faxial, 'r-')
ax2.set_ylabel('$f_{axisym}$', color='r', fontsize=16)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_xticks(t2s)
ax2.set_xticklabels(l2s)
ax3 = fig.add_subplot(212)
ax3.plot(time, symenergy, '-k')
ax3.set_ylabel('Symmetric Energy', fontsize=14)
ax3.set_xlabel('Time (years)', fontsize=14)
ax4 = ax3.twinx()
ax4.plot(time, antisymenergy, 'r-')
ax4.set_ylabel('Antisymmetric Energy', color='r', fontsize=16)
for tl in ax4.get_yticklabels():
    tl.set_color('r')
ax4.set_xticks(t2s)
ax4.set_xticklabels(l2s)

plt.show()

#Figure

fig=plt.figure()
ax1 = fig.add_subplot(311)
ax1.plot(times, total_surface_energy, '-k')
ax1.set_ylabel('$\sum_{\ell,m} \\alpha_{\ell,m}^2$ (erg/cm)', fontsize=14)
#ax1.set_xlabel('Time (years)', fontsize=14)
ax2 = ax1.twinx()
ax2.plot(times, faxial, 'r-')
ax2.set_ylabel('$f_{axisym}$', color='r', fontsize=16)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
#ax2.set_xticks(t2s)
#ax2.set_xticklabels(l2s)
ax3 = fig.add_subplot(312)
ax3.plot(times, total_surface_energy, '-k')
ax3.set_ylabel('$\sum_{\ell,m} \\alpha_{\ell,m}^2$ (erg/cm)', fontsize=14)
ax4 = ax3.twinx()
ax4.plot(times, parity, 'r-')
ax4.set_ylabel('Parity factor (P)', color='r', fontsize=16)
for tl in ax4.get_yticklabels():
    tl.set_color('r')
#ax4.set_xticks(t2s)
#ax4.set_xticklabels(l2s)
ax5 = fig.add_subplot(313)
ax5.plot(times, total_surface_energy, '-k')
ax5.set_ylabel('$\sum_{\ell,m} \\alpha_{\ell,m}^2$ (erg/cm)', fontsize=14)
ax5.set_xlabel('Time (years)', fontsize=16)
ax6 = ax5.twinx()
ax6.plot(times, complexity, 'r-')
ax6.set_ylabel('$l_c$', color='r', fontsize=16)
for tl in ax6.get_yticklabels():
    tl.set_color('r')
#ax6.set_xticks(t2s)
#ax6.set_xticklabels(l2s)

plt.show()

#plt.plot(time, total_surface_energy, '-k')
#plt.plot(time, fdipole, '-r')
#plt.xticks(t2s)
#plt.xticklabels(l2s)


#f= open("dataset_axial_40.txt", "w")
#f.write("# time total_surface_energy  fdipole  total_axial_surface_energy  fdipole(axial) \n")
#np.savetxt(f, np.array([time, total_surface_energy, fdipole, total_axial_surface_energy, fdipoles]).T)



#Figure 1
#f, axarr = plt.subplots(3, sharex=True)
#axarr[0].plot(time, -axial_dipole_coef)
#axarr[0].axhline(y=0)
#axarr[0].set_ylabel('Axial Dipole Coefficient (G)')
#axarr[1].plot(time, equatorial_dipole_coef)
#axarr[1].set_ylabel('Equatorial Dipole Magnitude (G)')
#axarr[2].semilogy(time, abs(equatorial_axial))
#axarr[2].set_ylabel('Equatorial to Axial Ratio')
#ax1.get_shared_x_axes().join(ax1, ax2,ax3)
#axarr[2].set_xticks(t2s)
#axarr[2].set_xticklabels(l2s)
#axarr[2].set_xlabel('Time (Yrs)')
#plt.savefig(path+'fig1.png',dpi=500)
#plt.show()


