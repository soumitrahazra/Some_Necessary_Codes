from __future__ import print_function
from platform import python_version
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from pylab import *
import pyPLUTO as pp
import glob
from mpl_toolkits.axes_grid.inset_locator import BboxPatch, inset_axes
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.image import NonUniformImage
import copy as cp
from matplotlib.collections import Collection 
from matplotlib.artist import allow_rasterization 
from subprocess import call
from scipy.special import sph_harm
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.optimize import curve_fit
from matplotlib.image import NonUniformImage
import copy as cp
from matplotlib.collections import Collection 
from matplotlib.artist import allow_rasterization 
from subprocess import call
import scipy.interpolate
from scipy.interpolate import RectBivariateSpline

def read_Bfield(myfile,dir="/home/soumitra/solar_wind/potential_field_extrapolation/",version=1):

    filename = dir+myfile
    f = open(filename,'r')
    tmp = f.readline()
    params = f.readline().split()
    nharms = int(params[0]) ; ncomps = params[1]
    nl = int((-3+np.sqrt(9+8*nharms))/2.) # 5 #nharms    
    alpha = np.zeros(nharms,dtype=complex)
    ii = 0
    for n in r_[1:nl+1]:
        for m in range(n+1):
            vals = f.readline().split()
            alpha[ii] = complex(float(vals[2]),-float(vals[3]))
            ii = ii + 1
    tmp=f.readline()
    beta = np.zeros(nharms,dtype=complex) 
    ii = 0
    for n in r_[1:nl+1]:
        for m in range(n+1):
            vals = f.readline().split()
            beta[ii] = -complex(float(vals[2]),-float(vals[3]))
            if (version == 2):
                beta[ii] -= alpha[ii]
            ii = ii + 1
    tmp=f.readline()
    gamma = np.zeros(nharms,dtype=complex) 
    ii = 0
    for n in r_[1:nl+1]:
        for m in range(n+1):
            vals = f.readline().split()
            gamma[ii] = complex(float(vals[2]),-float(vals[3]))
            ii = ii + 1
    f.close()
    return alpha,beta,gamma


def extrapol_B(myfile,theta,phi,dir="/home/soumitra/solar_wind/potential_field_extrapolation/",psfile='None',R=1.5,Rb=0.7,\
                   cclim=40.,liml=-1):

    ###############################
    ##### Potential extrapolation #
    ###############################
    pi = np.pi ; cos = np.cos ; sin = np.sin

    Rss = 50.;
    alpha,beta,gamma = read_Bfield(myfile,dir=dir)
    if (liml < 0):
        nl = int((-3+np.sqrt(9+8*len(alpha)))/2.)
    else:
        nl = liml
    print(nl) 
    Alm = np.zeros(np.shape(alpha),dtype=complex)
    Almm = np.zeros(np.shape(alpha),dtype=complex)
    Almp = np.zeros(np.shape(alpha),dtype=complex)
    Blm = np.zeros(np.shape(alpha),dtype=complex)
    Blmm = np.zeros(np.shape(alpha),dtype=complex)
    Blmp = np.zeros(np.shape(alpha),dtype=complex)
    ii =0 
    for n in r_[1:nl+1]:
        for m in range(n+1):
            Blm[ii] = alpha[ii]/( (1.+n) * (Rb**(-(n+2.))) + n * (Rss**(-(2.*n+1.))) * (Rb**(n-1.)) )
            Alm[ii] = - (Rss**(-(2.*n+1))) * Blm[ii]
            ii=ii+1
    print(ii)
    ii =0 
    for n in r_[1:nl+1]:
        for m in range(n+1):
            if ((n == 1) or (n == m)):
                Blmm[ii] = 0.  ; Almm[ii] = 0.
            else:
                Blmm[ii] = Blm[ii-n] ; Almm[ii] = Alm[ii-n]
            if (n == 5):
                Blmp[ii] =0. ; Almp[ii] = 0.
            elif ii+n+1 < np.shape(alpha)[0]-1:
                Blmp[ii] = Blm[ii+n+1] ; Almp[ii] = Alm[ii+n+1]
            ii=ii+1
            
    br = np.zeros(np.shape(theta))
    bt = np.zeros(np.shape(theta))
    bp = np.zeros(np.shape(theta))
    ii=0
    I = complex(0.,1.)
    for n in r_[1:nl+1]:
        for m in r_[0:n+1]:
            yy = ylm(m,n,phi,theta)
            if ((n == 1) and (m == 0)):
                yym = sqrt(1./(4.*pi))
            if (n==m):
                yym = 0.
            else:
                yym = ylm(m,n-1,phi,theta)
            if (m == 0):
                coeff = 1.0
            else:
                coeff = 1.
            yyp = ylm(m,n+1,phi,theta)
            l = 1.0*n ; mm = 1.0*m ; lp = 1.0*(n+1) ; lm = 1.0*(n-1)
            
            rml  = sqrt((l*l-mm*mm)/(4.*l*l-1.))
            rmlp = sqrt((lp*lp-mm*mm)/(4.*lp*lp-1.))
            
            tmp = - yy * coeff * ( Alm[ii]*n*(R**lm) - Blm[ii]*lp*(R**(-(l+2.))) )
            br = br + tmp.real
            
            tmp = - coeff * ( Alm[ii]*(R**l ) + Blm[ii]*(R**(-lp))) * ( l*rmlp*yyp - lp*rml*yym)
            bt = bt + tmp.real/(R*sin(theta))
            
            tmp = - yy * coeff * mm * I * ( Alm[ii]*(R**l) + Blm[ii]*(R**(-lp)) )
            bp = bp + tmp.real/(R*sin(theta))
            ii = ii + 1

    if (psfile != 'None'):
        plot_fields(br,bt,bp,theta,phi,psfile=psfile)

    return br,bt,bp


def plot_fields(br, bt, bp,theta,phi,psfile=None,cclim=40.,list_lats=[90,0,-90],tits=None,mycmap='RdBu_r'):
    list_bs=[br, bt, bp]
    fig=figure(figsize=(1+2*len(list_bs),2*len(list_lats)))
    subplots_adjust(hspace=0.1,wspace=0.05)
    titles = ['Lat ='+str(v) for v in list_lats]
    nf = len(list_bs)
    nl = len(list_lats)
    nlev=10
    for ii,ilat in enumerate(list_lats):
        for iff,ff in enumerate(list_bs):
            subplot(nf,nl,1+iff*nl+ii)
            m = Basemap(resolution='c',projection='ortho',lon_0=0.,lat_0=ilat)
            im1 = m.contourf(phi*180./pi,-theta*180/pi+90.,ff,list(np.linspace(-cclim,cclim,nlev)),\
                                  latlon=True,cmap=mycmap,extend='both')
            m.drawmapboundary()
            parallels = np.arange(-80.,80,20.)
            m.drawparallels(parallels)
            if (iff == 0):
                title(titles[ii])
            if (tits != None) and (ii==0):
                annotate(tits[iff],xycoords='axes fraction',xy=(0.01,0.01))
            if (ii == nl-1):
                pos1=gca().get_position()
                xx0 = pos1.x0+pos1.width
                cax=gcf().add_axes([xx0,pos1.y0,0.02,pos1.height])
                norm = mpl.colors.Normalize(vmin=-cclim,vmax=cclim)
                cb1=mpl.colorbar.ColorbarBase(cax,cmap=mycmap,norm=norm,\
                                                  orientation='vertical',ticks=[-cclim,0,cclim],\
                                                  extend='both')
                cb1.set_label('[G]')
    
    if (psfile != None):
        fig.savefig(psfile,bbox_inches='tight')


def ylm(m,n,phi,theta):
    ll = double(n) ; mm = double(m)
    return sph_harm(m,n,phi,theta)
def xlm(m,n,phi,theta):
    ll = double(n) ; mm = double(m)
    return sph_harm(m,n,phi,theta)*1j*mm/(sin(theta)*(ll+1))
def zlm(m,n,phi,theta):
    ll = double(n) ; mm = double(m)
    cc = ((ll+1.-mm) / (ll+1.) )/sqrt( ((2.*(ll+1.)+1.0)*(ll+1-mm)) / ((2.*ll+1.0)*(ll+1.+mm)) )
    return (-sph_harm(m,n,phi,theta)*cos(theta) + cc*sph_harm(m,n+1,phi,theta))/sin(theta)


#bx, by, bz= extrapol_B('WSO_1980A', np.pi, 2*np.pi,dir="/home/soumitra/solar_wind/potential_field_extrapolation/",psfile='Non',R=1.5,Rb=0.7, cclim=40.,liml=-1)
myfile='WSO_1980A'
dir="/home/soumitra/solar_wind/potential_field_extrapolation/"
pi=22/7
alpha,beta,gamma = read_Bfield(myfile,dir=dir, version=1)
theta = np.linspace(0.,np.pi, 512)
phi = np.linspace(0.,2.0*np.pi, 512)
bx, by, bz= extrapol_B('WSO_1980A', theta=theta, phi=phi,dir="/home/soumitra/solar_wind/potential_field_extrapolation/",psfile='shaz.eps',R=1.5,Rb=0.7, cclim=40.,liml=-1)














