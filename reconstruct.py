#!/usr/bin/env python
import numpy as np
import os,sys
import matplotlib
matplotlib.use('Agg')
#import geopandas
from pylab import *
import pyPLUTO as pp
from subprocess import call
from AStools import *
from scipy.special import sph_harm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def read_Bfield(myfile,dir="/home/soumitra/solar_wind/potential_field_extrapolation/"):
    filename = dir+myfile
    f = open(filename,'r')
    tmp = f.readline()
    params = f.readline().split()
    nharms = int(params[0]) ; ncomps = params[1]
    nl = int((-3+np.sqrt(9+8*nharms))/2.) #5 #nharms    
    alpha = np.zeros(nharms,dtype=complex)
    ii = 0
    for n in r_[1:nl+1]:
        for m in range(n+1):
            vals = f.readline().split()
            alpha[ii] = complex(float(vals[2]),float(vals[3]))
            ii = ii + 1
    tmp=f.readline()
    beta = np.zeros(nharms,dtype=complex) 
    ii = 0
    for n in r_[1:nl+1]:
        for m in range(n+1):
            vals = f.readline().split()
            beta[ii] = complex(float(vals[2]),float(vals[3]))
            ii = ii + 1
    tmp=f.readline()
    gamma = np.zeros(nharms,dtype=complex) 
    ii = 0
    for n in r_[1:nl+1]:
        for m in range(n+1):
            vals = f.readline().split()
            gamma[ii] = complex(float(vals[2]),float(vals[3]))
            ii = ii + 1
    f.close()
    return alpha,beta,gamma

def plot_fields(br,bt,bp,theta,phi,psfile='testa.png',myfile='File not specified',cclim=40.):
    fig=figure()
    titles = ['North Pole','South Pole']
    for ii,ilat in enumerate([90,-90]):
        subplot(3,2,1+ii)
        m = Basemap(resolution='c',projection='ortho',lon_0=0.,lat_0=ilat)
        im1 = m.contourf(-phi*180./pi,-theta*180/pi+90.,br,list(np.linspace(-cclim,cclim,10)),\
                              latlon=True,cmap='RdBu_r',extend='both')
        m.drawmapboundary()
        parallels = np.arange(-80.,80,20.)
        m.drawparallels(parallels)
        title(titles[ii])
        annotate(r'$B_r$',xycoords='axes fraction',xy=(0.01,0.01))
        if (ii == 1):
            ax=gca() 
            divider=make_axes_locatable(ax)
            cax = divider.append_axes("right",size="5%",pad=0.05)
            colorbar(cax=cax)
        
    for ii,ilat in enumerate([90,-90]):
        subplot(3,2,3+ii)
        m = Basemap(resolution='c',projection='ortho',lon_0=0.,lat_0=ilat)
        im1 = m.contourf(-phi*180./pi,-theta*180/pi+90.,bp,list(np.linspace(-cclim,cclim,10)),\
                              latlon=True,cmap='RdBu_r',extend='both')
        m.drawmapboundary()
        parallels = np.arange(-80.,80,20.)
        m.drawparallels(parallels)
        annotate(r'$B_\varphi$',xycoords='axes fraction',xy=(0.01,0.01))
        if (ii == 1):
            ax=gca() 
            divider=make_axes_locatable(ax)
            cax = divider.append_axes("right",size="5%",pad=0.05)
            colorbar(cax=cax)
        
    for ii,ilat in enumerate([90,-90]):
        subplot(3,2,5+ii)
        m = Basemap(resolution='c',projection='ortho',lon_0=0.,lat_0=ilat)
        im1 = m.contourf(-phi*180./pi,-theta*180/pi+90.,bt,list(np.linspace(-cclim,cclim,10)),\
                              latlon=True,cmap='RdBu_r',extend='both')
        m.drawmapboundary()
        parallels = np.arange(-80.,80,20.)
        m.drawparallels(parallels)
        annotate(r'$B_\theta$',xycoords='axes fraction',xy=(0.01,0.01))
        if (ii == 1):
            ax=gca() 
            divider=make_axes_locatable(ax)
            cax = divider.append_axes("right",size="5%",pad=0.05)
            colorbar(cax=cax)
        
    suptitle(myfile)
    fig.savefig(psfile)

    fig=figure()
    subplot(1,3,1)
    pcolormesh(phi,theta,br) ; clim([-cclim,cclim])
    xlabel('phi') ; ylabel('theta') ; axis('tight') ; title(r'$B_r$')
    subplot(1,3,2)
    pcolormesh(phi,theta,bp) ; clim([-cclim,cclim])
    xlabel('phi') ; axis('tight') ; title(r'$B_\varphi$')
    subplot(1,3,3)
    pcolormesh(phi,theta,bt) ; clim([-cclim,cclim])
    xlabel('phi') ; axis('tight') ; title(r'$B_\theta$')
    colorbar()
    fig.savefig(psfile.split('.')[0]+'_cart.'+psfile.split('.')[1])

def extrapol_B(myfile,theta,phi,dir="/home/soumitra/solar_wind/potential_field_extrapolation/",psfile='None',R=1.5,Rb=0.7,\
                   cclim=40.):

    ###############################
    ##### Potential extrapolation #
    ###############################
    pi = np.pi ; cos = np.cos ; sin = np.sin

    Rss = 50.;
    nl = 5 #nharms    
    alpha,beta,gamma = read_Bfield(myfile,dir=dir)

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
    ii =0 
    for n in r_[1:nl+1]:
        for m in range(n+1):
            if ((n == 1) or (n == m)):
                Blmm[ii] = 0.  ; Almm[ii] = 0.
            else:
                Blmm[ii] = Blm[ii-n] ; Almm[ii] = Alm[ii-n]
            if (n == 5):
                Blmp[ii] =0. ; Almp[ii] = 0.
            else:
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
        plot_fields(br,bt,bp,theta,phi,psfile=psfile,myfile=myfile)

    return br,bt,bp

def reconstruct_B(myfile,theta,phi,dir="/home/soumitra/solar_wind/potential_field_extrapolation/",psfile='None',reg=True, cclim=5):
    ###########################################################################################
    ##### Reconstruct the components of the magnetic field on the given grid of point theta,phi 
    ###########################################################################################
    pi = np.pi ; cos = np.cos ; sin = np.sin

    nl = 5 #nharms    
    alpha,beta,gamma = read_Bfield(myfile,dir=dir)
    
    br = np.zeros(np.shape(theta))
    bt = np.zeros(np.shape(theta))
    bp = np.zeros(np.shape(theta))
    if (reg):
        ii=0
        for n in r_[1:nl+1]:
            for m in range(n+1):
                xx = xlm(m,n,phi,theta)
                yy = ylm(m,n,phi,theta)
                zz = zlm(m,n,phi,theta)
                tmp = alpha[ii]*yy
                br = br + tmp.real
                tmp = beta[ii]*zz - gamma[ii]*xx
                bt = bt + tmp.real
                tmp = beta[ii]*xx + gamma[ii]*zz
                bp = bp + tmp.real
                ii = ii + 1 
    else:    
        ii=0 ; Nt,Np = np.shape(theta)
        for n in r_[1:nl+1]:
            for m in range(n+1):
                for i in range(Nt):
                    for j in range(Np):
                        xx = xlm(m,n,[phi[i,j]],[theta[i,j]])
                        yy = ylm(m,n,[phi[i,j]],[theta[i,j]])
                        zz = zlm(m,n,[phi[i,j]],[theta[i,j]])
                        tmp = alpha[ii]*yy
                        br[i,j] = br[i,j] + tmp.real
                        tmp = beta[ii]*zz + gamma[ii]*xx
                        bt[i,j] = bt[i,j] + tmp.real
                        tmp = beta[ii]*xx - gamma[ii]*zz
                        bp[i,j] = bp[i,j] + tmp.real
                ii = ii + 1 
                        
    if (psfile != 'None'):
        plot_fields(br,bt,bp,theta,phi,psfile=psfile,myfile=myfile,cclim=cclim)

    return br,bt,bp

def compare_fields(theta,phi,BB,BB_e,cclim=40.,psfile='test_compared.png',myfile='File not specified'):

    load_j('None')
    fig=figure(figsize=(12,8))

    subplot(3,3,1)
    pcolormesh(phi,theta,BB[0],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; ylabel(r'$\theta$') ; title(r'$B_r$ (obs)')
    subplot(3,3,2)
    pcolormesh(phi,theta,BB_e[0],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; title(r'$B_r$ (pot)')
    subplot(3,3,3)
    pcolormesh(phi,theta,BB[0]-BB_e[0],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; title(r'$B_r$ (non pot)')
    ax=gca() ; divider=make_axes_locatable(ax) ;cax = divider.append_axes("right",size="5%",pad=0.05)
    colorbar(cax=cax)

    subplot(3,3,4)
    pcolormesh(phi,theta,BB[1],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; ylabel(r'$\theta$') ; title(r'$B_\theta$ (obs)')
    subplot(3,3,5)
    pcolormesh(phi,theta,BB_e[1],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; title(r'$B_\theta$ (pot)')
    subplot(3,3,6)
    pcolormesh(phi,theta,BB[1]-BB_e[1],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; title(r'$B_\theta$ (non pot)')
    ax=gca() ; divider=make_axes_locatable(ax) ;cax = divider.append_axes("right",size="5%",pad=0.05)
    colorbar(cax=cax)

    subplot(3,3,7)
    pcolormesh(phi,theta,BB[2],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; ylabel(r'$\theta$') ; xlabel(r'$\varphi$') ; title(r'$B_\varphi$ (obs)')
    subplot(3,3,8)
    pcolormesh(phi,theta,BB_e[2],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; xlabel(r'$\varphi$') ; title(r'$B_\varphi$ (pot)')
    subplot(3,3,9)
    pcolormesh(phi,theta,BB[2]-BB_e[2],cmap='RdBu_r') ; clim([-cclim,cclim])
    axis('tight') ; xlabel(r'$\varphi$') ; title(r'$B_\varphi$ (non pot)')
    ax=gca() ; divider=make_axes_locatable(ax) ;cax = divider.append_axes("right",size="5%",pad=0.05)
    colorbar(cax=cax)

    suptitle(myfile)
    fig.savefig(psfile,dpi=300)

#myfile = "hd189733V_mean_jul08_gt.c1"
myfile='WSO_1980A'
myfile_tex = myfile.replace('_','\_')
theta, phi = np.mgrid[1.e-3:pi-1.e-3:101j, 0:2*pi:201j]
#theta = np.linspace(0.,np.pi, 512)
#phi = np.linspace(0.,2.0*np.pi, 5122
br,bt,bp = reconstruct_B(myfile,theta,phi,psfile='test_reconstuct3.png')

#br_e,bt_e,bp_e = extrapol_B(myfile,theta,phi,R=1.0,Rb=1.0)

#BB_e = np.sqrt(br_e**2+bt_e**2+bp_e**2)
#BB   = np.sqrt(br**2+bt**2+bp**2)

#print 'Bmax',np.amax(BB)
#print 'Bmax (pot)',np.amax(BB_e)

#compare_fields(theta,phi,[br,bt,bp],[br_e,bt_e,bp_e],cclim=40.,myfile=myfile_tex,psfile='test_compare.png')

#[Br,Bt,Bp]=reconstruct_B(myfile,theta,phi,psfile='test_reg.png')

#tt = np.sort(np.random.uniform(0.,np.pi,(100)))
#pp = np.sort(np.random.uniform(0.,2.*np.pi,(100)))
#theta,phi = np.meshgrid(tt,pp)
#[Br,Bt,Bp]=reconstruct_B(myfile,theta,phi,psfile='test_rand.png')

#tt = np.random.uniform(0.,np.pi,(100,100))
#pp = np.random.uniform(0.,2.*np.pi,(100,100))
#for i in range(100):
#    tt[:,i] = np.sort(tt[:,i])
#    pp[i,:] = np.sort(pp[i,:])
#theta,phi = tt,pp
#[Br,Bt,Bp]=reconstruct_B(myfile,theta,phi,psfile='test_rand2.png')

#theta = np.random.uniform(0.,np.pi,(100,100))
#phi   = np.random.uniform(0.,2.*np.pi,(100,100))
#[Br,Bt,Bp]=reconstruct_B(myfile,theta,phi,psfile='test_rand.png')
#[Br,Bt,Bp]=reconstruct_B(myfile,theta,phi,psfile='test_rand2.png',reg=False)



