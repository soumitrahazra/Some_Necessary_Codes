;  This routine generates plots of spherical harmonic modes as a function of
;  time for all Wilcox data that were available as of 2/21/2014.

;  Movie of photosphere B is stored in brbot.

;  Movie of source surface B is stored in brtop.

;  needs SSW for PFSS package

;  IDL commands to quickly do downloading (change to Wilcox directory first!!!)
;  (Note last CR for harmonic coefficients is usually empty.)

;    CRstart=2150  &  CRend=2179
;    for i=CRstart,CRend do spawn,'wget http://wso.stanford.edu/synoptic/WSO.'+string(i,f="(i4)")+'.F.txt'
;    for i=CRstart,CRend do spawn,'wget --output-document=CR'+string(i,f="(i4)")+'_001.txt http://wso.stanford.edu/Harmonic.rad/CR'+string(i,f="(i4)")

;  plots only flag
plotsonly=0

if not plotsonly then begin

  ;  PFSS common block
  @pfss_data_block

  ;  preliminaries
  nlat=60
  wilcoxpath='/home/soumitra/Final_Study/harmonic-analysis/Wilcox/'
  wilcox2gauss=0.01  ;  conversion factor from Wilcox's microTesla to Gauss
  cr1=1642e
  cr2=2209
  ncr=cr2-cr1+1
  xr=[1975,2020]

  ;  dimension arrays
  lcentarr=fltarr(ncr)
  crlist=lonarr(ncr)
  datelist=lon64arr(ncr)
  brcoeff=complexarr([nlat+1,nlat+1,ncr])
  magmovie=fltarr(nlat*2,nlat,ncr)
  brbot=fltarr(nlat*2,nlat,ncr)
  brtop=fltarr(nlat*2,nlat,ncr)
  hcgarr=fltarr(10,10,ncr)
  hcharr=fltarr(10,10,ncr)

  ;  compute array of l values
  nax=[nlat+1,nlat+1]
  larr=dindgen(nax(0))#replicate(1,nax(1))
  marr=replicate(1,nax(0))#dindgen(nax(1))

  ;  get list of valid files
  flist=file_search(wilcoxpath+'WSO.????.F.txt')
  hclist=file_search(wilcoxpath+'CR????_001.txt')
  if flist(0) eq '' then stop,'  No valid files found.'
  crlist1=long(strmid(flist,strlen(wilcoxpath)+4,4))
  crlist2=long(strmid(hclist,strlen(wilcoxpath)+2,4))

  ;  loop through files
  buffer=fltarr(10,10)
  for i=0l,ncr-1 do begin

    ;  print time left
   ; print_time,'  ',i+1,ncr,tst,slen;  &  print

    ;  create magnetogram
    whfile=where(crlist1 eq (cr1+i),nwh)
    if nwh eq 0 then stop,'  Can''t find map for CR '+strcompress(cr1+i,/r)
    pfss_mag_create,magout,1,nlat,file=flist(whfile),/quiet

    ;  get harmonic coefficients
    pfss_get_potl_coeffs,magout*wilcox2gauss,rtop=2.5,/quiet

    ;  do extrapolation
    pfss_potl_field,2.5,2,/trunc,/quiet

    ;  load in harmonic coefficients
    whfile=where(crlist2 eq (cr1+i),nwh)
    if nwh eq 0 then stop,'  Can''t find coeffs for CR '+strcompress(cr1+i,/r)
    openr,lun,hclist(whfile),/get
    for j=0,9 do begin
      minibuffer=fltarr(j+2)
      readf,lun,minibuffer
      buffer(0:j,j)=minibuffer(1:j+1)
    endfor
    hcgarr(*,*,i)=buffer*wilcox2gauss
    for j=0,9 do begin
      minibuffer=fltarr(j+2)
      readf,lun,minibuffer
      buffer(0:j,j)=minibuffer(1:j+1)
    endfor
    hcharr(*,*,i)=buffer*wilcox2gauss
    close,lun
    free_lun,lun

    ;  save data  
    magmovie(*,*,i)=magout
    crlist(i)=cr1+i
    brcoeff(*,*,i)=phibt*(larr+1)-phiat*(larr)
    brbot(*,*,i)=br(*,*,0)
    brtop(*,*,i)=br(*,*,nr-1)
    datelist(i)=long64(strmid(strjoin(strsplit( $
      anytim(carr2ex(cr1+i),/ccs),'-',/extract)),0,8))

  endfor

  ;  restore same variables from forecaster
  ;restore,'harmonic_power_forecaster.sav'

endif

;  compute phimask (originating from the way IDL does the FFT: m=0
;  modes have a weight of 1; |m|>0 modes have a weight of 0.5).  See
;  harmonic_power_test.pro for more info (and some tests...).
phimask=make_array(dim=size(phiat,/dim),val=0.5)
phimask(*,0)=1.0

;  compute factor used to convert from my normalization to
;  Wilcox's (called k in the annotated version of Xudong's notes) =
;  sqrt( (1/2) * ((l+m)!/(l-m)!) ) for m>0 (otherwise factor=1)
factorlmax=5
twol=dindgen(2*(factorlmax+1))
factor=make_array(dim=[factorlmax+1,factorlmax+1],val=1d0)
for i=1,factorlmax do for j=1,i do $
  factor(i,j)=sqrt(0.5*product(twol((i-j+1)>0:i+j)))

;  compute decimal dates from datelist
yearlist=long(datelist/10000ll)
monthlist=long((datelist mod 10000ll)/100ll)
daylist=((datelist mod 10000ll)/1e2 - monthlist)/0.31
decyear=yearlist+(monthlist+daylist-1)/12

;  determine latitude and longitude of dipole - note latitude has minus sign
;  because magnetograms are ordered from south to north
br2energy=make_array(dim=size(brcoeff,/dim),/float)
for i=0l,ncr-1 do br2energy(*,*,i)=phimask*abs(brcoeff(*,*,i))^2
diplat=reform(90-(atan(sqrt(br2energy(1,1,*)),double(brcoeff(1,0,*))))*!radeg)
diplon=reform(atan(imaginary(brcoeff(1,1,*)),double(brcoeff(1,1,*)))*!radeg)
diplon=(diplon+360) mod 360

;  get a gaussian smoothing kernel
gwid=4  ;  full-width at half-max of the gaussian
thresh=1e-16  ;  threshold to avoid underflow errors
sigma=gwid/(2.*sqrt(alog(2)))  ;  sigma in gaussian formula
khwid=long(1+(sqrt(-alog(thresh)/2)*sigma))
kwid=khwid*2+1  ;  size of kernel array, guaranteed to be odd
kernel=exp(-((findgen(kwid)-khwid)/sigma)^2)
kernel=kernel/total(kernel)  ;  normalize

;  add dipole axes vectorially
dipz=sin(diplat*!dtor)
dipx=cos(diplat*!dtor)*cos(diplon*!dtor)
dipy=cos(diplat*!dtor)*sin(diplon*!dtor)
dipz_sm=convol(dipz,kernel,/edge_zero)
dipx_sm=convol(dipx,kernel,/edge_zero)
dipy_sm=convol(dipy,kernel,/edge_zero)
diplat_sm=atan(dipz_sm,sqrt(dipx_sm*dipx_sm+dipy_sm*dipy_sm))*!radeg
diplon_sm=atan(dipy_sm,dipx_sm)*!radeg
diplon_sm=(diplon_sm+360) mod 360

dipole_energy= total(br2energy(1,0:1,*),2)
quadrapole_energy= total(br2energy(2,0:2,*),2)

;  plot some stuff
;pspage,8,10,/fonts
;mydevice= !D.NAME
set_plot, 'ps'
device,file='harmonic_power_wilcox.ps',/color
loadct,12

;  plot some stuff to match Todd's plots, similar to the top two panels of
;  figure 2 of Hoeksema 2010 (IAU Symp 264 proceedings)
!p.multi=[0,1,3]
plot,decyear,-float(brcoeff(1,0,*)),tit='VARIATION OF POLAR DIPOLE COEFFICIENT',xr=xr,/xsty,yr=[-5,5],ytit='[Gauss]',/ysty
;oplot,decyear_forecaster,float(brcoeff_forecaster(1,0,*)),col=64
oplot,decyear,-float(brcoeff(1,0,*))
oplot,xr,[0,0]
;xyouts,mean(xr),-4.75,align=0.5,'BLACK = WSO synoptic data / GREEN = LMSAL forecaster'
plot,decyear,sqrt(br2energy(1,1,*)),tit='VARIATION OF EQ DIPOLE MAGNITUDE',xr=xr,/xsty,yr=[0,5],ytit='[Gauss]'
;oplot,decyear_forecaster,sqrt(br2energy_forecaster(1,1,*)),col=64
oplot,decyear,sqrt(br2energy(1,1,*))
plot,decyear,sqrt(br2energy(1,1,*)/br2energy(1,0,*)),tit='RATIO OF EQ/POLAR DIPOLE MAGNITUDES',/ylog,xr=xr,/xsty
;oplot,decyear_forecaster,sqrt(br2energy_forecaster(1,1,*)/br2energy_forecaster(1,0,*)),col=64
oplot,decyear,sqrt(br2energy(1,1,*)/br2energy(1,0,*))
!p.multi=-1
oplot,xr,[1,1]

;  plot dipole component comparison
!p.multi=[0,1,3]
plot,decyear,sqrt(br2energy(1,0,*)),tit='VARIATION OF POLAR DIPOLE MAGNITUDE',yr=[0,5],xr=xr,/xsty,ytit='[Gauss]'
;oplot,decyear_forecaster,sqrt(br2energy_forecaster(1,0,*)),col=64
oplot,decyear,sqrt(br2energy(1,0,*))
plot,decyear,sqrt(br2energy(1,1,*)),tit='VARIATION OF EQ DIPOLE MAGNITUDE',yr=[0,5],xr=xr,/xsty,ytit='[Gauss]'
;oplot,decyear_forecaster,sqrt(br2energy_forecaster(1,1,*)),col=64
oplot,decyear,sqrt(br2energy(1,1,*))
plot,decyear,sqrt(total(br2energy(1,0:1,*),2)),tit='VARIATION OF TOTAL DIPOLE MAGNITUDE',yr=[0,5],xr=xr,/xsty,ytit='[Gauss]'
;oplot,decyear_forecaster,sqrt(total(br2energy_forecaster(1,0:1,*),2)),col=64
oplot,decyear,sqrt(total(br2energy(1,0:1,*),2))
!p.multi=-1

;  energy plots
!p.multi=[0,1,4]
plot,decyear,total(br2energy(1,0:1,*),2),tit='VARIATION OF DIPOLE ENERGY',xr=xr,/xsty,yr=[0,2.5e1]
;oplot,decyear_forecaster,total(br2energy_forecaster(1,0:1,*),2),col=64
oplot,decyear,total(br2energy(1,0:1,*),2)
plot,decyear,total(br2energy(2,0:2,*),2),tit='VARIATION OF QUADRUPOLE ENERGY',xr=xr,/xsty,yr=[0,2.5e1]
;oplot,decyear_forecaster,total(br2energy_forecaster(2,0:2,*),2),col=64
oplot,decyear,total(br2energy(2,0:2,*),2)
plot,decyear,total(br2energy(2,0:2,*),2)/total(br2energy(1,0:1,*),2),tit='RATIO OF QUADRUPOLE/DIPOLE ENERGY',xr=xr,/xsty,/ylog
;oplot,decyear_forecaster,total(br2energy_forecaster(2,0:2,*),2)/total(br2energy_forecaster(1,0:1,*),2),col=64
oplot,decyear,total(br2energy(2,0:2,*),2)/total(br2energy(1,0:1,*),2)
oplot,xr,[1,1]
plot,total(total(br2energy,3)/ncr,2),tit='TIME-AVERAGED ENERGY IN EACH DEGREE L',/xsty,/ylog,psym=-4,xtit='spherical harmonic degree L',yr=[1e-2,1e3]
;oplot,total(total(br2energy_forecaster,3)/ncr,2),col=64,psym=-4
oplot,total(total(br2energy,3)/ncr,2)
!p.multi=-1

;  plot comparison for (l,m)=(1,0) mode
!p.multi=[0,1,2]
plot,decyear,-sqrt(0.75/!dpi)*float(brcoeff(1,0,*))/factor(1,0),title='g_(l=1,m=0)',yr=[-2.5,2.5],xr=xr,/xsty
oplot,decyear,hcgarr(0,1,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-2.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
plot,decyear,((hcgarr(0,1,*))/(-sqrt(0.75/!dpi)*float(brcoeff(1,0,*)))),yr=[0.6,1.4],title='RATIO OF ABOVE',xr=xr,/xsty,psym=4
oplot,xr,[1,1]
!p.multi=-1

;  plot comparison for (l,m)=(1,1) mode
!p.multi=[0,1,4]
loadct,12
plot,decyear,sqrt(0.375/!dpi)*float(brcoeff(1,1,*))/factor(1,1),title='g_(l=1,m=1)',yr=[-2,2],xr=xr,/xsty
oplot,decyear,hcgarr(1,1,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-1.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
plot,decyear,(hcgarr(1,1,*))/(sqrt(0.375/!dpi)*float(brcoeff(1,1,*))/factor(1,1)),yr=[0.6,1.4],title='RATIO OF ABOVE',xr=xr,/xsty,psym=4
oplot,xr,[1,1]
plot,decyear,-sqrt(0.375/!dpi)*imaginary(brcoeff(1,1,*))/factor(1,1),title='h_(l=1,m=1)',xr=xr,/xsty
oplot,decyear,hcharr(1,1,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-1.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
plot,decyear,(hcharr(1,1,*))/(-sqrt(0.375/!dpi)*imaginary(brcoeff(1,1,*))/factor(1,1)),yr=[0.6,1.4],title='RATIO OF ABOVE',xr=xr,/xsty,psym=4
oplot,xr,[1,1]

!p.multi=-1

;  plot comparison for (l,m)=(2,[0,1,2]) modes
!p.multi=[0,1,5]
loadct,12
plot,decyear,-sqrt(1.25/!dpi)*float(brcoeff(2,0,*))*factor(2,0),title='g_(l=2,m=0)',yr=[-2,2],xr=xr,/xsty
oplot,decyear,hcgarr(0,2,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-1.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
plot,decyear,sqrt(5./(24.*!dpi))*float(brcoeff(2,1,*))*factor(2,1),title='g_(l=2,m=1)',yr=[-2,2],xr=xr,/xsty
oplot,decyear,hcgarr(1,2,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-1.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
plot,decyear,-sqrt(5./(24.*!dpi))*imaginary(brcoeff(2,1,*))*factor(2,1),title='h_(l=2,m=1)',yr=[-2,2],xr=xr,/xsty
oplot,decyear,hcharr(1,2,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-1.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
plot,decyear,-sqrt(5./(96.*!dpi))*float(brcoeff(2,2,*))*factor(2,2),title='g_(l=2,m=2)',yr=[-2,2],xr=xr,/xsty
oplot,decyear,hcgarr(2,2,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-1.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
plot,decyear,sqrt(5./(96.*!dpi))*imaginary(brcoeff(2,2,*))*factor(2,2),title='h_(l=2,m=2)',yr=[-2,2],xr=xr,/xsty
oplot,decyear,hcharr(2,2,*),col=128
oplot,xr,[0,0]
xyouts,mean(xr),-1.9,align=0.5,'BLACK = my reconstruction from WSO synoptic data, RED = from WSO website'
!p.multi=-1

device,/cl
set_plot,'x'

spawn,'ps2pdf harmonic_power_wilcox.ps'

;  save some stuff
save,file='harmonic_power_wilcox.sav',decyear,brcoeff,br2energy,phiat,phibt, $
  hcgarr,hcharr,diplat,diplon,magmovie,crlist,diplat_sm,diplon_sm,phi,theta, phimask, dipole_energy, quadrapole_energy

end
