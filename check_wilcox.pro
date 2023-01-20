
wilcoxpath='/home/soumitra/Final_Study/harmonic-analysis/Wil/'
wilcoxfinalpath='/home/soumitra/Final_Study/harmonic-analysis/Wilcox_Magnetogram/'
file=file_search(wilcoxpath+'WSO.????.F.txt')
ncr= n_elements(file)
nlat=60
nlon=nlat*2
gaussquad_legendre,nlat,cth,weights
theta=acos(cth)  ;  radians
lat=90-theta*180/!dpi
lon=(linrange(nlon+1,0,360))(0:nlon-1)
phi=lon*!dtor

for j=0l,ncr-1 do begin
;  parse file
    openr,lun,file(j),/g
    fs=fstat(lun)
    tablebyte=bytarr(fs.size)
    readu,lun,tablebyte
    table=string(tablebyte)
    ix=strpos(table,'CT')
    posix=strpos(table,'CT',ix+1)
    repeat begin
      ix=[ix,posix]
 ;     print, ix, posix
      posix=strpos(table,'CT',posix+1)
    endrep until posix eq -1
    longitudes=round(float(strmid(table,ix+7,3)))
 ;   print, ix
    ix=ix+18
 ;   print, ix
    point_lun,lun,0
    data=fltarr(72,30)  ;  assumes 72 lon bins, 30 slat bins
    for i=0,71 do begin
      point_lun,lun,ix(i)
      buff1=fltarr(6)
      readf,lun,buff1,f='(6f9.3)'
      buff2=fltarr(8)
      readf,lun,buff2,f='(8f9.3)'
      buff3=fltarr(8)
      readf,lun,buff3,f='(8f9.3)'
      buff4=fltarr(8)
      readf,lun,buff4,f='(8f9.3)'
      data(i,*)=[buff1,buff2,buff3,buff4]
    endfor
    free_lun,lun
    data=shift(data,-(where(longitudes eq 360))(0),0)  ;  align longitudes to 0
    data=reverse(data,1)  ;  make longitude increasing
    data=reverse(data,2)  ;  make latitude (instead of colatitude) increasing
    extract= strmid(file(j), 9,4, /reverse_offset)

    ;  convert from line-of-sight to radial field
    dlonix=linrange(72,0,355)+2.5  ;  72 longitude bins
    dslatix=linrange(30,14.5,-14.5)/15  ;  30 slat bins
    dlatix=asin(dslatix)*180/!dpi
    slatgrid=replicate(1,72)#dslatix
    data=data/sqrt(1-slatgrid*slatgrid)

    ;  remap onto our grid
    dlatinterp=get_interpolation_index(reverse(dlatix),lat)
    dloninterp=get_interpolation_index(dlonix,lon)
    magout=interpolate(data,dloninterp,dlatinterp,/grid)
    save, filename=wilcoxfinalpath+'WSO.'+extract+'.F.sav', magout
    fits_write, wilcoxfinalpath+'WSO.'+extract+'.F.fits', magout
endfor

end

