load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRFUserARW.ncl"

  dir="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_pl/v01/z/"
  year1=1980
  year2=1980

; First, constants

  dtr=get_pi("double")/180.
  rot=86160.
  omg=2*get_pi("double")/rot
  sc=-1e5

; Then, lat/lon details

  FILES = systemfunc (" ls -1 " + dir + "z_6hrs_ERAI_historical_an-pl_" + year1 + "*.nc")
  f1 = addfiles(FILES,"r")

  lat=f1[0]->lat
  ff=2.*omg*sin(lat*dtr)
  dlat=lat(2)-lat(1)
  lon=f1[0]->lon
  dlon=lon(2)-lon(1)

  dx = new((/dimsizes(lat),dimsizes(lon)/),"double")
  dy=dx

  do ix=0,dimsizes(lon)-1
   do iy=0,dimsizes(lat)-1
    dx(iy,ix) = 111.2 * 1000. * dlon * cos(lat(iy)*dtr) 
    dy(iy,ix) = 111.2 * 1000. * dlat
   end do
  end do

  lev=f1[0]->lev
  ilev=ind(lev.eq.50000.)

; Next, do the analysis

  delete([/FILES,f1/])

year=year1
do while(year.le.year2)

  FILES = systemfunc (" ls -1 " + dir + "z_6hrs_ERAI_historical_an-pl_" + year + "*.nc")
  f1 = addfiles(FILES,"r")
  time=f1[:]->time

  ntimes=dimsizes(time)-1
  lap = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),"double")

  do it = 0,ntimes
    z=f1[:]->z(it,ilev,:,:)
    z2=smth9(smth9( smth9(z,0.5,0.25,False),0.5,0.25,False),0.5,0.25,False) 
    z2=smth9(smth9( smth9(z,0.5,0.25,False),0.5,0.25,False),0.5,0.25,False) 
    do ix=2,dimsizes(lon)-3
     do iy=2,dimsizes(lat)-3
      lapx=(z2(iy+1,ix)+z2(iy-1,ix)-2*z2(iy,ix))/(1*dx(iy,ix)^2.)
      lapy=(z2(iy,ix+1)+z2(iy,ix-1)-2*z2(iy,ix))/(1*dy(iy,ix)^2.)
      if(ff(iy).ne.0) then
        lap(it,iy,ix)=(lapx+lapy)*sc/ff(iy)
      end if
     end do
    end do
    delete([/z,lapx,lapy/])
  end do

; Finally, calculate what we want - the maximum GV between 

  I=ind(lon.ge.141.and.lon.le.156.)
  J=ind(lat.ge.-36.and.lat.le.-25)

  maxGV=dim_max_n(lap(:,J,I),(/1,2/))
  outfile="ERAI_GV500_"+year+"_v3.txt"
  asciiwrite(outfile, maxGV)
  year=year+1
 
  delete([/FILES,f1,I,J,outfile,maxGV,lap,time,ntimes/])
end do

