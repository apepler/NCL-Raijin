load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRFUserARW.ncl"

  dir="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_pl/v01/z/"
  outdir="/g/data/eg3/asp561/GV/"
  year1=1980
  year2=2016

; First, constants

  dtr=get_pi("double")/180.
  rot=86160.
  omg=2*get_pi("double")/rot
  sc=-1e5

; Then, lat/lon details

  FILES = systemfunc (" ls -1 " + dir + "z_6hrs_ERAI_historical_an-pl_" + year1 + "*.nc")
  f1 = addfiles(FILES,"r")

  lat=f1[0]->lat
  dlat=lat(2)-lat(1)
  ilat=ind(lat.ge.-40.and.lat.le.-20)
  lat2=lat(ilat)
  ff=2.*omg*sin(lat2*dtr)

  lon=f1[0]->lon
  dlon=lon(2)-lon(1)
  ilon=ind(lon.ge.135.and.lon.le.165.)
  lon2=lon(ilon)

  dx = new((/dimsizes(lat2),dimsizes(lon2)/),"double")
  dy=dx

  do ix=0,dimsizes(lon2)-1
   do iy=0,dimsizes(lat2)-1
    dx(iy,ix) = 111.2 * 1000. * dlon * cos(lat2(iy)*dtr) 
    dy(iy,ix) = 111.2 * 1000. * dlat
   end do
  end do

  lev=f1[0]->lev
  levlist=(/1000.,925.,850.,700.,500.,300./)*100
  ilev=get1Dindex(lev,levlist)

; Next, do the analysis

  delete([/FILES,f1/])

year=year1
do while(year.le.year2)

  FILES = systemfunc (" ls -1 " + dir + "z_6hrs_ERAI_historical_an-pl_" + year + "*.nc")
  f1 = addfiles(FILES,"r")
  time=f1[:]->time

  ntimes=dimsizes(time)-1
  lap = new((/dimsizes(time)/4,dimsizes(ilev),dimsizes(lat2),dimsizes(lon2)/),"double")

  do it = 0,ntimes/4
   do il = 0,dimsizes(ilev)-1
    z=dim_avg_n_Wrap(f1[:]->z((it*4):(it*4+3),ilev(il),ilat,ilon),0)
;    z2=smth9(smth9( smth9(z,0.5,0.25,False),0.5,0.25,False),0.5,0.25,False) 
;    z2=smth9(smth9( smth9(z,0.5,0.25,False),0.5,0.25,False),0.5,0.25,False) 
    do ix=1,dimsizes(lon2)-2
     do iy=1,dimsizes(lat2)-2
      lapx=(z(iy+1,ix)+z(iy-1,ix)-2*z(iy,ix))/(1*dx(iy,ix)^2.)
      lapy=(z(iy,ix+1)+z(iy,ix-1)-2*z(iy,ix))/(1*dy(iy,ix)^2.)
      if(ff(iy).ne.0) then
        lap(it,il,iy,ix)=(lapx+lapy)*sc/ff(iy)
      end if
     end do
    end do
    delete([/z,lapx,lapy/])
   end do
  end do

; Now, save as a netcdf file for late

lap!0="time"
lap&time=yyyymmdd_time(year,year,"integer")
lap!1="lev"
lap&lev=lev(ilev)
lap!2="lat"
lap&lat=lat2
lap!3="lon"
lap&lon=lon2
 
  outfile=outdir+"ERAI_dailyGV_"+year+".nc"
  fout=addfile(outfile,"c")
fatt = True
fatt@title = "Daily geostrophic vorticity for ERA Interim"
fileattdef(fout,fatt)
fout->GV=lap

  year=year+1
 
  delete([/FILES,f1,fout,fatt,outfile,lap,time,ntimes/])
end do

