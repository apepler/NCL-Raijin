load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRFUserARW.ncl"

  dir="/g/data2/ub7/access-s1/hc/raw_model/atmos/zg/daily/"
  outdir="/g/data/eg3/asp561/GV/ACCESS/"
  year1=1990
  year2=1990
  members=(/"e01","e02","e03","e04","e05","e06","e07","e08","e09","e10","e11"/)
  start="0501"
  numdays=125

; First, constants

  dtr=get_pi("double")/180.
  rot=86160.
  omg=2*get_pi("double")/rot
  sc=-1e5
  G=9.8065

; Then, lat/lon details

  f1 = addfile(dir+"e01/da_zg_19970801_e01.nc","r")

  lat=f1->lat
  dlat=lat(2)-lat(1)
  ilat=ind(lat.ge.-40.and.lat.le.-20)
  lat2=lat(ilat)
  ff=2.*omg*sin(lat2*dtr)

  lon=f1->lon
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

  lev=f1->z1_p_level
  levlist=(/1000.,925.,850.,700.,500.,300./)
  ilev=get1Dindex(lev,levlist)

; Next, do the analysis

  delete(f1)

year=year1
do while(year.le.year2)
;do em=0,dimsizes(members)-1
em=0
  f1 = addfile(dir+members(em)+"/da_zg_"+year+start+"_"+members(em)+".nc","r")

  lap = new((/numdays,dimsizes(ilev),dimsizes(lat2),dimsizes(lon2)/),"double")

  do it = 0,numdays-1
   do il = 0,dimsizes(ilev)-1
    z=(f1->zg(it,ilev(il),ilat,ilon))*G
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

lap!0="leadtime"
lap&leadtime=ispan(0,numdays-1,1)
lap!1="lev"
lap&lev=lev(ilev)
lap!2="lat"
lap&lat=lat2
lap!3="lon"
lap&lon=lon2
 
  outfile=outdir+"ACCESS_"+members(em)+"_dailyGV_"+year+start+".nc"
  fout=addfile(outfile,"c")
fatt = True
fatt@title = "Daily geostrophic vorticity for ACCESS-S1"
fileattdef(fout,fatt)
fout->GV=lap

  delete([/f1,fout,fatt,outfile,lap,time,ntimes/])
;end do
  year=year+1
end do
