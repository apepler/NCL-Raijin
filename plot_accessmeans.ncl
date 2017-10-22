load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/wrf/WRFUserARW.ncl"

;  dir="/g/data2/ub7/access-s1/hc/raw_model/atmos/zg/daily/emn/"
  dir="/g/data2/ub7/access-s1/hc/raw_model/atmos/slp/daily/emn/"
  outdir="/short/eg3/asp561/"
  year1=1990
  year2=2012
  start="0501"
  numdays=125
;  lev1=500.
; First, constants

; Then, lat/lon details

  f1 = addfile(dir+"emn/da_slp_19900501_emn.nc","r")

  lat=f1->lat
  dlat=lat(2)-lat(1)
  ilat=ind(lat.ge.-50.and.lat.le.-10)
  lat2=lat(ilat)
  ff=2.*omg*sin(lat2*dtr)

  lon=f1->lon
  dlon=lon(2)-lon(1)
  ilon=ind(lon.ge.110.and.lon.le.160.)
  lon2=lon(ilon)

  stats = new((/dimsizes(lat2),dimsizes(lon2),year2-year1+1,2/),"double")
  names=(/"Mean: May (Days 1-31)","Mean: July (Days 62-92)"/) 

;  lev=f1->z1_p_level
;  ilev=get1Dindex(lev,lev1)

; Next, do the analysis

  delete(f1)

year=year1
n=0
do while(year.le.year2)
  f1 = addfile(dir+members(em)+"/da_zg_"+year+start+"_emn.nc","r")
  stats(:,:,n,1)=dim_avg_n_Wrap(f1->slp(0:30,ilon,ilat),0)
  stats(:,:,n,2)=dim_avg_n_Wrap(f1->slp(61:91,ilon,ilat),0)
  n=n+1
  year=year+1
end do



