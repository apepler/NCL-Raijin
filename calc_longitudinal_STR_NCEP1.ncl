load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

rdir="/g/data/eg3/asp561/NCEP1/"
ylist=ispan(1979,2016,1)

fout=rdir+"NCEP1.str.mon.mean.longitude2.nc"

a=addfile(rdir+"slp.mon.mean.nc","r")
lon=tofloat(a->lon)
lat=tofloat(a->lat)
time=a->time
slp=a->slp

winwid=1 
if(winwid.gt.0) then
  slp=runave_n_Wrap(slp,1+2*winwid,-1,2)
end if

; Get STRI

SH=ind(lat.ge.-45.and.lat.le.-10)
NH=ind(lat.ge.10.and.lat.le.45)

STRI_SH=dim_max_n_Wrap(slp(:,SH,:),1)
STRI_NH=dim_max_n_Wrap(slp(:,NH,:),1)

STRP_SH=new(dimsizes(STRI_SH),"float")
STRP_NH=STRP_SH

; Get STRP

do n=0,dimsizes(time)-1
 do i=0,dimsizes(lon)-1
  indmax=maxind(slp(n,SH,i))
  if(.not.ismissing(indmax) ) then 
    STRP_SH(n,i) = (/lat(SH(indmax))/)  
  end if

  indmax=maxind(slp(n,NH,i))
  if(.not.ismissing(indmax) ) then
    STRP_NH(n,i) = (/lat(NH(indmax))/)
  end if
 end do
end do  

; Format and save

STRI_SH!0="time"
STRI_SH&time=time
STRI_SH!1="lon"
STRI_SH&lon=lon
STRI_SH@name = "intensity"
STRI_SH@standard_name = "STRI"
STRI_SH@long_name = "subtropical ridge intensity by longitude "
STRI_SH@units = "hPa"
copy_VarMeta(STRI_SH,STRI_NH)

STRP_SH!0="time"
STRP_SH&time=time
STRP_SH!1="lon"
STRP_SH&lon=lon
STRP_SH@name = "position"
STRP_SH@standard_name = "STRP"
STRP_SH@long_name = "subtropical ridge position by longitude "
STRP_SH@units = lat@units
copy_VarMeta(STRP_SH,STRP_NH)

system("rm -rf " + fout )
out=addfile(fout,"c")
out->STRI_SH=STRI_SH
out->STRP_SH=STRP_SH
out->STRI_NH=STRI_NH
out->STRP_NH=STRP_NH


