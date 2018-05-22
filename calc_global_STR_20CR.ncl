load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

rdir="/g/data/eg3/asp561/20CRv2c/EnsMean/"
ylist=ispan(1851,2014,1)

fout=rdir+"20CR.str.mon.mean.global.nc"

a=addfile(rdir+"prmsl.mon.mean.nc","r")
lon=a->lon
lat=a->lat
time=a->time
slp=dim_avg_n_Wrap((a->prmsl)/100,2)

; Get STRI

SH=ind(lat.ge.-60.and.lat.le.-10)
NH=ind(lat.ge.10.and.lat.le.60)

STRI_SH=dim_max_n_Wrap(slp(:,SH),1)
STRI_NH=dim_max_n_Wrap(slp(:,NH),1)

STRP_SH=new(dimsizes(STRI_SH),"float")
STRP_NH=STRP_SH

; Get STRP

do n=0,dimsizes(time)-1
  indmax=maxind(slp(n,SH))
  if(.not.ismissing(indmax) ) then 
    STRP_SH(n) = (/lat(SH(indmax))/)  
  end if

  indmax=maxind(slp(n,NH))
  if(.not.ismissing(indmax) ) then
    STRP_NH(n) = (/lat(NH(indmax))/)
  end if
end do  

; Format and save

STRI_SH!0="time"
STRI_SH&time=time
STRI_SH@name = "intensity"
STRI_SH@standard_name = "STRI"
STRI_SH@long_name = "subtropical ridge intensity"
STRI_SH@units = "hPa"
copy_VarMeta(STRI_SH,STRI_NH)

STRP_SH!0="time"
STRP_SH&time=time
STRP_SH@name = "position"
STRP_SH@standard_name = "STRP"
STRP_SH@long_name = "subtropical ridge position"
STRP_SH@units = lat@units
copy_VarMeta(STRP_SH,STRP_NH)

system("rm -rf " + fout )
out=addfile(fout,"c")
out->STRI_SH=STRI_SH
out->STRP_SH=STRP_SH
out->STRI_NH=STRI_NH
out->STRP_NH=STRP_NH


