load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Edir="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_pl/v01/"
outdir="/g/data/eg3/asp561/CHaines_ERAI/"
do year=1979,2016

tfile=addfiles( systemfunc (" ls -1 " + Edir + "ta/ta_6hrs_ERAI_historical_an-pl_" + year + "*.nc"), "r")
hfile=addfiles( systemfunc (" ls -1 " + Edir + "hur/hur_6hrs_ERAI_historical_an-pl_" + year + "*.nc"), "r")

time=tfile[:]->time
lev=tfile[0]->lev
lat=tfile[0]->lat
lon=tfile[0]->lon

L=ind(lev.eq.85000..or.lev.eq.70000.) ; Because ERAI is in Pa not hPa

; Now, need our three variables

AirTemp=tfile[:]->ta(:,L,:,:)
RH=hfile[:]->hur(:,L,:,:)
DewPoint=dewtemp_trh(AirTemp,RH)

CA=0.5*(AirTemp(:,0,:,:)-AirTemp(:,1,:,:)) - 2 ; 850 hPa - 750 hPa
CB=0.3333*(AirTemp(:,0,:,:)-DewPoint(:,0,:,:)) - 1 ; 850 hPa T - DP
CB=where(CB.ge.9.,9.,CB)
CB=where(CB.ge.5,5+0.5*(CB-5.),CB)
CHaines=CB+CA

CHaines!0="time"
CHaines&time=time
CHaines!1="lat"
CHaines&lat=lat
CHaines!2="lon"
CHaines&lon=lon

outfile=outdir+"CHaines_"+year+".nc"

system("/bin/rm -f "+outfile)   ; remove any pre-existing file
ncdf = addfile(outfile ,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@title         = "C-Haines index calculated as per Eqs. 1-5 in Mills & McCaw (2010 : CAWCR Tech Report No. 20)"
fAtt@source        = "ERAI-Interim 6hrly data - ta & hur at 850 & 700 hPa, 0.75x0.75 grid"
fAtt@software      = "chaines_erai.ncl"
fAtt@author        = "Acacia Pepler <acacia.pepler@bom.gov.au>"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
filedimdef(ncdf,"time",-1,True) 

copy_VarMeta(CHaines,CA)
copy_VarMeta(CHaines,CB)

CHaines@standard_name="C_Haines"
CHaines@long_name = "Continuous Haines Index"

CA@standard_name="C_Haines_A"
CA@long_name="C-Haines stability component: 850 - 700 hPa temperature difference"

CB@standard_name="C_Haines_B"
CB@long_name="C-Haines moisture component: Dewpoint depression at 850 hPa"

ncdf->CHaines=CHaines
ncdf->CA=CA
ncdf->CB=CB

delete([/outfile,fAtt,tfile,hfile,L,lon,lat,time,lev/])
delete([/AirTemp,RH,DewPoint,CHaines,CA,CB/])

end do





