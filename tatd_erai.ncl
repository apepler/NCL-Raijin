load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Edir="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_pl/v01/"
outdir="/g/data/eg3/asp561/CHaines_ERAI/"
do year=1979,2016

print(year)
tfile=addfiles( systemfunc (" ls -1 " + Edir + "ta/ta_6hrs_ERAI_historical_an-pl_" + year + "*.nc"), "r")
hfile=addfiles( systemfunc (" ls -1 " + Edir + "hur/hur_6hrs_ERAI_historical_an-pl_" + year + "*.nc"), "r")

time=tfile[:]->time
lev=tfile[0]->lev
lat=tfile[0]->lat
lon=tfile[0]->lon

I=ind(lon.ge.110..and.lon.le.155.)
J=ind(lat.ge.-45..and.lat.le.-10.)
L=ind(lev.eq.85000..or.lev.eq.70000.) ; Because ERAI is in Pa not hPa

; Now, need our three variables

AirTemp=tfile[:]->ta(:,L,J,I)
RH=hfile[:]->hur(:,L,J,I)
DewPoint=dewtemp_trh(AirTemp,RH)
copy_VarMeta(AirTemp,DewPoint)
DewPoint@long_name="Dewpoint Temperature"
DewPoint@standard_name="dewpoint_temperature"

outfile=outdir+"TaTd_"+year+".nc"

system("/bin/rm -f "+outfile)   ; remove any pre-existing file
ncdf = addfile(outfile ,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@source        = "ERAI-Interim 6hrly data - ta & hur at 850 & 700 hPa, 0.75x0.75 grid"
fAtt@software      = "chaines_erai.ncl"
fAtt@author        = "Acacia Pepler <acacia.pepler@bom.gov.au>"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
filedimdef(ncdf,"time",-1,True) 

ncdf->ta=AirTemp
ncdf->hur=RH
ncdf->td=DewPoint

delete([/outfile,fAtt,tfile,hfile,L,lon,lat,time,lev/])
delete([/AirTemp,RH,DewPoint/])

end do





