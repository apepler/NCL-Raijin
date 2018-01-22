load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

Edir="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_pl/v01/"
outdir="/g/data/eg3/asp561/CHaines_ERAI/"

; Regrid setup

afile=addfile("/g/data/eg3/ajd548/fire/NCCfire_grids/data_from_nina/FFDI/2000/FFDI20001002.nc","r")
lat2=afile->lat
lon2=afile->lon

Opt                = True
Opt@SrcRegional    = True             ; These are necessary
Opt@DstRegional    = True
Opt@DstGridLat = lat2
Opt@DstGridLon = lon2
Opt@ForceOverwrite = True

do year=1980,1989

tfile=addfiles( systemfunc (" ls -1 " + Edir + "ta/ta_6hrs_ERAI_historical_an-pl_" + year + "*.nc"), "r")
hfile=addfiles( systemfunc (" ls -1 " + Edir + "hur/hur_6hrs_ERAI_historical_an-pl_" + year + "*.nc"), "r")

time=tfile[:]->time
lev=tfile[0]->lev
lat=tfile[0]->lat
lon=tfile[0]->lon

I=ind(lon.ge.100..and.lon.le.170.)
J=ind(lat.ge.-60..and.lat.le.-0.)
L=ind(lev.eq.85000..or.lev.eq.70000.) ; Because ERAI is in Pa not hPa

; Now, need our three variables

ta=tfile[:]->ta(:,L,J,I)
RH=hfile[:]->hur(:,L,J,I)
td=dewtemp_trh(ta,RH)
copy_VarMeta(ta,td)

AirTemp = ESMF_regrid(ta,Opt)
AirTemp!0="time"
AirTemp!1="lev"
AirTemp!2="lat"
AirTemp!3="lon"
AirTemp&time=time
AirTemp&lev=ta&lev
AirTemp&lat=lat2
AirTemp&lon=lon2
AirTemp@units="K"
AirTemp@standard_name="air_temperature"
AirTemp@long_name="Temperature"

DewPoint = ESMF_regrid(td,Opt) 
copy_VarMeta(AirTemp,DewPoint)
DewPoint@standard_name="dew_point_temperature"
DewPoint@long_name="Dewpoint Temperature, calculated in ncl using dewtemp_trh" 

; Now, need to regrid them

CA=0.5*(AirTemp(:,0,:,:)-AirTemp(:,1,:,:)) - 2 ; 850 hPa - 750 hPa
CB=0.3333*(AirTemp(:,0,:,:)-DewPoint(:,0,:,:)) - 1 ; 850 hPa T - DP
CB=where(CB.ge.9.,9.,CB)
CB=where(CB.ge.5,5+0.5*(CB-5.),CB)
CHaines=CB+CA

CHaines!0="time"
CHaines&time=time
CHaines!1="lat"
CHaines&lat=lat2
CHaines!2="lon"
CHaines&lon=lon2

outfile=outdir+"CHaines_"+year+"_awapgrid.nc"

system("/bin/rm -f "+outfile)   ; remove any pre-existing file
ncdf = addfile(outfile ,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@title         = "C-Haines index calculated as per Eqs. 1-5 in Mills & McCaw (2010 : CAWCR Tech Report No. 20)"
fAtt@source        = "ERAI-Interim data - ta & hur at 850 & 700 hPa, bilinearly regridded to the 0.05 degree AWAP grid using ESMF in ncl"
fAtt@software      = "chaines_erai_awapgrid.ncl"
fAtt@author        = "Acacia Pepler <acacia.pepler@bom.gov.au>"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
filedimdef(ncdf,"time",-1,True) 

copy_VarMeta(CHaines,CA)
copy_VarMeta(CHaines,CB)

CHaines@standard_name="C_Haines"
CHaines@long_name = "Continuous Haines index"

CA@standard_name="C_Haines_A"
CA@long_name="C-Haines stability component: 850 - 700 hPa temperature difference"

CB@standard_name="C_Haines_B"
CB@long_name="C-Haines moisture component: Dewpoint depression at 850 hPa"


ncdf->ta=AirTemp
ncdf->td=DewPoint
ncdf->CHaines=CHaines
ncdf->CA=CA
ncdf->CB=CB

delete([/outfile,fAtt,tfile,hfile,I,J,L,lon,lat,time,lev/])
delete([/AirTemp,RH,DewPoint,CHaines,CA,CB,ta,td/])

end do





