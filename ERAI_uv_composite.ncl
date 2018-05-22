load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

DATADIR="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_sfc/v01/"
mdir="/g/data/eg3/asp561/ERAI/"

dir="/short/eg3/asp561/cts.dir/gcyc_out/"
reanal="ERAI"
proj="proj100"
thresh="rad10cv0.075"
type="highs"
edir=dir+reanal+"/"+proj+"_"+type+"_"+thresh+"/"
filename=edir+"UM_"+type+"_"+reanal+"_"+proj+"_"+thresh+"_bigaust_fixes.csv"
fout=edir+"UM_"+type+"_"+reanal+"_"+proj+"_"+thresh+"_bigaustwind.nc"

latlim=(/-50,-10/)
lonlim=(/90,180/)

; This is the dimensions of the output grid

  olon=fspan(-21,21,57)
  olon!0="lon"
  olon@units="degrees_E"
  olat=fspan(-21,21,57)
  olat!0="lat"
  olat@units="degrees_N"

; Set up ERAI

fin=addfile(DATADIR+"tas/tas_6hrs_ERAI_historical_an-sfc_20171201_20171231.nc","r")
ilat=fin->lat
ilon=fin->lon

  read_data= asciiread(filename,-1,"string")
  header    = read_data(0)        ; Header. Use for variable names.
  data      = read_data(1:)       ; Get rid of first line which is a header.
  quote = inttochar(34) 
  delim=",:"+quote
  lat1 = tofloat(str_get_field(data,9,delim))
  lon1 = tofloat(str_get_field(data,8,delim))
  I=ind(lon1.ge.min(lonlim).and.lon1.le.max(lonlim).and.lat1.ge.min(latlim).and.lat1.le.max(latlim))
  delete([/lat1,lon1/])

  date = tointeger(str_get_field(data(I),4,delim))
  month= tointeger(mod(floor(date/100.),100))
  time = tointeger(str_get_field(data(I),5,delim))
  lat = tofloat(str_get_field(data(I),9,delim))
  lon= tofloat(str_get_field(data(I),8,delim))
  lon=where(lon.gt.180.,lon-360,lon)
  len=dimsizes(date)  
  tmp=new(len,"float")
  tmp=0
  date2=cd_inv_calendar(floor(date/10000.),tointeger(mod(floor(date/100.),100)),mod(date,100),tointeger(str_get_field(data(I),5,delim)),tmp,tmp,"hours since 1900-01-01 00:00:00",0)
  delete(tmp)

  U10=new((/len,dimsizes(olat),dimsizes(olon)/),"float")
  V10=U10

  I=ind_nearest_coord(lat,ilat,0)
  J=ind_nearest_coord(lon,ilon,0)
  
  FILE2=addfiles(systemfunc (" ls -1 " + DATADIR + "uas/uas_6hrs_ERAI_historical_an-sfc_*.nc"),"r")
  edate=FILE2[:]->time
  FILE3=addfiles(systemfunc (" ls -1 " + DATADIR + "vas/vas_6hrs_ERAI_historical_an-sfc_*.nc"),"r")


; Now, load in the data for each of our time steps.

  do ip=0,len-1
      if(mod(ip,100).eq.0)
         print(date(ip))
      end if
      if(I(ip).gt.30.and.I(ip).lt.210.and.J(ip).gt.30.and.J(ip).lt.450.) then
        K=ind(edate.eq.date2(ip))
        U10(ip,:,:)=FILE2[:]->uas(K,ispan(I(ip)-28,I(ip)+28,1),ispan(J(ip)-28,J(ip)+28,1))
        V10(ip,:,:)=FILE3[:]->vas(K,ispan(I(ip)-28,I(ip)+28,1),ispan(J(ip)-28,J(ip)+28,1))
      end if    
  end do

  U10!0="Time"
  U10!1="lat"
  U10!2="lon"
  U10&lat=olat
  U10&lon=olon
  U10@long_name="Zonal wind"
  U10@units="m/s"
  U10@description="10m zonal wind at observation time within 2000km radius of centre"

  copy_VarMeta(U10,V10)
  V10@long_name="Meridional wind"
  V10@description="10m meridional wind at observation time within 2000km radius of centre"

; Finally, write me a netcdf file

system("/bin/rm -f "+fout)
ncdf = addfile(fout ,"c")  ; open output netCDF file

fAtt               = True            ; assign file attributes
fAtt@title         = "System-centred rainfall composite"
fAtt@source        = "ERAI-Interim data - tas, uas and vas at t=0, 0.75x0.75 grid and track file: "+filename
fAtt@author        = "Acacia Pepler <acacia.pepler@bom.gov.au>"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
filedimdef(ncdf,"time",-1,True)

ncdf->U10=U10
ncdf->V10=V10
