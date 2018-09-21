load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/500hPa_z/proj240_lows_rad2cv10_v2/"
;filename=ECLDIR+"UM_lows_ERAI_500hPa_proj240_rad2cv10_bigaust_fixes.csv"
;fout=ECLDIR+"composite_ERAI_500hPa_proj240_rad2cv10_aust.nc"

ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/"
filename=ECLDIR+"UM_lows_ERAI_proj100_rad5cv0.15_bigaust_fixes.csv"
;filename=ECLDIR+"ECLfixes_UM_ERAI_proj240_rad5cv0.15.csv"
fout=ECLDIR+"composite_ERAI_proj100_rad5cv0.15_aust.nc"

lonlim=(/110,160/)
latlim=(/-45,-10/)

syear=1979
eyear=2017
edir="/g/data1/ub4/erai/netcdf/"

lon0=fspan(-10.5,10.5,29)
lon0!0="lon"
lon0@units="degrees_E"
lat0=fspan(-10.5,10.5,29)
lat0!0="lat"
lat0@units="degrees_N"

fin=addfile(edir+"6hr/atmos/oper_an_sfc/v01/uas/uas_6hrs_ERAI_historical_an-sfc_19980701_19980731.nc","r")
elat=fin->lat
elon=fin->lon
delete(fin)

  read_data= asciiread(filename,-1,"string")
  header    = read_data(0)        ; Header. Use for variable names.
  data      = read_data(1:)       ; Get rid of first line which is a header.
  quote = inttochar(34) 
  delim=",:"+quote
  date1 = tointeger(str_get_field(data,4,delim))
  lat1 = tofloat(str_get_field(data,9,delim))
  lon1= tofloat(str_get_field(data,8,delim))
  loc1= tointeger(str_get_field(data,12,delim))

  I=ind(date1.ge.syear*10000.and.date1.le.(eyear+1)*10000.and.lat1.ge.min(latlim).and.lat1.le.max(latlim).and.lon1.ge.min(lonlim).and.lon1.le.max(lonlim))
;  I=ind(date1.ge.syear*10000.and.date1.le.(eyear+1)*10000)
;    I=ind(date1.ge.syear*10000.and.date1.le.(eyear+1)*10000.and.loc1.eq.1)
  delete([/date1,lat1,lon1,loc1/])

  date = tointeger(str_get_field(data(I),4,delim))
  time = tointeger(str_get_field(data(I),5,delim))/6
  lat = tofloat(str_get_field(data(I),9,delim))
  lon= tofloat(str_get_field(data(I),8,delim))
  len=dimsizes(date)  
  year=floor(date/10000.)
  month=tointeger(mod(floor(date/100.),100))
  day=mod(date,100)
  delete(I)

; Output files

  ECL_U10=new((/len,dimsizes(lat0),dimsizes(lon0)/),"double")
  ECL_U10!0="Time"
  ECL_U10!1="lat"
  ECL_U10!2="lon"
  ECL_U10&lat=lat0
  ECL_U10&lon=lon0
  ECL_V10=ECL_U10
  ECL_SLP=ECL_U10
;  ECL_PRCP=ECL_U10
  ECL_U10@description="6-hour instantaneous 10m zonal wind within 1000km radius of low centre"
  ECL_V10@description="6-hour instantaneous 10m meridional wind within 1000km radius of low centre"
  ECL_SLP@description="6-hour instantaneous slp within 1000km radius of low centre"
;  ECL_PRCP@description="6-hour total rainfall within 1000km radius of low centre and +- 3 hours"

  I=ind_nearest_coord(lat,elat,0)
  J=ind_nearest_coord(lon,elon,0)

  do ip=0,len-1
    if(mod(ip,100).eq.0)
      print(date(ip))
    end if

    n=(day(ip)-1)*4+time(ip)

    Data=addfile(systemfunc("ls -1 "+edir+"/6hr/atmos/oper_an_sfc/v01/uas/uas_6hrs_ERAI_historical_an-sfc_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
    u=Data->uas(n,:,:)
    delete(Data)

    Data=addfile(systemfunc("ls -1 "+edir+"/6hr/atmos/oper_an_sfc/v01/vas/vas_6hrs_ERAI_historical_an-sfc_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
    v=Data->vas(n,:,:)
    delete(Data)

    Data=addfile(systemfunc("ls -1 "+edir+"/6hr/atmos/oper_an_sfc/v01/psl/psl_6hrs_ERAI_historical_an-sfc_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
    p=Data->psl(n,:,:)
    p=p/100. 
    delete(Data)

;    Data=addfile(systemfunc("ls -1 "+edir+"/3hr/atmos/oper_fc_sfc/v01/tp/tp_3hrs_ERAI_historical_fc-sfc_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
;    n=(day(ip)-1)*8+time(ip)
;    prcp=dim_sum_n_Wrap(Data->tp(n:(n+1),:,:),0)*1000.
;    delete(Data)
 
    if(I(ip).gt.14.and.J(ip).gt.14.and.I(ip).lt.227.and.J(ip).lt.466) then
      ECL_U10(ip,:,:)=(/u((I(ip)+14):(I(ip)-14),(J(ip)-14):(J(ip)+14))/)
      ECL_V10(ip,:,:)=(/v((I(ip)+14):(I(ip)-14),(J(ip)-14):(J(ip)+14))/)
      ECL_SLP(ip,:,:)=(/p((I(ip)+14):(I(ip)-14),(J(ip)-14):(J(ip)+14))/)
;      ECL_PRCP(ip,:,:)=(/prcp((I(ip)+14):(I(ip)-14),(J(ip)-14):(J(ip)+14))/)
    end if

;    delete([/u,v,p,prcp/])
    delete([/u,v,p/])
  end do

; Finally, write me a netcdf file
ECL_WS10=(ECL_U10^2+ECL_V10^2)^0.5
print(max(ECL_WS10))
copy_VarMeta(ECL_U10,ECL_WS10)
ECL_WS10@description="6-hour instantaneous wind speed within 1000km radius of low centre"

system("/bin/rm -f "+fout)
out = addfile(fout,"c")
out->ECL_U10=ECL_U10
out->ECL_V10=ECL_V10
out->ECL_WS10=ECL_WS10
out->ECL_SLP=ECL_SLP
;out->ECL_PRCP=ECL_PRCP
