load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

yname=2015
ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/BARRA/proj240_lows_rad2cv1/"
filename=ECLDIR+"UM_lows_BARRA_proj240_rad2cv1_bigaust_fixes.csv"
fout=ECLDIR+"ECLcomposite_BARRA_proj240_rad2cv1_aust_"+yname+".nc"

lonlim=(/110,160/)
latlim=(/-45,-10/)

syear=yname
eyear=yname
bdir="/g/data/ma05/BARRA_R/v1/"

lon0=fspan(-10.01,10.01,183)
lon0!0="lon"
lon0@units="degrees_E"
lat0=fspan(-10.01,10.01,183)
lat0!0="lat"
lat0@units="degrees_N"

fin=addfile(bdir+"analysis/spec/mslp/2015/01/mslp-an-spec-PT0H-BARRA_R-v1-20150131T1800Z.nc","r")
elat=fin->latitude
elon=fin->longitude
delete(fin)

  read_data= asciiread(filename,-1,"string")
  header    = read_data(0)        ; Header. Use for variable names.
  data      = read_data(1:)       ; Get rid of first line which is a header.
  quote = inttochar(34) 
  delim=",:"+quote
  date1 = tointeger(str_get_field(data,4,delim))
  lat1 = tofloat(str_get_field(data,9,delim))
  lon1= tofloat(str_get_field(data,8,delim))

  I=ind(date1.ge.syear*10000.and.date1.le.(eyear+1)*10000.and.lat1.ge.min(latlim).and.lat1.le.max(latlim).and.lon1.ge.min(lonlim).and.lon1.le.max(lonlim))
  date = tointeger(str_get_field(data(I),4,delim))
  time = tointeger(str_get_field(data(I),5,delim))
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
  ECL_PRCP=ECL_U10
  ECL_U10@description="6-hour instantaneous 10m zonal wind within 1000km radius of low centre"
  ECL_V10@description="6-hour instantaneous 10m meridional wind within 1000km radius of low centre"
  ECL_SLP@description="6-hour instantaneous slp within 1000km radius of low centre"
  ECL_PRCP@description="6-hour total rainfall within 1000km radius of low centre and +- 3 hours"

  I=ind_nearest_coord(lat,elat,0)
  J=ind_nearest_coord(lon,elon,0)

  do ip=0,len-1
    if(mod(ip,100).eq.0)
      print(date(ip))
    end if

    date3=date(ip)+"T"+sprinti("%2.2i",time(ip))+"00Z"
    Data=addfile(bdir+"analysis/spec/uwnd10m/"+year(ip)+"/"+sprinti("%2.2i",month(ip))+"/uwnd10m-an-spec-PT0H-BARRA_R-v1-"+date3+".nc","r")
    u=Data->uwnd10m
    delete(Data)

    Data=addfile(bdir+"analysis/spec/vwnd10m/"+year(ip)+"/"+sprinti("%2.2i",month(ip))+"/vwnd10m-an-spec-PT0H-BARRA_R-v1-"+date3+".nc","r")
    v=Data->vwnd10m
    delete(Data)

    Data=addfile(bdir+"analysis/spec/mslp/"+year(ip)+"/"+sprinti("%2.2i",month(ip))+"/mslp-an-spec-PT0H-BARRA_R-v1-"+date3+".nc","r")
    p=Data->mslp
    p=p/100.
    delete(Data)

    Data=addfile(bdir+"/forecast/spec_proc/accum_prcp/"+year(ip)+"/"+sprinti("%2.2i",month(ip))+"/accum_prcp-fc-spec_proc-PT1H-BARRA_R-v1-"+date3+".nc","r")
    prcp=Data->accum_prcp(2,:,:)
    delete(Data)

    if(I(ip).gt.91.and.J(ip).gt.91.and.I(ip).lt.670.and.J(ip).lt.1108) then
      ECL_U10(ip,:,:)=(/u((I(ip)-91):(I(ip)+91),(J(ip)-91):(J(ip)+91))/)
      ECL_V10(ip,:,:)=(/v((I(ip)-91):(I(ip)+91),(J(ip)-91):(J(ip)+91))/)
      ECL_SLP(ip,:,:)=(/p((I(ip)-91):(I(ip)+91),(J(ip)-91):(J(ip)+91))/)
      ECL_PRCP(ip,:,:)=(/prcp((I(ip)-91):(I(ip)+91),(J(ip)-91):(J(ip)+91))/)
    end if

    delete([/u,v,p,prcp/])
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
out->ECL_PRCP=ECL_PRCP
