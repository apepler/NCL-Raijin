load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj240_lows_rad2cv1/"
;filename=ECLDIR+"UM_lows_ERAI_proj240_rad2cv1_bigaust_fixes.csv"
;fout=ECLDIR+"moiststable_ERAI_proj240_rad2cv1_aust.nc"

ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/500hPa_z/proj240_lows_rad2cv10_v2/"
filename=ECLDIR+"UM_lows_ERAI_500hPa_proj240_rad2cv10_bigaust_fixes.csv"
fout=ECLDIR+"moiststable_ERAI_500hPa_proj240_rad2cv10_aust.nc"

;ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj240_lows_rad5cv0.15/"
;filename=ECLDIR+"UM_lows_ERAI_proj240_rad2cv1_bigaust_fixes.csv"
;filename=ECLDIR+"ECLfixes_UM_ERAI_proj240_rad5cv0.15.csv"
;fout=ECLDIR+"moiststable_ERAI_proj240_rad5cv0.15.nc"

lev=700
lonlim=(/110,160/)
latlim=(/-45,-10/)

syear=1979
eyear=2017
edir="/g/data1/ub4/erai/netcdf/"
cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/erainterim_reversible_nomix/"

lon0=fspan(-10.5,10.5,29)
lon0!0="lon"
lon0@units="degrees_E"
lat0=fspan(-10.5,10.5,29)
lat0!0="lat"
lat0@units="degrees_N"

fin=addfile(edir+"6hr/atmos/oper_an_pl/v01/z/z_6hrs_ERAI_historical_an-pl_19980701_19980731.nc","r")
elat=fin->lat
elon=fin->lon
elev=fin->lev
ll=ind(elev.eq.lev*100.)
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

  ECL_TCW=new((/len,dimsizes(lat0),dimsizes(lon0)/),"double")
  ECL_TCW!0="Time"
  ECL_TCW!1="lat"
  ECL_TCW!2="lon"
  ECL_TCW&lat=lat0
  ECL_TCW&lon=lon0
  ECL_BV700=ECL_TCW
  ECL_MUCAPE=ECL_TCW

  ECL_TCW@description="6-hour Total column water within 1000km radius of low centre"
  ECL_BV700@description="6-hour 700hPa Brunt-Vaisala frequency within 1000km radius of low centre"
  ECL_MUCAPE@description="6-hour most unstable CAPE within 1000km radius of low centre - calculated reversibly without entrainment from ADowdy files in /g/data/eg3/ajd548/NESP_EWP/CAPE/erainterim_reversible_nomix/. Only available until 2016"
  ECL_TCW@units="kg m**-2"
  ECL_BV700@units="1/s"
  ECL_MUCAPE@units="J/kg"

  I=ind_nearest_coord(lat,elat,0)
  J=ind_nearest_coord(lon,elon,0)

  do ip=0,len-1
    if(mod(ip,100).eq.0)
      print(date(ip))
    end if

    n=(day(ip)-1)*4+time(ip)

    Data=addfile(systemfunc("ls -1 "+edir+"/6hr/atmos/oper_an_sfc/v01/tcw/tcw_6hrs_ERAI_historical_an-sfc_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
    tcw=Data->tcw(n,:,:)
    delete(Data)

    if(year(ip).le.2016) then
     Data=addfile(systemfunc("ls -1 "+cdir+"mucape_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
     cape=Data->mucape(n,:,:)
     delete(Data)
    end if

    ; Need a lot for brunt vasiala

    Data=addfile(systemfunc("ls -1 "+edir+"/6hr/atmos/oper_an_pl/v01/z/z_6hrs_ERAI_historical_an-pl_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
    z=Data->z(n,:,:,:)/9.8065
    Data2=addfile(systemfunc("ls -1 "+edir+"/6hr/atmos/oper_an_pl/v01/ta/ta_6hrs_ERAI_historical_an-pl_"+year(ip)+sprinti("%2.2i",month(ip))+"*"),"r")
    ta=Data2->ta(n,:,:,:)
    th=pot_temp(elev,ta,0,False)
    bv=brunt_vaisala_atm(th,z,0,0)
    bv=where(bv.lt.0,0,bv) 

    delete([/Data,Data2/])
 
    if(I(ip).gt.14.and.J(ip).gt.14.and.I(ip).lt.227.and.J(ip).lt.466) then
      ECL_TCW(ip,:,:)=(/tcw((I(ip)+14):(I(ip)-14),(J(ip)-14):(J(ip)+14))/)
      ECL_BV700(ip,:,:)=(/bv(ll,(I(ip)+14):(I(ip)-14),(J(ip)-14):(J(ip)+14))/)
      if(year(ip).le.2016) then
       ECL_MUCAPE(ip,:,:)=(/cape((I(ip)+14):(I(ip)-14),(J(ip)-14):(J(ip)+14))/)
      end if
    end if

    delete([/tcw,ta,th,bv,z/])
    if(year(ip).le.2016) then
     delete(cape)
    end if
  
  end do

; Finally, write me a netcdf file

system("/bin/rm -f "+fout)
out = addfile(fout,"c")
out->ECL_TCW=ECL_TCW
out->ECL_BV700=ECL_BV700
out->ECL_MUCAPE=ECL_MUCAPE
