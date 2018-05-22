load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

DATADIR="/g/data/eg3/asp561/TRMM/"
syear=1998

dir="/short/eg3/asp561/cts.dir/gcyc_out/"

reanal="ERAI"
proj="proj100"
thresh="rad10cv0.075"
type="highs"
edir=dir+reanal+"/"+proj+"_"+type+"_"+thresh+"/"
filename=edir+"UM_"+type+"_"+reanal+"_"+proj+"_"+thresh+"_bigaust_fixes.csv"
fout=edir+"UM_"+type+"_"+reanal+"_"+proj+"_"+thresh+"_austrain.nc"

latlim=(/-45,-10/)
lonlim=(/110,155/)

; This is the dimensions of the output grid
  lon0=fspan(-10,10,81)
  lon0!0="lon"
  lon0@units="degrees_E"
  lat0=fspan(-10,10,81)
  lat0!0="lat"
  lat0@units="degrees_N"

  read_data= asciiread(filename,-1,"string")
  header    = read_data(0)        ; Header. Use for variable names.
  data      = read_data(1:)       ; Get rid of first line which is a header.
  quote = inttochar(34) 
  delim=",:"+quote
  date1 = tointeger(str_get_field(data,4,delim))
  lat1 = tofloat(str_get_field(data,9,delim))
  lon1 = tofloat(str_get_field(data,8,delim))
  I=ind(date1.ge.syear*10000.and.lon1.ge.min(lonlim).and.lon1.le.max(lonlim).and.lat1.ge.min(latlim).and.lat1.le.max(latlim))
  delete([/date1,lat1,lon1/])

  date = tointeger(str_get_field(data(I),4,delim))
  time = tointeger(str_get_field(data(I),5,delim))/3
  lat = tofloat(str_get_field(data(I),9,delim))
  lon= tofloat(str_get_field(data(I),8,delim))
  len=dimsizes(date)  
  year=tointeger(floor(date/10000.))

  PRCP=new((/len,dimsizes(lat0),dimsizes(lon0)/),"float")
  PRCP!0="Time"
  PRCP!1="lat"
  PRCP!2="lon"
  PRCP@description="6-hour instantaneous rainfall rate within 2000km radius of centre"

; Get info about the 
  f1=addfile(DATADIR+"1998/3B42.19981003.7.nc","r")
  lat1=f1->latitude ; +-50
  lon1=f1->longitude ; -180  to 180, need to reshape? Naaa, don't care about 170E lows.

  lon1=where(lon1.lt.0,lon1+360.,lon1)

  I=ind_nearest_coord(lat,lat1,0)
  J=ind_nearest_coord(lon,lon1,0)

; Now, load in the data for each of our time steps.

  do ip=0,len-1
      if(mod(ip,100).eq.0)
         print(date(ip))
      end if
      if(I(ip).gt.40.and.I(ip).lt.360.and.J(ip).gt.20.and.J(ip).lt.1400.) then

      if(year(ip).ge.2000.and.year(ip).le.2010) then
         FILE1=addfile(DATADIR+year(ip)+"/3B42."+date(ip)+".7A.nc","r")
      else
         FILE1=addfile(DATADIR+year(ip)+"/3B42."+date(ip)+".7.nc","r")
      end if

      p=FILE1->pcp(time(ip),:,:)   
      if(I(ip).gt.40.and.J(ip).gt.40.and.I(ip).lt.360.and.J(ip).lt.1400.) then
          PRCP(ip,:,:)=(/p((I(ip)-40):(I(ip)+40),(J(ip)-40):(J(ip)+40))/)
      end if
      
      delete([/FILE1,p/])
      end if
  end do
  PRCP!1="lat"
  PRCP!2="lon"

; Finally, write me a netcdf file

system("/bin/rm -f "+fout)
out = addfile(fout,"c")
out->ECL_PRCP=ECL_PRCP

