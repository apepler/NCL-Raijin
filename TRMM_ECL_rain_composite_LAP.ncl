load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/500hPa_z/proj240_lows_rad2cv10_v2/"
;filename=ECLDIR+"UM_lows_ERAI_500hPa_proj240_rad2cv10_bigaust_fixes.csv"
;fout=ECLDIR+"rain_TRMM_ERAI_500hPa_proj240_rad2cv10_aust_19982015.nc"

ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/"
filename=ECLDIR+"UM_lows_ERAI_proj100_rad5cv0.15_bigaust_fixes.csv"
fout=ECLDIR+"rain_TRMM_ERAI_proj100_rad5cv0.15_aust_19982015.nc"

;ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj240_lows_rad5cv0.15/"
;filename=ECLDIR+"UM_lows_ERAI_proj240_rad2cv1_bigaust_fixes.csv"
;filename=ECLDIR+"ECLfixes_UM_ERAI_proj240_rad5cv0.15.csv"
;fout=ECLDIR+"ECLrain_TRMM_ERAI_proj240_rad5cv0.15_19982015.nc"
DATADIR="/g/data/eg3/asp561/TRMM/"

syear=1998
eyear=2015
lonlim=(/110,160/)
latlim=(/-45,-10/)

namx="_tavg"

; This is the dimensions of the output grid
; Unlike mLDB, where all rounded to nearest 0.5-1 degree so points never match, for this one they can be centred anywhere so centre on the closest point
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
  lon1= tofloat(str_get_field(data,8,delim))
  loc1= tointeger(str_get_field(data,12,delim))

  I=ind(date1.ge.syear*10000.and.date1.le.(eyear+1)*10000.and.lat1.ge.min(latlim).and.lat1.le.max(latlim).and.lon1.ge.min(lonlim).and.lon1.le.max(lonlim))
;  I=ind(date1.ge.syear*10000.and.date1.le.(eyear+1)*10000.and.loc1.eq.1)
  delete([/date1,lat1,lon1,loc1/])

  date = tointeger(str_get_field(data(I),4,delim))
  time = tointeger(str_get_field(data(I),5,delim))/3
  lat = tofloat(str_get_field(data(I),9,delim))
  lon= tofloat(str_get_field(data(I),8,delim))
  len=dimsizes(date)  
  year=tointeger(floor(date/10000.))

  ECL_PRCP=new((/len,dimsizes(lat0),dimsizes(lon0)/),"float")
  ECL_PRCP!0="Time"
  ECL_PRCP!1="lat"
  ECL_PRCP!2="lon"
  ECL_PRCP@description="6-hour instantaneous rainfall rate within 1000km radius of low centre"

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

      if(year(ip).ge.2000.and.year(ip).le.2010) then
         dfile=DATADIR+year(ip)+"/3B42."+date(ip)+".7A.nc"
         if (fileexists(dfile)) then
          FILE1=addfile(dfile,"r")
          p=FILE1->pcp(time(ip),:,:)
         end if
      else if(year(ip).gt.2013) then
         dfile=DATADIR+year(ip)+"/3B42."+date(ip)+"."+sprinti("%2.2i",3*time(ip))+".7.nc"
         if (fileexists(dfile)) then
          FILE1=addfile(dfile,"r")
          p=FILE1->pcp(0,:,:)
         end if
      else
         dfile=DATADIR+year(ip)+"/3B42."+date(ip)+".7.nc"
         if (fileexists(dfile)) then
          FILE1=addfile(dfile,"r")
          p=FILE1->pcp(time(ip),:,:)
         end if
      end if
      end if

      if(I(ip).gt.40.and.J(ip).gt.40.and.I(ip).lt.360.and.J(ip).lt.1400.and.fileexists(dfile)) then
          ECL_PRCP(ip,:,:)=(/p((I(ip)-40):(I(ip)+40),(J(ip)-40):(J(ip)+40))/)
      else if(I(ip).ge.0.and.fileexists(dfile)) then
          x=I(ip)-40
          ECL_PRCP(ip,(40-I(ip)):80,:)=(/p(0:(I(ip)+40),(J(ip)-40):(J(ip)+40))/)
      end if
      end if
      
  end do
  ECL_PRCP!1="lat"
  ECL_PRCP!2="lon"

; Finally, write me a netcdf file

system("/bin/rm -f "+fout)
out = addfile(fout,"c")
out->ECL_PRCP=ECL_PRCP

