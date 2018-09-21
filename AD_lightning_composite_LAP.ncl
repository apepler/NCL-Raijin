load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/500hPa_z/proj240_lows_rad2cv10_v2/"
filename=ECLDIR+"UM_lows_ERAI_500hPa_proj240_rad2cv10_bigaust_fixes.csv"
fout=ECLDIR+"lightning_ERAI_500hPa_proj240_rad2cv10_aust_20052015.nc"

;ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj240_lows_rad5cv0.15/"
;filename=ECLDIR+"UM_lows_ERAI_proj240_rad2cv1_bigaust_fixes.csv"
;filename=ECLDIR+"ECLfixes_UM_ERAI_proj240_rad5cv0.15.csv"
;fout=ECLDIR+"ECLlightning_ERAI_proj240_rad5cv0.15_20052015.nc"

DATADIR="/g/data/eg3/asp561/CattoData_20052015/new_lightning_Australasia/"

syear=2005
eyear=2015
lonlim=(/110,160/)
latlim=(/-45,-10/)

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
  time = tointeger(str_get_field(data(I),5,delim))
  lat = tofloat(str_get_field(data(I),9,delim))
  lon= tofloat(str_get_field(data(I),8,delim))
  len=dimsizes(date)  
  year=tointeger(floor(date/10000.))

  ECL_lightning=new((/len,dimsizes(lat0),dimsizes(lon0)/),"float")
  ECL_lightning!0="Time"
  ECL_lightning!1="lat"
  ECL_lightning!2="lon"
  ECL_lightning@description="6-hour lightning (+- 3 hour accumulation) within 1000km radius of low centre"

; Get info about the 
  f1=addfile(DATADIR+"lightning_Australasia0.250000degree_6.00000hr_2005.nc","r")
  lat1=f1->lat ; -44.5 to 15.5
  lon1=f1->lon ; 112-215

  I=ind_nearest_coord(lat,lat1,0)
  J=ind_nearest_coord(lon,lon1,0)

; Now, load in the data for each of our time steps.
  do yy=syear,eyear
    Y=ind(year.eq.yy)
    len2=dimsizes(Y)

    dfile=DATADIR+"lightning_Australasia0.250000degree_6.00000hr_"+yy+".nc"
    FILE1=addfile(dfile,"r")
    ltime1=FILE1->time
    ltime1@units="hours since "+yy+"-01-01 00:00:00"
    ltime2=cd_calendar(ltime1,0)
    ldate=tointeger(ltime2(:,0))*10000+tointeger(ltime2(:,1))*100+tointeger(ltime2(:,2))
    ltime=tointeger(ltime2(:,3))

    do ip=0,len2-1
      if(mod(ip,100).eq.0)
         print(date(Y(ip)))
      end if

      ip2=ind(ldate.eq.date(Y(ip)).and.ltime.eq.time(Y(ip)))
      p=FILE1->Lightning_observed(ip2,:,:)

      if(I(Y(ip)).gt.40.and.J(Y(ip)).gt.40.and.I(Y(ip)).lt.200.and.J(Y(ip)).lt.320.and.fileexists(dfile)) then
          ECL_lightning(Y(ip),:,:)=(/p((I(Y(ip))-40):(I(Y(ip))+40),(J(Y(ip))-40):(J(Y(ip))+40))/)
;      print(date(Y(ip))+" "+ldate(ip2))
;      else if(I(Y(ip)).gt.0.and.J(Y(ip)).gt.40.and.I(Y(ip)).lt.200.and.J(Y(ip)).lt.320.and.fileexists(dfile)) then
;          x=I(Y(ip))-40
;          ECL_lightning(Y(ip),(40-I(Y(ip))):80,:)=(/p(0:(I(Y(ip))+40),(J(Y(ip))-40):(J(Y(ip))+40))/)
;          print(date(Y(ip))+" "+ldate(ip2))
;      end if
      end if
    end do
   delete([/Y,len2,dfile,FILE1,ltime1,ltime2,ltime,ldate/])
      
  end do
  ECL_lightning!1="lat"
  ECL_lightning!2="lon"

; Finally, write me a netcdf file

system("/bin/rm -f "+fout)
out = addfile(fout,"c")
out->ECL_lightning=ECL_lightning

