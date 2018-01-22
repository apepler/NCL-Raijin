load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

ECLDIR="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/"
;DATADIR="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_sfc/v01/psl/"
DATADIR="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_pl/v01/z/"

; Set up my lat-lon grids as 15x15 degrees

  olon=fspan(-15,15,41)
  olon!0="lon"
  olon@units="degrees_E"
  olat=fspan(-15,15,41)
  olat!0="lat"
  olat@units="degrees_N"

fin=addfile(DATADIR+"z_6hrs_ERAI_historical_an-pl_20160901_20160930.nc","r")
;  fin=addfile(DATADIR+"psl_6hrs_ERAI_historical_an-sfc_20160901_20160930.nc","r")
  ilev=fin->lev
  L=ind(ilev.eq.50000)
  ilat=fin->lat
  ilon=fin->lon
  delete(fin)

do year=2009,2016 
  print(year)
  read_data= asciiread(ECLDIR+"ERAI_UM_highs_ERAI_proj100_rad10cv0.075_500km_"+year+".csv",-1,"string")
  header    = read_data(0)        ; Header. Use for variable names.
  data      = read_data(1:)       ; Get rid of first line which is a header.
  quote = inttochar(34)
  delim=",:"+quote
  date1 = tointeger(str_get_field(data,4,delim))
  len=dimsizes(date1)
  lat = tofloat(str_get_field(data,9,delim))
  lon= tofloat(str_get_field(data,8,delim))  ; To make same as ERAI
  lon=where(lon.gt.180.,lon-360,lon) ; To make same as ERAI

  tmp=new(len,"float")
  tmp=0
  date2=cd_inv_calendar(floor(date1/10000.),tointeger(mod(floor(date1/100.),100)),mod(date1,100),tointeger(str_get_field(data,5,delim)),tmp,tmp,"hours since 1900-01-01 00:00:00",0)

; Output files

  SLP_composite=new((/len,dimsizes(olat),dimsizes(olat)/),"double")

  I=ind_nearest_coord(lat,ilat,0)
  J=ind_nearest_coord(lon,ilon,0)

FILE1=addfiles(systemfunc (" ls -1 " + DATADIR + "z_6hrs_ERAI_historical_an-pl_" + year + "*.nc"),"r")
;  FILE1=addfiles(systemfunc (" ls -1 " + DATADIR + "psl_6hrs_ERAI_historical_an-sfc_" + year + "*.nc"),"r")
  edate=FILE1[:]->time
;  slp=FILE1[:]->psl/100
   slp=FILE1[:]->z(:,L,:,:)/9.80665
 
  do ip=0,len-1
;    if(mod(ip,1000).eq.0)
;       print(date1(ip))
;    end if
    K=ind(edate.eq.date2(ip))

    if(I(ip).ge.21.and.I(ip).lt.220) then ; Don't bother if too far S or N
      if(J(ip).lt.20) then
      ; First half comes from the right hand side of the matrix

        SLP_composite(ip,:,ispan(0,19-J(ip),1))=slp(K,ispan(I(ip)-20,I(ip)+20,1),ispan(460+J(ip),479,1))
        SLP_composite(ip,:,ispan(20-J(ip),40,1))=slp(K,ispan(I(ip)-20,I(ip)+20,1),ispan(0,J(ip)+20,1))

      else if(J(ip).ge.460) then
      ; Second half comes from left hand size od the matrix

        SLP_composite(ip,:,ispan(0,499-J(ip),1))=slp(K,ispan(I(ip)-20,I(ip)+20,1),ispan(J(ip)-20,479,1))
        SLP_composite(ip,:,ispan(500-J(ip),40,1))=slp(K,ispan(I(ip)-20,I(ip)+20,1),ispan(0,J(ip)-460,1))

      else
       SLP_composite(ip,:,:)=slp(K,ispan(I(ip)-20,I(ip)+20,1),ispan(J(ip)-20,J(ip)+20,1))       
      end if
      end if
    end if
  end do
 SLP_composite!0="time"
 SLP_composite&time=date2
 SLP_composite!1="lat"
 SLP_composite&lat=olat
 SLP_composite!2="lon"
 SLP_composite&lon=olon


 SLP_composite@description="6-hour instantaneous 500 hPa geopotential height within 15 degree radius of system centre"
 SLP_composite@units="m"

 fout=ECLDIR+"ERAI_UM_highs_ERAI_proj100_rad10cv0.075_500km_"+year+"_slp.nc"
; system("/bin/rm -f "+fout)
 out = addfile(fout,"w")
 out->z500=SLP_composite
 delete([/SLP_composite,I,J,K,date1,date2,edate,header,data,slp,lat,lon,tmp,FILE1,read_data/])

end do

