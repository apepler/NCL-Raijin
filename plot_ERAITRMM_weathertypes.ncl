load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;load "./WRFUserARW.ncl"

edir="/g/data1/ub4/erai/netcdf/6hr/atmos/"
tdir="/g/data/eg3/asp561/TRMM/"
odir=""

year=2011
month=12

do day=1,31

date2=(/year,month,day,0/)
date=date2(0)+sprinti("%2.2i",date2(1))+sprinti("%2.2i",date2(2))+sprinti("%2.2i",date2(3))

fdir="/g/data/eg3/asp561/Fronts/"
fname="ERAI_fronts_"

jdir="/g/data/eg3/ajd548/vicci/"
jfname="front_data_JC/cold_fronts_Aus"

cdir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/"
cname="Catto_cyclones_proj100_rad5cv0.15_"
tsname="TS_environments/TS_environments_"

; Set up plotting resources

  res = True
  res@mpFillOn = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is Lowres
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@gsnDraw         = False
  res@gsnFrame        = False
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res@mpGridLineDashPattern = 2 ; lat/lon lines as dashed
  res@mpGridAndLimbOn = False
  res@mpGridSpacingF = 10.
  res@mpProjection = "CylindricalEquidistant"
  res@mpGeophysicalLineThicknessF = 2.0
  res@mpMinLonF = 120
  res@lbOrientation     = "Vertical"
  res@mpMaxLonF = 160
  res@mpMinLatF = -50
  res@mpMaxLatF = -20
  res@mpGridLineThicknessF = 0.5
  res@cnFillOn        = True
  res@cnLinesOn       = False
;  res@lbLabelBarOn        = False
  res@gsnStringFontHeightF = 0.02
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/ 0.1, 0.2, 0.5, 1., 2., 3., 4., 5., 7.5, 10. /)

  res@gsnLeftString = ""
  date3=date2(0)+"-"+date2(1)+"-"+date2(2)+" "+date2(3)+" UTC"
  res@gsnRightString=date3
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

  vres                         = True             ; plot mods desired
  vres@vcRefLengthF            = 0.05            ; define length of vec ref
  vres@vcMinDistanceF          = 0.017
  vres@gsnDraw         = False
  vres@gsnFrame        = False
  vres@vcLineArrowThicknessF = 1.1
  vres@vcRefMagnitudeF           = 10.0             ; define vector ref mag
  vres@tiMainFontHeightF = 0.02
  vres@gsnLeftString = ""
  vres@gsnRightString=""
  vres@vcLineArrowColor="black"
  vres@vcVectorDrawOrder         = "PostDraw"
  vres@vcRefAnnoOrthogonalPosF = -1   ; move ref vector into plot
  vres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  vres@vcRefAnnoString1 = "$VMG$ m/s"
  vres@vcRefAnnoString2On = False

levels=ispan(990,1030,2)
opts_psl = True
     opts_psl@cnLineColor       = "Black"
     opts_psl@cnInfoLabelOn     = False
     opts_psl@cnLineLabelFontHeightF = 0.01
     opts_psl@cnLineLabelPerimOn = False
     opts_psl@gsnContourLineThicknessesScale = 2.5
     opts_psl@gsnDraw         = False
     opts_psl@gsnFrame        = False
     opts_psl@cnLevelSelectionMode = "ExplicitLevels"
     opts_psl@cnLevels             = levels
     opts_psl@cnLineDrawOrder="Draw"
     opts_psl@cnLineLabelPlacementMode = "constant"
     opts_psl@cnLineDashSegLenF = 0.4
     opts_psl@cnLineLabelFontThicknessF = 2.
  opts_psl@gsnLeftString = ""
  opts_psl@gsnRightString=""

  nlevels = dimsizes(levels)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(mod(levels(i),4).ne.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

; Set up panel

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 6
wks_type@wkPaperHeightF = 5

wks = gsn_open_wks(wks_type,odir+"ERAI_TRMM_weather_"+date)
  gsn_define_colormap(wks,"precip_11lev")

; ERAI plotting

  Data=addfile(systemfunc("ls -1 "+edir+"oper_an_sfc/v01/psl/psl_6hrs_ERAI_historical_an-sfc_"+date2(0)+sprinti("%2.2i",date2(1))+"*"),"r")
  time=Data->time
  time2=cd_calendar(time,0)
  I=ind(time2(:,2).eq.date2(2).and.time2(:,3).eq.date2(3))
  
  p=Data->psl(I,:,:)
  p=p/100.
  delete([/Data,I,time,time2/])

; TRMM plotting
  date4=date2(0)+sprinti("%2.2i",date2(1))+sprinti("%2.2i",date2(2))
  if(date2(0).ge.2000.and.date2(0).le.2010) then
    dfile=tdir+date2(0)+"/3B42."+date4+".7A.nc"
    if (fileexists(dfile)) then
      FILE1=addfile(dfile,"r")
      prcp=FILE1->pcp(date2(3)/3,:,:)
    end if
  else if(date2(0).gt.2013) then
      dfile=tdir+date2(0)+"/3B42."+date4+"."+sprinti("%2.2i",date2(3))+".7.nc"
      if (fileexists(dfile)) then
       FILE1=addfile(dfile,"r")
       prcp=FILE1->pcp(0,:,:)
      end if
  else
    dfile=tdir+date2(0)+"/3B42."+date4+".7.nc"
    if (fileexists(dfile)) then
     FILE1=addfile(dfile,"r")
     prcp=FILE1->pcp(date2(3)/3,:,:)
    end if
  end if
  end if
; print(dfile)
  res@gsnLeftString = "TRMM rain rate (mm/h)"
  
  plots = gsn_csm_contour_map(wks,prcp,res)
  contour_psl = gsn_csm_contour(wks,p,opts_psl)
  overlay(plots,contour_psl)

; Now add some weather types

; Thunderstorms first
  res2=True
  res2@gsnDraw  = False                          ; don't draw
  res2@gsnFrame = False
  res2@cnLevels             = 0.9
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 1., 0.5, 0., 0.5 /)/)
  res2@cnFillOn        = True
  res2@cnLinesOn  = False
  res2@gsnLeftString=""
  res2@gsnRightString=""
  res2@cnInfoLabelOn     = False
  res2@lbLabelBarOn = False
  res2@cnLineLabelsOn  = False

a=addfile(jdir+tsname+date2(0)+".nc","r")
time=(a->time)*6
time@units="hours since "+date2(0)+"-1-1 00:00:00"
time2=cd_calendar(time,0)
I1=ind(time2(:,1).eq.date2(1).and.time2(:,2).eq.date2(2).and.time2(:,3).eq.date2(3))
var="Thunderstorm_environments_produced_by_ADowdy_May_2018:"
delete([/time,time2/])

tmp=a->$var$
dowdyTS=tmp(I1,:,:)
dowdyTS&lat@units="degrees north"
dowdyTS&lon@units="degrees east"
delete(tmp)

if(max(dowdyTS).gt.0) then
  contour2 = gsn_csm_contour(wks,dowdyTS,res2)
  overlay(plots,contour2)
end if

; Add Jen's fronts
a=addfile(jdir+jfname+date2(0)+".nc","r")
var="Cold fronts for use in collaboration with A. Dowdy & J. Catto:"

tmp=a->$var$
cattof=tmp(I1,:,:)
cattof&lat@units="degrees north"
cattof&lon@units="degrees east"
delete(tmp)
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 0.5, 0., 1., 1. /)/)

if(max(cattof).gt.0) then
  contour4 = gsn_csm_contour(wks,cattof,res2)
  overlay(plots,contour4)
end if
; Then cyclones

  res3=True
  res3@gsnDraw  = False                          ; don't draw
  res3@gsnFrame = False
  res3@cnLevels             = 0.9
  res3@cnFillOn        = False
  res3@cnLinesOn  = True
  res3@gsnLeftString=""
  res3@gsnRightString=""
  res3@cnInfoLabelOn     = False
  res3@lbLabelBarOn = False
  res3@cnLineLabelsOn  = False
  res3@cnLineColor = "blue"

a=addfile(cdir+cname+date2(0)+".nc","r")
cattoc=a->CYCFLAG(I1,:,:)
cattoc=where(cattoc.eq.1,1,0)

if(max(cattoc).gt.0) then
  contour3 = gsn_csm_contour(wks,cattoc,res3)
  overlay(plots,contour3)
end if

; Add Irina's fronts
; Load in Irina's fronts

data=asciiread(fdir+fname+date2(0)+".dat",-1,"string")
delim=" "

idate=tointeger(str_get_field(data,3,delim))
itime=tointeger(str_get_field(data,4,delim))/100
tmp1=itime
tmp1=0
tmp=cd_inv_calendar(idate/10000,mod(idate/100,100),mod(idate,100),itime,tmp1,tmp1,"hours since 1-1-1 00:00:00",0)
tmp=tmp+6
tmp2=cd_calendar(tmp,0)
I=ind(tmp2(:,0).eq.date2(0).and.tmp2(:,1).eq.date2(1).and.tmp2(:,2).eq.date2(2).and.tmp2(:,3).eq.date2(3))

id = tointeger(str_get_field(data(I),1,delim))
ilat = tofloat(str_get_field(data(I),5,delim))
ilon= tofloat(str_get_field(data(I),6,delim))
delete([/tmp,tmp1,tmp2,I/])

; Add Irina lines to plot
   lres=True
   lres@gsLineColor = "red"
   lres@gsLineThicknessF = 6.0

   ids2=ispan(min(id),max(id),1) ; Possible only for fronts because all of same day together

 dum=new(dimsizes(ids2),graphic)
 do i=0,dimsizes(ids2)-1
  J=ind(id.eq.ids2(i))
  dum(i)=gsn_add_polyline(wks,plots,ilon(J),ilat(J),lres)
  delete(J)
 end do
 draw(plots)
 frame(wks)

delete([/wks,plots,id,ids2,ilat,ilon,data,idate,itime,dum/])

end do
