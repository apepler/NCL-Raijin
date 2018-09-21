load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

edir="/g/data1/ub4/erai/netcdf/6hr/atmos/"
odir=""
date2=(/2012,10,11,0/)
date=date2(0)+sprinti("%2.2i",date2(1))+sprinti("%2.2i",date2(2))+sprinti("%2.2i",date2(3))

cdir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/"
cname="austlows_5deg_proj100_rad5cv0.15_"

fdir="/g/data/eg3/asp561/Fronts/"
fname="ERAI_fronts_"

jdir="/g/data/eg3/ajd548/vicci/"
jfname="front_data_JC/cold_fronts_Aus"
jcname="cyclone_data_JC/cyclones_"

; Jen's fronts

a=addfile(jdir+jfname+date2(0)+".nc","r")
var="Cold fronts for use in collaboration with A. Dowdy & J. Catto:"
time=(a->time)*6
time@units="hours since "+date2(0)+"-1-1 00:00:00"
time2=cd_calendar(time,0)
I1=ind(time2(:,1).eq.date2(1).and.time2(:,2).eq.date2(2).and.time2(:,3).eq.date2(3))

tmp=a->$var$
cattof=tmp(I1,:,:)
cattof&lat@units="degrees north"
cattof&lon@units="degrees east"


delete(tmp)

; Jen's cyclones

a=addfile(jdir+jcname+date2(0)+"_0.75.nc","r")
cattoc=a->IDCLUST(I1,:,:)
cattoc=where(cattoc.gt.0,1,0)
;copy_VarMeta(cattof,cattoc)

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

; Set up code to make a plot of fronts & cyclones for selected day
; In the lambert conformal projection

  wks_type = "pdf"
  wks_type@wkPaperWidthF = 6
  wks_type@wkPaperHeightF = 4

   res                 = True                    ; plot mods desired
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@gsnPaperOrientation = "portrait"
  res@gsnMaximize      = True
   res@cnLinesOn       = True
   res@gsnAddCyclic = False
   res@cnFillOn        = False               ; color plot desired
  res@lbLabelBarOn        = False
  res@mpFillOn              = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = True                ; Turn on lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res@tmXTOn = False
  res@tmYROn = False
  res@mpGridSpacingF = 10.
  res@mpProjection = "CylindricalEquidistant"
  res@mpGeophysicalLineThicknessF = 2.0
  res@mpMinLonF = 110
  res@lbOrientation     = "Vertical"
  res@mpMaxLonF = 160
  res@mpMinLatF = -50
  res@mpMaxLatF = -10
  res@mpGridLineThicknessF = 0.5
  res@lbLabelBarOn = False

levels=ispan(990,1030,4)
  res@cnLineColor       = "Black"
  res@cnInfoLabelOn     = False
  res@cnLineLabelFontHeightF = 0.01
  res@cnLineLabelPerimOn = False
  res@gsnContourLineThicknessesScale = 2.5
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels
  res@cnLineDrawOrder="Draw"
  res@cnLineLabelPlacementMode = "constant"
  res@cnLineDashSegLenF = 0.4
  res@cnLineLabelFontThicknessF = 2.
  res@gsnLeftString = "ERAI"

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

  res2=True
  res2@gsnDraw  = False                          ; don't draw
  res2@gsnFrame = False
  res2@cnLevels             = 0.9
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 1., 0., 0., 0.25 /)/)
  res2@cnFillOn        = True
  res2@cnLinesOn  = False
  res2@gsnLeftString=""
  res2@gsnRightString=""
  res2@cnInfoLabelOn     = False
  res2@lbLabelBarOn = False
  res2@cnLineLabelsOn  = False

   lres=True
   lres@gsLineColor = "red"
   lres@gsLineThicknessF = 3.0

  txres               = True
  txres@txFontHeightF = 0.02

  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.

  date3=date2(0)+"-"+date2(1)+"-"+date2(2)+" "+date2(3)+" UTC"
   wks = gsn_open_wks(wks_type,"FrontsCyclones_"+date+"_v2")          ; send graphics to PNG file

   res@gsnRightString=date3
  Data=addfile(systemfunc("ls -1 "+edir+"oper_an_sfc/v01/psl/psl_6hrs_ERAI_historical_an-sfc_"+date2(0)+sprinti("%2.2i",date2(1))+"*"),"r")
   time3=cd_calendar(Data->time,0)
   I2=ind(time3(:,2).eq.date2(2).and.time3(:,3).eq.date2(3))
  p=Data->psl(I2,:,:)
  p=p/100.
  plots = gsn_csm_contour_map(wks,p,res)

  FILES1=addfile(cdir+cname+date2(0)+".nc","r")
  systems=FILES1->systems(I1,:,:)
  contour2 = gsn_csm_contour(wks,systems,res2)
  overlay(plots,contour2)

; Add Jen's cyclones
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 0., 0., 1., 0.25 /)/)
  contour3 = gsn_csm_contour(wks,cattoc,res2)
  overlay(plots,contour3)

; Add Jen's fronts
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 0., 0., 1., 1. /)/)
  contour4 = gsn_csm_contour(wks,cattof,res2)
  overlay(plots,contour4)

; Add Irina's fronts

; Add Irina lines to plot
   ids2=ispan(min(id),max(id),1) ; Possible only for fronts because all of same day together

 dum=new(dimsizes(ids2),graphic)
 do i=0,dimsizes(ids2)-1
  J=ind(id.eq.ids2(i))
  dum(i)=gsn_add_polyline(wks,plots,ilon(J),ilat(J),lres)
  delete(J)
 end do
 draw(plots)

; Add legend to plot somehow?

  gsres@gsMarkerColor = "blue"
  gsn_polymarker_ndc(wks,0.23,0.05,gsres)
  gsn_text_ndc      (wks,"Catto",0.3,0.05,txres)

  gsres@gsMarkerColor = "red"
  gsn_polymarker_ndc(wks,0.4,0.05,gsres)          
  gsn_text_ndc      (wks,"Simmonds",0.5,0.05,txres)

  frame(wks)
  delete([/wks,plots/])
