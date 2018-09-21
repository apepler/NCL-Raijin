load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

year=2010
daylist=(/0904/)     ;(/0215,0531,0710/) ; Dates I want
;numlist=ispan(0,1400,100)

a=addfile("cold_fronts_Aus"+year+".nc","r")

; Getting file with spaces is a pain
var="Cold fronts for use in collaboration with A. Dowdy & J. Catto:"
catto=a->$var$
catto&lat@units="degrees north"
catto&lon@units="degrees east"

time=(a->time)*6
time@units="hours since "+year+"-1-1 00:00:00"
time2=cd_calendar(time,0)
jdate=100*time2(:,1)+time2(:,2)

; Load in Irina's fronts

data=asciiread("ERAI_fronts_"+year+".dat",-1,"string")
delim=" "
tmp = 20000000+tointeger(str_get_field(data,3,delim))
I=ind((tmp/10000).eq.year)
idate=mod(tmp(I),10000)
itime=tointeger(str_get_field(data(I),4,delim))
id = tointeger(str_get_field(data(I),1,delim))
ilat = tofloat(str_get_field(data(I),5,delim))
ilon= tofloat(str_get_field(data(I),6,delim))
delete(I)

; Set up code to make a plot of Catto fronts for selected day
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
   res@cnLineColor = "blue"
   res@gsnAddCyclic = False
   res@cnFillOn        = True               ; color plot desired
  res@lbLabelBarOn        = False
   res@cnLineLabelsOn  = False              ; turn off contour lines
  res@cnLevelSelectionMode = "ExplicitLevels"

  res@cnLevels             = (/0.9,1.2/)
  res@cnFillPalette = (/"white","blue","white"/)
  res@mpFillOn              = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = True                ; Turn on lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks   
  res@tmXTOn = False
  res@tmYROn = False
   res@mpLimitMode       = "Corners"            ; choose range of map
   res@mpLeftCornerLatF  = -45.
   res@mpLeftCornerLonF  = 85.
   res@mpRightCornerLatF = -5.
   res@mpRightCornerLonF = 170.
;  res@cnFillMode = "RasterFill"   
   res@mpProjection        = "LambertConformal"
   res@mpLambertParallel1F = -18
   res@mpLambertParallel2F = -36
   res@mpLambertMeridianF  = 136
   res@gsnLeftString="Cold Fronts"

   lres=True
   lres@gsLineColor = "red"
   lres@gsLineThicknessF = 3.0

  txres               = True
  txres@txFontHeightF = 0.02

  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.


;do p=0,dimsizes(daylist)-1
;   n=ind(jdate.eq.daylist(p).and.time2(:,3).eq.0)
;   date=daylist(p)

do p=0,dimsizes(numlist)-1
   n=numlist(p)
   date=time2(n,0)+"-"+time2(n,1)+"-"+time2(n,2)+"_"+time2(n,3)+"UTC"
   wks = gsn_open_wks(wks_type,"FrontComp_"+date)          ; send graphics to PNG file
   res@gsnRightString=date
   plot=gsn_csm_contour_map(wks,catto(n,:,:), res)

; Add Irina lines to plot
   I=ind(idate.eq.jdate(n).and.itime.eq.time2(n,3))
   ids2=ispan(min(id(I)),max(id(I)),1) ; Possible only for fronts because all of same day together

 dum=new(dimsizes(ids2),graphic)
 do i=0,dimsizes(ids2)-1
  J=ind(id.eq.ids2(i))
  dum(i)=gsn_add_polyline(wks,plot,ilon(J),ilat(J),lres)
  delete(J)
 end do
 draw(plot)

; Add legend to plot somehow?

  gsres@gsMarkerColor = "blue"
  gsn_polymarker_ndc(wks,0.23,0.15,gsres)
  gsn_text_ndc      (wks,"Catto",0.3,0.15,txres)

  gsres@gsMarkerColor = "red"
  gsn_polymarker_ndc(wks,0.4,0.15,gsres)          
  gsn_text_ndc      (wks,"Simmonds",0.5,0.15,txres)

; draw(plot)
  frame(wks)
  delete([/wks,plot,dum,I,ids2/])
end do
