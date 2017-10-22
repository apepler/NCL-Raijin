load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

startdate=20070601
keydate=20070703

; Load POAMA data

f1=addfile("/g/data2/rr8/POAMA2/hindcast/e24/e24a/emn/dac_"+startdate+".nc","r")
lat1=f1->lat
lon1=f1->lon
time1=f1->time
time1a=cd_calendar(time1,2)

x=ind(time1a.eq.keydate)
print(x)
slp1=f1->mslp(x,:,:)
prcp1=f1->hr24_prcp(x,:,:)
lev=f1->lvl
y=ind(lev.eq.500)
gph1=f1->geop_ht(x,y,:,:) ; 500hPa

slp1@lat=lat1
slp1@lon=lon1
prcp1@lat=lat1
prcp1@lon=lon1
gph1@lat=lat1
gph1@lon=lon1

delete([/time1,time1a,x,y,lev,f1/])

; Load ACCESS data

f1=addfile("/g/data2/ub7/access-s1/hc/raw_model/atmos/psl/daily/emn/da_psl_"+startdate+"_emn.nc","r")
lat2=f1->lat
lon2=f1->lon
time2=f1->time
time2a=cd_calendar(time2,2)

x=ind(time2a.eq.keydate)
print(x)
slp2=f1->psl(x,:,:)/100.

f2=addfile("/g/data2/ub7/access-s1/hc/raw_model/atmos/pr/daily/emn/da_pr_"+startdate+"_emn.nc","r")
prcp2=f2->pr(x,:,:)*60*60*24 ; To convert to Kg/m2/24 hours


; If it exists - currently doesn't for 0601 start times
tmp="/g/data2/ub7/access-s1/hc/raw_model/atmos/zg/daily/emn/da_zg_"+startdate+"_emn.nc"
if(isfilepresent(tmp))
  f3 = addfile(tmp,"r")
  lev=f3->z1_p_level
  y=ind(lev.eq.500)
  gph2=f3->zg(x,y,:,:)
  delete([/f3,lev,y/])

  gph2&lat=where(gph2&lat.gt.90.,90.,gph2&lat)
  print(max(gph2&lat))
else
;  gph2=new(dimsizes(slp2),typeof(gph1),gph1@_FillValue)
end if

slp2!0="lat"
slp2!1="lon"
slp2&lat=lat2
slp2&lon=lon2

prcp2!0="lat"
prcp2!1="lon"
prcp2&lat=lat2
prcp2&lon=lon2
prcp2@long_name=prcp1@long_name
prcp2@units=prcp1@units

;gph2!0="lat"
;gph2!1="lon"
;gph2&lat=lat2
;gph2&lon=lon2

;delete([/time2,time2a,x,f1,f2/])

; Plotting info 

 res                 = True
; Set some contouring resources.
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn   = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/  1., 2., 5., 10., 20., 30., 40., 50. /)
  res@cnFillColors         = (/2,11,10,6,7,8,4,5,17/)
  res@gsnDraw         = False
  res@gsnFrame        = False
; Add map resources
  res@mpFillOn              = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res@mpLimitMode        = "Corners"               ; Portion of map to zoom
  res@mpLeftCornerLatF   = -50.
  res@mpLeftCornerLonF   = 110.
  res@mpRightCornerLatF  = -10.
  res@mpRightCornerLonF  = 180.
  res@mpGridLineDashPattern = 2
  res@mpGridLineThicknessF = 0.5
  res@mpProjection = "CylindricalEquidistant"
; Add label bar resources
  res@lbLabelAutoStride = True
  res@lbBoxMinorExtentF        = 0.13
  res@lbLabelFontHeightF       = 0.012
  res@gsnLeftStringFontHeightF = 0.01
  res@gsnRightStringFontHeightF = 0.008
  res@gsnMaximize     = False    ; Maximize plot in frame
  res@gsnLeftString="24 hour rainfall (mm)"
  res@gsnRightString=""

  opts_psl = True
  opts_psl@cnLineColor       = "Blue"
  opts_psl@cnInfoLabelOn     = False
  opts_psl@cnLineLabelFontHeightF = 0.01
  opts_psl@cnLineLabelPerimOn = False
  opts_psl@gsnContourLineThicknessesScale = 2.0
  opts_psl@gsnDraw         = False
  opts_psl@gsnFrame        = False
  opts_psl@cnLevelSpacingF = 2.0
  opts_psl@gsnLeftString = "Sea level pressure (hPa)"
  opts_psl@gsnRightString = ""

  opts_gph=opts_psl
  opts_gph@cnLineColor       = "Red"
  opts_gph@cnLevelSpacingF = 20.0
  opts_gph@gsnLeftString = ""
  opts_gph@gsnRightString = "500 hPa height (m)"

; Plotting

wks = gsn_open_wks("pdf","ACCESS_POAMA_prslp_"+startdate+"_"+keydate+"_member1")
gsn_define_colormap(wks,"prcp_1")                ; choose colormap
plots = new(2,graphic)

res@tiMainString="POAMA"
contour_psl = gsn_csm_contour(wks,slp1,opts_psl)
plots(0) = gsn_csm_contour_map(wks, prcp1, res)
overlay(plots(0), contour_psl)
;contour_gph=gsn_csm_contour(wks,gph1,opts_gph)
;overlay(plots(0),contour_gph)

res@tiMainString="ACCESS"
contour_psl = gsn_csm_contour(wks,slp2,opts_psl)
plots(1) = gsn_csm_contour_map(wks, prcp2, res)
overlay(plots(1), contour_psl)

tmp="/g/data2/ub7/access-s1/hc/raw_model/atmos/zg/daily/emn/da_zg_"+startdate+"_emn.nc"
if(isfilepresent(tmp))
 contour_gph=gsn_csm_contour(wks,gph2,opts_gph)
 overlay(plots(1),contour_gph)
end if

    res2                   = True              ; modify the panel plot
    res2@gsnFrame          = False             ; don't advance frame yet
    res2@gsnPanelTop       = 0.95
    res2@gsnPanelBottom    = 0.05              ; top 2/3 of frame
    res2@gsnPanelFigureStrings= (/"a)","b)"/)
    res2@amJust   = "TopLeft"
    res2@gsnPaperOrientation = "landscape"
;    res@gsnPanelLabelBar  = True              ; add common colorbar
;    res@lbOrientation     = "vertical"
;    res@pmLabelBarWidthF  = 0.075             ; make thinner
    res2@gsnMaximise=True

    res2@gsnMaximize = False
    gsn_panel(wks,plots,(/1,2/),res2)
    frame(wks)

