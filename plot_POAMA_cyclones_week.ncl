load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

startdate=20130101
keydate1=20130122
keydate2=20130128
; Load POAMA data

f1=addfile("/g/data2/rr8/POAMA2/hindcast/e24/e24a/emn/dac_"+startdate+".nc","r")
lat1=f1->lat
lon1=f1->lon
time1=f1->time
time1a=cd_calendar(time1,2)

x=ind(time1a.ge.keydate1.and.time1a.le.keydate2)
print(x)
slp1=dim_avg_n_Wrap(f1->mslp(x,:,:),0)
prcp1=dim_sum_n_Wrap(f1->hr24_prcp(x,:,:),0)

slp1@lat=lat1
slp1@lon=lon1
prcp1@lat=lat1
prcp1@lon=lon1


; Plotting info 

 res                 = True
; Set some contouring resources.
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn   = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/  1., 5., 10., 25., 50., 100., 200., 300. /)
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

wks = gsn_open_wks("pdf","POAMA_prslp_"+startdate+"_"+keydate1+"_"+keydate2+"_cyclones")
gsn_define_colormap(wks,"prcp_1")                ; choose colormap

res@tiMainString="POAMA"
contour_psl = gsn_csm_contour(wks,slp1,opts_psl)
plot = gsn_csm_contour_map(wks, prcp1, res)
overlay(plot, contour_psl)

pmres=True
pmres@gsMarkerColor = "magenta"
pmres@gsMarkerIndex = 5
pmres@gsMarkerThicknessF = 2
pmres@gsMarkerSizeF = 0.01
dum = new(11,graphic)
do em=0,10
f2="/short/eg3/asp561/cts.dir/gcyc_out/POAMA/proj240_rad5cv0.15/e24a/e"+sprinti("%2.2i",em)+"/tracks_20130101_4month.dat"
data=readAsciiTable(f2,14,"float",0)
I=ind(data(:,6).ge.140.and.data(:,6).le.155.and.data(:, 7).ge.-40.and.data(:,7).le.-10.and.data(:,2).ge.mod(keydate1,1000000).and.data(:,2).le.mod(keydate2,1000000))
print(data(I(0),6)+" "+data(I(0),7))
dum(em)=gsn_add_polymarker(wks,plot,data(I,6),data(I,7),pmres)
delete([/data,I,f2/])
end do


draw(plot)

