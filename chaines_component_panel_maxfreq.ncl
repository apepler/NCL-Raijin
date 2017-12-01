load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

levels=ispan(0,20,2)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon

data=new((/3,dimsizes(time),dimsizes(lat),dimsizes(lon)/),"float")
data(0,:,:,:)=(/ifile->CA/)
data(1,:,:,:)=(/ifile->CB/)
data(2,:,:,:)=(/ifile->CHaines/)

years=ispan(1979,2016,1)

CHp=new((/3,dimsizes(lat),dimsizes(lon)/),"float")

do x=0,2
  CHp(x,:,:)=(/dim_max_n(data(x,:,:,:),0)/)
end do

a=dimsizes(CHp)
CH95p=new((/3,dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!3="lat"
CH95p&lat=lat
CH95p!4="lon"
CH95p&lon=lon

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  do x=0,2
    CH95p(x,y,m,:,:)=dim_num_n(data(x,I,:,:).eq.conform_dims(dimsizes(data(x,I,:,:)),CHp(x,:,:),(/1,2/)),0)
  end do
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  do x=0,2
    CH95p(x,y,0,:,:)=dim_num_n(data(x,I,:,:).eq.conform_dims(dimsizes(data(x,I,:,:)),CHp(x,:,:),(/1,2/)),0)
  end do
  delete(I)
end do

CH95pave=dim_avg_n_Wrap(CH95p,0)

  wks = gsn_open_wks("pdf","chaines_components_eqmax")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels
  res@gsnAddCyclic = False
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpFillOn              = True
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res@mpLimitMode        = "Corners"               ; Portion of map to zoom
  res@mpLeftCornerLatF   = -45.
  res@mpLeftCornerLonF   = 110.
  res@mpRightCornerLatF  = -10.
  res@mpRightCornerLonF  = 155.
  res@mpGridLineDashPattern = 2
  res@mpGridLineThicknessF = 0.5
  res@mpProjection = "CylindricalEquidistant"
;  res@cnFillMode = "RasterFill"
;  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False
;  res@gsnMaximize=True
  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "Draw"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"
  res@gsnStringFontHeightF = 0.035
  res@gsnRightString=""
  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02    
;  res@cnLabelBarEndStyle="IncludeMinMaxLabels"
  plots1 = new(4,graphic)
  plots2 = new(4,graphic)
  plots3 = new(4,graphic)

res2=res
res2@tmYLLabelFontColor = "white"

labs1=(/"a","b","c","d"/)
labs2=(/"e","f","g","h"/)
labs3=(/"i","j","k","l"/)

  do m=0,3

;  Make the array as we go
  tmp=dim_avg_n_Wrap(CH95p(0,:,m,:,:),0)
  tmp!0="lat"
  tmp&lat=lat
  tmp!1="lon"
  tmp&lon=lon

  if(m.lt.3) then
   res@tmXBLabelFontColor = "white"
   res2@tmXBLabelFontColor = "white"
  else
   res@tmXBLabelFontColor = "black"
   res2@tmXBLabelFontColor = "black"
  end if

  res@gsnLeftString=labs1(m)+") "+snames(m)
  plots1(m) = gsn_csm_contour_map(wks,tmp, res)     ; create plot

  tmp=(/dim_avg_n_Wrap(CH95p(1,:,m,:,:),0)/)
  res2@gsnLeftString=labs2(m)+") "+snames(m)
  plots2(m) = gsn_csm_contour_map(wks,tmp, res2)     ; create plot

  tmp=(/dim_avg_n_Wrap(CH95p(2,:,m,:,:),0)/)

  res2@gsnLeftString=labs3(m)+") "+snames(m)
  plots3(m) = gsn_csm_contour_map(wks,tmp, res2)     ; create plot
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize=True
  resP@gsnPaperOrientation = "portrait"
  
  resP@lbTitleString="CA == max (days)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
;  resP@cnLabelBarEndStyle="IncludeMinMaxLabels"
;  resP@lbLabelStrings=levels1a
  resP@pmLabelBarOrthogonalPosF = -0.02

  resP@gsnPanelLabelBar    = True
;  resP@gsnPanelYWhiteSpacePercent = 3
  resP@amJust   = "TopLeft"
  resP@gsnPanelLeft=0.15
  resP@gsnPanelRight=0.4
  resP@gsnPanelBottom=0.05
  resP@gsnFrame = False ; don't advance frame yet
  resP@gsnPanelMainPosXF  = 0.28
  resP@gsnPanelFigureStringsFontHeightF = 0.012
  resP@gsnPanelMainFontHeightF = 0.018
  resP@pmLabelBarWidthF = 0.2
  resP@lbLabelFontHeightF = 0.012
  resP@lbTitleFontHeightF = 0.012
  gsn_panel(wks,plots1,(/4,1/),resP)

  resP2=resP
  resP2@lbTitleString="CB == max (days)"
;  delete(resP2@lbLabelStrings)
;  resP2@lbLabelStrings=levels2a
  resP2@gsnPanelLeft=0.4
  resP2@gsnPanelRight=0.65
  resP2@gsnPanelMainPosXF  = 0.52
 gsn_panel(wks,plots2,(/4,1/),resP2)

  resP2=resP

  resP2@lbTitleString= "CHaines == max (days)"
  resP2@gsnPanelLeft=0.65
  resP2@gsnPanelRight=0.9
  resP2@gsnPanelMainPosXF  = 0.77
 gsn_panel(wks,plots3,(/4,1/),resP2)

frame(wks)


