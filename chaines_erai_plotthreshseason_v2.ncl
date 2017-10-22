load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

thresh1=8
levels1=ispan(0,80,10)
thresh2=90
levels2=ispan(0,30,3)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)

CH=ifile->CHaines
CHp=ifile->CH_percentile

years=ispan(1979,2016,1)

a=dimsizes(CHp)
CH95p=new((/2,dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!3="lat"
CH95p&lat=CHp&lat
CH95p!4="lon"
CH95p&lon=CHp&lon

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(0,y,m,:,:)=dim_num_n(CH(I,:,:).ge.thresh1,0)
  CH95p(1,y,m,:,:)=dim_num_n(CHp(I,:,:).ge.thresh2,0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(0,y,0,:,:)=dim_num_n(CH(I,:,:).ge.thresh1,0)
  CH95p(1,y,0,:,:)=dim_num_n(CHp(I,:,:).ge.thresh2,0)
  delete(I)
end do

CH95pave=dim_avg_n_Wrap(CH95p,1)

  wks = gsn_open_wks("pdf","chaines_erai_ge"+thresh1+"_ge"+thresh2+"p_seasons_19792016_raster")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels             = levels
  res@gsnAddCyclic = False
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpFillOn              = False
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
  res@gsnRightString="" 
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False
  res@gsnLeftStringFontHeightF = 0.035
;  res@tmYRLabelsOn         = False              ; no right labels
;  res@tmYROn               = False              ; no right tickmarks
;  res@tmXBLabelsOn         = False              ; do not draw bottom labels
;  res@tmXBOn               = False              ; no bottom tickmarks

; now change the size of the tickmark labels
  res@tmXBLabelFontHeightF = 0.025                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.025    

; change label spacing to avoid overlap  
  res@tmXBTickSpacingF = 10.                      ; every 60 degrees
  res@tmYLTickSpacingF = 10.                      ; every 60 degrees


  plots1 = new(4,graphic)
  plots2 = new(4,graphic)

  res2=res
  res@cnLevels             = levels1
  res2@cnLevels             = levels2
;  res2@tmYLLabelsOn         = False              ; do not draw bottom labels
;  res2@tmYLOn               = False              ; no bottom tickmarks

  do p=0,3
   res@gsnLeftString=snames(p)
   res2@gsnLeftString=snames(p)
   plots1(p) = gsn_csm_contour_map(wks,CH95pave(0,p,:,:), res)     ; create plot
   plots2(p) = gsn_csm_contour_map(wks,CH95pave(1,p,:,:), res2)

;   if(p.eq.3) then
;     res@tmXBLabelsOn = True
;     res2@tmXBLabelsOn = True
;     res@tmXBOn = True
;     res2@tmXBOn = True
;   end if
  end do

  resP                  = True                   ; modify the panel plot
;  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPanelMainString = "Days with C-Haines >="+thresh1
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelFigureStrings= (/"a)","c)","b)","d)"/) ; add strings to panel
  resP@amJust   = "TopLeft"
  resP@gsnPanelLeft=0.2
  resP@gsnPanelRight=0.55
  resP@gsnFrame = False ; don't advance frame yet
  resP@gsnPanelMainPosXF  = 0.39
  resP@gsnPanelFigureStringsFontHeightF = 0.012
  resP@gsnPanelMainFontHeightF = 0.015
  resP@lbLabelFontHeightF = 0.012
 gsn_panel(wks,plots1,(/4,1/),resP)

  resP2=resP
  resP2@gsnPanelMainString = "Days with C-Haines >="+thresh2+"p"
  resP2@gsnPanelFigureStrings= (/"e)","f)","g)","h)"/) ; add strings to panel
  resP2@gsnPanelLeft=0.45
  resP2@gsnPanelRight=0.8
  resP2@gsnPanelMainPosXF  = 0.64
 gsn_panel(wks,plots2,(/4,1/),resP2)

frame(wks)


