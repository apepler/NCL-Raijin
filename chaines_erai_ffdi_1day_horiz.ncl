load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

date=20090207
tit="Black Saturday: 7 February 2009"
names=(/"C-Haines","C-Haines percentile","FFDI","FFDI percentile"/)
levels=(/(/3,6,7,8,9,10,11,12,13/),(/0,25,50,60,70,80,90,95,99/),(/10,20,30,40,50,60,70,80,90/),(/0,25,50,60,70,80,90,95,99/)/)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","w")
time=ifile->time
time2=cd_calendar(time,2)
T=ind(time2.eq.date)
lat=ifile->lat
lon=ifile->lon

; Get both FFDI and CHaines as value & percentile

CH=new((/4,dimsizes(lat),dimsizes(lon)/),"float")

CH(0,:,:)=ifile->CHaines(T,:,:)
CH(1,:,:)=ifile->CH_percentile(T,:,:)

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
CH(2,:,:)=ffile->FFDI(T,:,:)
CH(3,:,:)=ffile->FFDI_percentile(T,:,:)
CH(3,:,:)=where(CH(2,:,:).eq.0..and.CH(3,:,:).gt.90,CH@_FillValue,CH(3,:,:))

wks_type = "pdf"
;wks_type@wkOrientation = "landscape" 
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 4
  wks = gsn_open_wks(wks_type,"chaines_ffdi_"+date+"_horizontal")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
;  res@tiMainString         = "C-Haines index: "+date    ; main title
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels             = fspan(0.,18.,10)
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
  res@gsnLeftString=""
  res@gsnRightString=""
  res@gsnCenterStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.025
  res@pmLabelBarOrthogonalPosF = .15
  res@gsnPaperOrientation = "landscape"

  plot = new(4,graphic)

  do p=0,3
   res@cnLevels = levels(p,:) 
   res@gsnCenterString=names(p)
   plot(p) = gsn_csm_contour_map(wks,CH(p,:,:), res)     ; create plot
   delete(res@cnLevels)
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
;  resP@txString           = tit
  resP@gsnPanelMainString = tit
  resP@gsnPaperOrientation = "portrait"
  resP@gsnPanelXWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/1,4/),resP)

