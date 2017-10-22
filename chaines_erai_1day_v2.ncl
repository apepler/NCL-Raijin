load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

date=20030118

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","w")
time=ifile->time
time2=cd_calendar(time,2)
T=ind(time2.eq.date)

CH=ifile->CHaines(T,:,:)
CHp=ifile->CH_percentile(T,:,:)

  wks = gsn_open_wks("pdf","chaines_"+date)             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@tiMainString         = "C-Haines index: "+date    ; main title
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = fspan(0.,18.,10)
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

  plot = new(2,graphic)
  plot(0) = gsn_csm_contour_map(wks,CH, res)     ; create plot

  delete(res@cnLevels)
  lev=fspan(10.,110.,11)
  lev(9)=95.
  lev(10)=99.
  res@cnLevels             = lev
  res@tiMainString         = "C-Haines index percentile: "+date    ; main title

  plot(1) = gsn_csm_contour_map(wks,CHp, res)     ; create plot
  gsn_panel(wks,plot,(/1,2/),False)
