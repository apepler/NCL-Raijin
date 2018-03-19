load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

ofile="ERAIcorr_MUCAP_CHaines_Australia_methodpanel"

years=ispan(2000,2016,1)
lonlim=(/110,160/)
latlim=(/-45,-10/)
levels=(/0,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)

; Setup - need to have the latitude and longitude
; Get CHaines to use for correlations

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_components.nc","r")
lat=ifile->lat
lon=ifile->lon

I=ind(lon.ge.min(lonlim).and.lon.le.max(lonlim))
J=ind(lat.ge.min(latlim).and.lat.le.max(latlim))

time=ifile->time
time2=cd_calendar(time,0)
Y=ind(time2(:,0).ge.min(years).and.time2(:,0).le.max(years))
CHaines=ifile->CHaines(Y,J,I)

delete([/ifile,time,time2/])

; Set up plot first

  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@gsnAddCyclic = False
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpFillOn              = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res@mpLimitMode        = "Corners"               ; Portion of map to zoom
  res@mpLeftCornerLatF   = latlim(0)
  res@mpLeftCornerLonF   = lonlim(0)
  res@mpRightCornerLatF  = latlim(1)
  res@mpRightCornerLonF  = lonlim(1)
  res@mpGridLineDashPattern = 2
  res@mpGridLineThicknessF = 0.5
  res@mpProjection = "CylindricalEquidistant"
  res@cnLevels             = levels
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

  wks_type = "pdf"
  ;wks_type@wkOrientation = "landscape"
  wks_type@wkPaperWidthF = 10
  wks_type@wkPaperHeightF = 8
  wks = gsn_open_wks(wks_type,ofile)             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")

  plot=new(5,graphic)

  corr=new((/dimsizes(types),dimsizes(J),dimsizes(I)/),"float")
  corr!1="lat"
  corr&lat=lat(J)
  corr!2="lon"
  corr&lon=lon(I)

  do p=0,dimsizes(types)-1
   print(types(p))
   ifile=addfiles( systemfunc (" ls -1 " + cdir+ types(p) + "/mucape_2?????.nc"), "r")
    mucape=ifile[:]->mucape(:,J,I)

    corr(p,:,:)=(/spcorr_n(CHaines,mucape,0)/)

   res@gsnRightString=tnames(p)
   plot(p)=gsn_csm_contour_map_ce(wks,corr(p,:,:), res)     ; create plot

   delete(mucape)
 end do

  resP                  = True                   ; modify the panel plot
resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="Corr between 6-hrly CHaines and Most Unstable CAPE"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelRowSpec = True
  gsn_panel(wks,plot,(/2,3/),resP)

