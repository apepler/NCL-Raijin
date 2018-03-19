load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)
ofile="MUCAPE_ERAI_Australia_methodpanel"

years=ispan(2000,2016,1)
lonlim=(/110,160/)
latlim=(/-45,-10/)
levels=(/10,25,50,75,100,150,200,250,300,400,500,600,700,800,900,1000/)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

ifile=addfile("mucape_20002016_monthly.nc","r")
mucape=ifile->mucape

;mucape=new((/dimsizes(types),dimsizes(years),12,dimsizes(lat),dimsizes(lon)/),"float")

;do t=0,dimsizes(types)-1
; do y=0,dimsizes(years)-1
;  do m=0,11
;   ifile=addfile(cdir+types(t)+"/mucape_"+years(y)+sprinti("%0.2i",m+1)+".nc","r")
;   mucape(t,y,m,:,:)=(/dim_avg_n_Wrap(ifile->mucape,0)/) 
;  end do
; end do
;end do

Cmean=dim_avg_n_Wrap(mucape,(/1,2/))
Cmean!1="lat"
Cmean&lat=lat
Cmean!2="lon"
Cmean&lon=lon

; Plotting

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
  res@gsnLeftString="ERAI"
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 6
  wks = gsn_open_wks(wks_type,ofile)             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")

  plot = new(dimsizes(types),graphic)

  do p=0,dimsizes(types)-1
   res@gsnRightString=tnames(p)
   plot(p) = gsn_csm_contour_map(wks,Cmean(p,:,:), res)     ; create plot
  end do

  resP                  = True                   ; modify the panel plot
resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="Average Most Unstable CAPE (J/Kg)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelRowSpec = True
  gsn_panel(wks,plot,(/2,3/),resP)

;mucape!0="type"
;mucape!1="year"
;mucape&year=years
;mucape!2="month"
;mucape&month=ispan(1,12,1)
;mucape!3="lon"
;mucape&lon=lon
;mucape!4="lat"
;mucape&lat=lat

mucape@standard_name = "Most Unstable CAPE"
mucape@units = "J/kg"

;a=addfile("mucape_20002016_monthly.nc","c")
;a->mucape=mucape
;a->method=tnames


