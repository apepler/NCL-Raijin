load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

odir="/g/data/eg3/asp561/ERAI/"

years=ispan(2000,2016,1)
lonlim=(/110,160/)
latlim=(/-45,-10/)
levels=(/10,25,50,75,100,150,200,250,300,400,500,600,700,800,900,1000/)
hours=(/0,6,12,18/)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

a=addfile(odir+"mucape_20002016_hourly.nc","r")
mucape=a->mucape

snames=(/"DJF","MAM","JJA","SON"/)
smonths=new((/4,3/),"integer")
smonths(0,:)=(/12,1,2/)-1
smonths(1,:)=ispan(3,5,1)-1
smonths(2,:)=ispan(6,8,1)-1
smonths(3,:)=ispan(9,11,1)-1

Cmean=new((/dimsizes(types),dimsizes(snames),dimsizes(hours),dimsizes(lat),dimsizes(lon)/),"float")

do s=0,3
  Cmean(:,s,:,:,:)=dim_avg_n_Wrap(mucape(:,:,smonths(s,:),:,:,:),(/1,2/))
end do

Cmean!0="variable"
Cmean!1="season"
Cmean!2="hour"
Cmean!3="lat"
Cmean&lat=lat
Cmean!4="lon"
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
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False
  res@gsnRightString=""

  wks_type = "pdf"
  ;wks_type@wkOrientation = "landscape"
  wks_type@wkPaperWidthF = 10
  wks_type@wkPaperHeightF = 8

do t=0,dimsizes(types)
  ofile="MUCAPE_ERAI_Australia_"+types(t)+"_seasonhour"
  wks = gsn_open_wks(wks_type,ofile)             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")

  lablets=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)
  plot = new(16,graphic)

  do m=0,3
   do h=0,3
    p=m*4+h

    res@gsnLeftString=lablets(h*4+m)+") "+snames(m)+" "+hours(h)+" UTC CAPE"
;    res@gsnRightString=tnames(t)

    plot(p) = gsn_csm_contour_map(wks,Cmean(t,m,h,:,:), res)     ; create plot
   end do
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
  gsn_panel(wks,plot,(/4,4/),resP)
  delete(wks)
end do

