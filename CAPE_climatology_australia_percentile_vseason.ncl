load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

pthresh=95
ofile="MUCAPE_ERAI_"+pthresh+"pc_Australia_methodpanel_vseason"

years=ispan(2000,2016,1)
lonlim=(/110,160/)
latlim=(/-45,-10/)

levels=ispan(1,15,1)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

a=addfile("mucape_"+pthresh+"pc_20002016_monthly.nc","r")
mucape=(a->mucape)/4.

;mucape=new((/dimsizes(types),dimsizes(years),12,dimsizes(lat),dimsizes(lon)/),"integer")
;
;do t=0,dimsizes(types)-1
; a=addfile("mucape_"+types(t)+"_95pc.nc","r")
; mucapepc=a->mucape_95pc
;
; do y=0,dimsizes(years)-1
;  do m=0,11
;   ifile=addfile(cdir+types(t)+"/mucape_"+years(y)+sprinti("%0.2i",m+1)+".nc","r")
;   dailymc=ifile->mucape
;   mucape(t,y,m,:,:)=(/dim_num_n(dailymc.gt.conform_dims(dimsizes(dailymc),mucapepc,(/1,2/)),0)/)
;   delete(dailymc)
;  end do
; end do
;end do

snames=(/"DJF","MAM","JJA","SON"/)
smonths=new((/4,3/),"integer")
smonths(0,:)=(/12,1,2/)-1
smonths(1,:)=ispan(3,5,1)-1
smonths(2,:)=ispan(6,8,1)-1
smonths(3,:)=ispan(9,11,1)-1

seasons1=(/12,3,6,9/)

Cmean=new((/dimsizes(types),dimsizes(snames),dimsizes(lat),dimsizes(lon)/),"float")

do s=1,3
  Cmean(:,s,:,:)=dim_avg_n_Wrap(dim_sum_n_Wrap(mucape(:,:,smonths(s,:),:,:),2),1)
end do

y=dimsizes(years)-1
Cmean(:,0,:,:)=dim_avg_n_Wrap(dim_sum_n_Wrap(mucape(:,1:y,0:1,:,:),2)+mucape(:,0:(y-1),11,:,:),1)

Cmean!0="variable"
Cmean!1="season"
Cmean!2="lat"
Cmean&lat=lat
Cmean!3="lon"
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

  wks_type = "pdf"
  ;wks_type@wkOrientation = "landscape"
  wks_type@wkPaperWidthF = 10
  wks_type@wkPaperHeightF = 8
  wks = gsn_open_wks(wks_type,ofile)             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")

  lablets=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)
  plot = new(20,graphic)

  do m=0,3
   do x=0,4
    p=m*5+x

    res@gsnLeftString=lablets(x*4+m)+") "+snames(m)+" CAPE"
    res@gsnRightString=tnames(x)

    plot(p) = gsn_csm_contour_map(wks,Cmean(x,m,:,:), res)     ; create plot
   end do
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="Average number of days with MUCAPE >= "+pthresh+"th percentile"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
  resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/4,5/),resP)

;mucape!0="type"
;mucape!1="year"
;mucape&year=years
;mucape!2="month"
;mucape&month=ispan(1,12,1)
;mucape!3="lat"
;mucape&lat=lat
;mucape!4="lon"
;mucape&lon=lon

;mucape@standard_name = "MUCAPE days above "+pthresh+"th percentile"
;mucape@units = "J/kg"

;a=addfile("mucape_"+pthresh+"pc_20002016_monthly.nc","c")
;a->mucape=mucape

