load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
dir2="/g/data/eg3/asp561/ERAI/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

pthresh=95
ofile="MUCAPE_ERAI_Australia_"+pthresh+"pc_methodcomp_CSI"

years=ispan(2000,2016,1)
lonlim=(/110,160/)
latlim=(/-45,-10/)
levels=fspan(0.1,0.9,9)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

I=ind(lon.ge.min(lonlim).and.lon.le.max(lonlim))
J=ind(lat.ge.min(latlim).and.lat.le.max(latlim))

delete(ifile)

ifile=addfiles( systemfunc (" ls -1 " + cdir+ types(4) + "/mucape_2?????.nc"), "r")
time=cd_calendar(ifile[:]->time,0)
tlen=dimsizes(time)
delete(ifile)

data=new((/dimsizes(types),tlen(0),dimsizes(J),dimsizes(I)/),"float")
data!0="Type"
data!1="time"
data!2="lat"
data&lat=lat(J)
data!3="lon"
data&lon=lon(I)

do t=0,dimsizes(types)-1
 print(types(t))
 ifile=addfiles( systemfunc (" ls -1 " + cdir+ types(t) + "/mucape_2?????.nc"), "r")
 mucape=ifile[:]->mucape(:,J,I) 
; print(dimsizes(mucape))
; print(max(mucape))
 a=addfile(dir2+"mucape_"+types(t)+"_95pc.nc","r")
 mucapepc=a->mucape_95pc(J,I)
; print(dimsizes(mucapepc))
; print(max(mucapepc))
 data(t,:,:,:)=where(mucape.gt.conform_dims(dimsizes(mucape),mucapepc,(/1,2/)),1,0)
 print(avg(data(t,:,:,:)))
end do

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
  res@cnLinesOn = True
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = True

  wks_type = "pdf"
  ;wks_type@wkOrientation = "landscape"
  wks_type@wkPaperWidthF = 10
  wks_type@wkPaperHeightF = 10
  wks = gsn_open_wks(wks_type,ofile)             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")

  plot=new(25,graphic)

  do x=0,dimsizes(types)-1
   do y=0,dimsizes(types)-1
    p=x*5+y
    res@gsnLeftString=tnames(x)
    res@gsnRightString=tnames(y)
    tmp1=1.0*dim_sum_n_Wrap(where(data(x,:,:,:).eq.1.and.data(y,:,:,:).eq.1,1,0),0)
    tmp2=1.0*dim_sum_n_Wrap(where(data(x,:,:,:).eq.1.or.data(y,:,:,:).eq.1,1,0),0)
    print(dimsizes(tmp1))
    print(dimsizes(tmp2))
    tmp2@_FillValue=-999999
    CSI=tmp1/where(tmp2.eq.0,tmp2@_FillValue,tmp2)
    print(p+": "+avg(tmp1)+" / "+avg(tmp2)+" = "+avg(CSI))
    CSI!0="lat"
    CSI&lat=lat(J)
    CSI!1="lon"
    CSI&lon=lon(I)

    plot(p) = gsn_csm_contour_map(wks,CSI, res)
   end do
 end do

  resP                  = True                   ; modify the panel plot
resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="CSI for hours with MUCAPE > 95th percentile"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/5,5/),resP)

