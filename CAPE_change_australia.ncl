load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_adiabatic_nomix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Adiabatic, No Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

years=ispan(1979,2016,1)
levels=ispan(-60,60,10)
I1=ind(years.le.1997)
I2=ind(years.ge.1998)

sname="djf"
seasons=(/12,2/)
ofile="MUCAPE_ERAI_Australia_methodpanel_change_1997_9818_"+sname


; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

lonlim=(/110,160/)
latlim=(/-45,-10/)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

ifile=addfile("mucape_19792016_monthly.nc","r")
mucape=ifile->mucape

Cmean=new((/dimsizes(types),dimsizes(years),dimsizes(lat),dimsizes(lon)/),"float")
Cmean!1="year"
Cmean!2="lat"
Cmean&lat=lat
Cmean!3="lon"
Cmean&lon=lon

if(seasons(1).ge.seasons(0)) then
 mlist=ispan(seasons(0),seasons(1),1)-1
 Cmean=(/dim_avg_n_Wrap(mucape(:,:,mlist,:,:),2)/)
else
 do y=1,dimsizes(years)-1
  list1=ispan(seasons(0)-1,11,1)
  list2=ispan(0,seasons(1)-1,1)
  if(dimsizes(list1).eq.1) then
    Cmean(:,y,:,:)=(mucape(:,y-1,list1,:,:)+dim_sum_n_Wrap(mucape(:,y,list2,:,:),1))/(dimsizes(list1)+dimsizes(list2))
  else
    Cmean(:,y,:,:)=(dim_sum_n_Wrap(mucape(:,y-1,list1,:,:),1)+dim_sum_n_Wrap(mucape(:,y,list2,:,:),1))/(dimsizes(list1)+dimsizes(list2))
 end if

 end do
end if 

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

opts_psl = True
     opts_psl@cnLineColor       = "Black"
     opts_psl@cnInfoLabelOn     = False
     opts_psl@cnLineLabelFontHeightF = 0.01
     opts_psl@cnLineLabelPerimOn = False
     opts_psl@gsnContourLineThicknessesScale = 4
     opts_psl@gsnDraw         = False
     opts_psl@gsnFrame        = False
     opts_psl@cnLevelSelectionMode = "ExplicitLevels"
     opts_psl@cnLevels             = (/ -100.,0.05,100./)
     opts_psl@cnLineDrawOrder="Draw"
;     opts_psl@gsnContourPosLineDashPattern = 2


wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 6
  wks = gsn_open_wks(wks_type,ofile)             ; send graphics to PNG file
  gsn_define_colormap(wks,"MPL_RdBu")
;gsn_reverse_colormap(wks)

  plot = new(dimsizes(types),graphic)

  do p=0,dimsizes(types)-1
   res@gsnRightString=tnames(p)

   tmp1=dim_avg_n_Wrap(Cmean(p,I1,:,:),0)
   tmp2=dim_avg_n_Wrap(Cmean(p,I2,:,:),0)
   tmp=100*((tmp2/where(tmp1.eq.0,Cmean@_FillValue,tmp1))-1)
   tmp!0="lat"
   tmp&lat=lat
   tmp!1="lon"
   tmp&lon=lon

  plot(p) = gsn_csm_contour_map(wks,tmp, res)     ; create plot

; And also the significance
  r1a=Cmean(p,I1,:,:)
  r2a=Cmean(p,I2,:,:)
  xAve = dim_avg_n(r1a,0)              ; calculate means at each grid point
  yAve = dim_avg_n(r2a,0)
  xVar = dim_variance_n(r1a,0)         ; calculate variances
  yVar = dim_variance_n(r2a,0)
  sigr = 0.05                        ; critical sig lvl for r
  xEqv = equiv_sample_size (r1a(lat|:,lon|:,year|:), sigr,0)
  yEqv = equiv_sample_size (r2a(lat|:,lon|:,year|:), sigr,0)
;  xN   = wgt_areaave (xEqv, 1., 1., 0)    ; wgty could be gaussian weights
;  yN   = wgt_areaave (yEqv, 1., 1., 0)
  iflag= False                        ; population variance similar
  pval = ttest(xAve,xVar,xEqv, yAve,yVar,yEqv, iflag, False)
  pval=where(dim_avg_n_Wrap(Cmean(p,:,:,:),0).lt.0.5,pval@_FillValue,pval)
  pval!0="lat"
  pval!1="lon"
  pval&lat = lat
  pval&lon = lon

  contour_psl = gsn_csm_contour(wks,pval,opts_psl)
  overlay(plot(p), contour_psl)

  end do

  resP                  = True                   ; modify the panel plot
resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="Long-term change in Most Unstable CAPE (%)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelRowSpec = True
  gsn_panel(wks,plot,(/3,3/),resP)



