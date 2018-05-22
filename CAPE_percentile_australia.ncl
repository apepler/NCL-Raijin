load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_adiabatic_nomix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Adiabatic, No Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

pthresh=95
ofile="MUCAPE_ERAI_Australia_"+pthresh+"pc_methodpanel_19792016"

years=ispan(1979,2016,1)
lonlim=(/110,160/)
latlim=(/-45,-10/)
levels=(/50,100,200,500,1000,1500,2000,2500,3000/)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

delete(ifile)

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

  plot=new(dimsizes(types),graphic)

  do p=0,dimsizes(types)-1
   print(p)
   if(p.eq.5) then
    ifile=addfiles( systemfunc (" ls -1 " + cdir+ types(p) + "/mucape_*.nc"), "r")
    mucape=ifile[:]->mucape
    counts=round(pthresh*dimsizes(mucape&time)/100.,3)
    mucapepc=new((/dimsizes(lat),dimsizes(lon)/),typeof(mucape))
    mucapepc!0="lat"
    mucapepc&lat=lat
    mucapepc!1="lon"
    mucapepc&lon=lon

    do jlat=0,dimsizes(lat)-1
     do ilon=0,dimsizes(lon)-1
      mucapesort=mucape(ind(.not.ismissing(mucape(:,jlat,ilon))),jlat,ilon)
      if(dimsizes(mucapesort).ge.0.95*dimsizes(mucape&time)) then
        qsort(mucapesort)
        tmp=round(counts*1.0*dimsizes(mucapesort)/dimsizes(mucape&time),3)
        mucapepc(jlat,ilon)=(/mucapesort(tmp)/)
      end if
      delete(mucapesort)
    end do
    end do

    mucapepc@standard_name = "95th percentile of Most Unstable CAPE"
    mucapepc@units = "J/kg"
 
    a=addfile("mucape_"+types(p)+"_95pc_19792016.nc","c")
    a->mucape_95pc=mucapepc
   else 
    a=addfile("mucape_"+types(p)+"_95pc_19792016.nc","r")
    mucapepc=a->mucape_95pc
   end if

   res@gsnLeftString=""
   res@gsnRightString=tnames(p)
   plot(p)=gsn_csm_contour_map_ce(wks,mucapepc, res)     ; create plot

 end do

  resP                  = True                   ; modify the panel plot
resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString=pthresh+"th pctile of Most Unstable CAPE (J/Kg)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelRowSpec = True
  gsn_panel(wks,plot,(/3,3/),resP)

