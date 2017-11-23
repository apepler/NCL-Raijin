load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

date=20030118
years=ispan(1979,2016,1)
pthresh=95
levels=ispan(1,13,1)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,2)
T=ind(time2.eq.date)
lat=ifile->lat
lon=ifile->lon
CH=ifile->CHaines

mfile=addfile("/g/data1/ub4/erai/netcdf/fx/ei_invariant_075x075_90N0E90S35925E.nc","r")
lat2=mfile->latitude
lon2=mfile->longitude
I=ind(lon2.ge.min(lon).and.lon2.le.max(lon))
J=ind(lat2.ge.min(lat).and.lat2.le.max(lat))
landmask=short2flt(mfile->lsm(0,J,I))

wks_type = "eps"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 5
wks_type@wkPaperHeightF = 4.6

  wks = gsn_open_wks(wks_type,"Fig1a_chaines_erai_p"+pthresh+"_forsub")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@gsnMaximize=True
  res@gsnPaperMargin = 0.2
  res@gsnPaperOrientation = "portrait" 
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels
  res@gsnAddCyclic = False
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpFillOn              = True
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "Draw"           ; Draw map outlines last
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
  res@gsnRightString=""
;  res@lbLabelBarOn        = False
  res@cnLineLabelPlacementMode = "constant"
  res@cnLineLabelsOn = True
  res@cnLineDashSegLenF = 0.3
  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@lbTitleOn        =  True                ; turn on title
  res@lbTitlePosition  = "Bottom"              ; title position
;  res@lbTitleFontHeightF= .015                ; make title smaller
  res@lbTitleDirection = "Across"             ; title direction
;  res@lbLabelPosition  = "Center"             ; label position
  res@cnLabelBarEndStyle="IncludeMinMaxLabels"
  res@lbLabelStrings=ispan(0,14,1)
  res@tmXBLabelFontHeightF = 0.015                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.015
  res@lbTitleFontHeightF = 0.015
res@pmLabelBarOrthogonalPosF = .15
res@pmLabelBarWidthF = 0.55
res@pmLabelBarHeightF = 0.05

  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "Draw"

;  res@mpAreaMaskingOn       = True
;  res@mpGridMaskMode        = "MaskLand"
;  res@mpMaskAreaSpecifiers = "water"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"
;
;  res@mpFillAreaSpecifiers        = (/"water",       "land" /)
;  res@mpSpecifiedFillColors       = (/"white","white"/)


                                               ; last.
  nlevels = dimsizes(levels)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(levels(i).lt.6.or.mod(levels(i),2).ne.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

  counts=round(pthresh*dimsizes(CH&time)/100.,3)

;   res@gsnLeftString="a) "+pthresh+"th percentile of CH based on ERAI reanalysis"
    res@lbTitleString=pthresh+"~S~th~N~ percentile of CH based on reanalysis"
    a=dimsizes(CH)-1
    CHpc=new((/dimsizes(lat),dimsizes(lon)/),typeof(CH))
    CHpc!0="lat"
    CHpc&lat=lat
    CHpc!1="lon"
    CHpc&lon=lon

    do jlat=0,dimsizes(lat)-1
    do ilon=0,dimsizes(lon)-1
      CHsort=CH(ind(.not.ismissing(CH(:,jlat,ilon))),jlat,ilon)
      if(dimsizes(CHsort).ge.0.95*dimsizes(CH&time)) then
        qsort(CHsort)
        tmp=round(counts*1.0*dimsizes(CHsort)/dimsizes(CH&time),3)
        CHpc(jlat,ilon)=(/CHsort(tmp)/)
      end if
      delete(CHsort)
    end do
    end do

    g6=where(CHpc.ge.6,1,0)
    print(sum(g6*landmask)/sum(landmask))
 
;   CHpc=mask(CHpc,landmask,0)

   plot=gsn_csm_contour_map(wks,CHpc, res)     ; create plot

   ptnames=(/"Alice Springs","Townsville","Perth","Wagga Wagga"/)
   ptlats=(/-23.7951,-19.2483,-31.9275,-35.1583/)
   ptlons=(/133.8890,146.7661,115.9764,147.4575/)

   gsres = True
   gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
   gsres@gsMarkerSizeF = 0.015

   txres               = True
   txres@txFontHeightF = 0.025
   txres@txJust="CenterRight"
   txres@txFontThicknessF = 3.   

  dum1=gsn_add_polymarker(wks,plot,ptlons,ptlats,gsres)
  dum2=gsn_add_text(wks,plot,ptnames,ptlons-1,ptlats,txres)

   txres2               = True
   txres2@txFontHeightF = 0.04
  dum3=gsn_add_text(wks,plot,"a)",113,-13,txres2)

   draw(plot)
   frame(wks)
