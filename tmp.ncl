load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

levels1=(/0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95/)
levels2=ispan(0,100,10)

Fthresh=90
Cthresh=90

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon
CH=ifile->CHaines
CHp=ifile->CH_percentile

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI
FFDIp=ffile->FFDI_percentile

;FFDI2=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDI@_FillValue,FFDI)

;FFDIp=where(FFDIp.gt.90.and.FFDI.eq.0,FFDIp@_FillValue,FFDIp)

CHcor=new((/2,dimsizes(lat),dimsizes(lon)/),"float")
CHcor!1="lat"
CHcor&lat=lat
CHcor!2="lon"
CHcor&lon=lon

  CHcor(0,:,:)=spcorr_n(FFDI,CH,0)
  CHcor(1,:,:)=100*dim_num_n(CHp.ge.Cthresh.and.FFDIp.ge.Fthresh,0)/where(dim_num_n(FFDIp.ge.Fthresh,0).eq.0,CHcor@_FillValue,dim_num_n(FFDIp.ge.Fthresh,0))
 
  CHcor(1,:,:)=where(CHcor(0,:,:).eq.CHcor@_FillValue,CHcor@_FillValue,CHcor(1,:,:)) 

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 6
wks_type@wkPaperHeightF = 4
  wks = gsn_open_wks(wks_type,"Fig3_FFDICHcorr_propF"+Fthresh+"pC"+Cthresh+"p")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels             = levels
  res@gsnAddCyclic = False
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpFillOn              = True
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
  res@gsnLeftStringFontHeightF = 0.03
  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02
res@pmLabelBarOrthogonalPosF = .15
  res@cnLineLabelPlacementMode = "constant"
  res@cnLineLabelsOn = True
  res@cnLineDashSegLenF = 0.3



  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "PostDraw"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"


  plot = new(2,graphic)

   res@cnLevels=levels1
  nlevels = dimsizes(levels1)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(mod(i,2).eq.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

   res@gsnLeftString="a) Correlation coefficient (r)"
   plot(0) = gsn_csm_contour_map(wks,CHcor(0,:,:), res)     ; create plot

    delete([/res@cnLevels,flags,res@cnLevelFlags/])

  nlevels = dimsizes(levels2)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(mod(i,2).eq.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

   res@tmYLLabelFontColor = "white"   
   res@cnLevels=levels2
   res@gsnLeftString="b) Compound events (%)"
   plot(1) = gsn_csm_contour_map(wks,CHcor(1,:,:), res)     ; create plot

  resP                  = True                   ; modify the panel plot
;  resP@gsnMaximize=True
;  resP@amJust   = "TopLeft"
  resP@gsnPanelBottom = 0.05
  gsn_panel(wks,plot,(/1,2/),resP)
