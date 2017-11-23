load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

levels1=(/0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95/)
levels2=ispan(10,90,10)

levels1a=(/0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1/)
levels2a=ispan(0,100,10)

pthresh=95

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon

data=new((/2,dimsizes(time),dimsizes(lat),dimsizes(lon)/),"float")
data(0,:,:,:)=(/ifile->CHaines/)

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
data(1,:,:,:)=(/ffile->FFDI/)
data(1,:,:,:)=where(data(1,:,:,:).ge.999,data@_FillValue,data(1,:,:,:))

CHcor=new((/2,dimsizes(lat),dimsizes(lon)/),"float")
CHcor!1="lat"
CHcor&lat=lat
CHcor!2="lon"
CHcor&lon=lon

  CHcor(0,:,:)=spcorr_n(data(0,:,:,:),data(1,:,:,:),0)

CHpc=new((/2,dimsizes(lat),dimsizes(lon)/),"float")

   do x=0,1
    do jlat=0,dimsizes(lat)-1
    do ilon=0,dimsizes(lon)-1
      tmp1=ind(.not.ismissing(data(x,:,jlat,ilon)))
      if(dimsizes(tmp1).ge.0.95*dimsizes(time)) then
        CHsort=data(x,tmp1,jlat,ilon)
        qsort(CHsort)
        tmp=round(pthresh*dimsizes(CHsort)/100.,3)
        CHpc(x,jlat,ilon)=(/CHsort(tmp)/)
        delete([/tmp,CHsort/])
      end if
      delete(tmp1)
    end do
    end do
   end do

atmp=dim_num_n(data(0,:,:,:).gt.conform_dims(dimsizes(data(0,:,:,:)),CHpc(0,:,:),(/1,2/)).and.data(1,:,:,:).gt.conform_dims(dimsizes(data(1,:,:,:)),CHpc(1,:,:),(/1,2/)),0)
btmp=dim_num_n(data(1,:,:,:).gt.conform_dims(dimsizes(data(1,:,:,:)),CHpc(1,:,:),(/1,2/)),0)

CHcor(1,:,:)=(/100*atmp/where(btmp.eq.0,CHcor@_FillValue,btmp)/)
  CHcor(1,:,:)=where(CHcor(0,:,:).eq.CHcor@_FillValue,CHcor@_FillValue,CHcor(1,:,:)) 

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 9
wks_type@wkPaperHeightF = 6
  wks = gsn_open_wks(wks_type,"Fig2_FFDICHcorr_propF"+pthresh+"pC"+pthresh+"p_forsub")             ; send graphics to PNG file

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
;  res@gsnLeftStringFontHeightF = 0.03
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
 res@pmLabelBarWidthF = 0.7

; Label bar title

  res@lbTitleOn        =  True                ; turn on title
  res@lbTitlePosition  = "Bottom"              ; title position
;  res@lbTitleFontHeightF= .015                ; make title smaller
  res@lbTitleDirection = "Across"             ; title direction
;  res@lbLabelPosition  = "Center"             ; label position
  res@cnLabelBarEndStyle="IncludeMinMaxLabels"


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

;   res@gsnLeftString="a) Correlation coefficient (r)"
  res@lbTitleString="Correlation coefficient (r)"
  res@lbLabelStrings=levels1a
   plot(0) = gsn_csm_contour_map(wks,CHcor(0,:,:), res)     ; create plot

    delete([/res@cnLevels,flags,res@cnLevelFlags,res@lbLabelStrings/])

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
  res@lbLabelStrings=levels2a
;   res@gsnLeftString="b) Compound events (%)"
  res@lbTitleString="%"

   plot(1) = gsn_csm_contour_map(wks,CHcor(1,:,:), res)     ; create plot

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize=True
  resP@gsnPaperOrientation = "portrait"
  resP@amJust   = "TopLeft"
  resP@gsnPanelFigureStrings = (/"a)","b)"/)
  resP@gsnPanelFigureStringsPerimOn = False
;  resP@gsnPanelBottom = 0.05
  gsn_panel(wks,plot,(/1,2/),resP)
