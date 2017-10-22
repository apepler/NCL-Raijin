load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

levels1=(/0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95/)
levels2=ispan(0,100,10)

Fthresh=90
Cthresh=8

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon
CH=ifile->CHaines

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI
FFDIp=ffile->FFDI_percentile

FFDI2=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDI@_FillValue,FFDI)

FFDIp=where(FFDIp.gt.90.and.FFDI.eq.0,FFDIp@_FillValue,FFDIp)

CHcor=new((/2,dimsizes(lat),dimsizes(lon)/),"float")
CHcor!1="lat"
CHcor&lat=lat
CHcor!2="lon"
CHcor&lon=lon

  CHcor(0,:,:)=spcorr_n(FFDI2,CH,0)
  CHcor(1,:,:)=100*dim_num_n(CH.ge.Cthresh.and.FFDIp.ge.Fthresh,0)/where(dim_num_n(FFDIp.ge.Fthresh,0).eq.0,CHcor@_FillValue,dim_num_n(FFDIp.ge.Fthresh,0))
 
  CHcor(1,:,:)=where(CHcor(0,:,:).eq.CHcor@_FillValue,CHcor@_FillValue,CHcor(1,:,:)) 

  wks = gsn_open_wks("pdf","FFDI_CH_corr_propF"+Fthresh+"pC"+Cthresh+"_19792016")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels             = levels
  res@gsnAddCyclic = False
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpFillOn              = False
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
  res@gsnLeftStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.025
res@pmLabelBarOrthogonalPosF = .15

  plot = new(2,graphic)

   res@cnLevels=levels1
   res@tiMainString="Correlation between daily FFDI and C-Haines"
   plot(0) = gsn_csm_contour_map(wks,CHcor(0,:,:), res)     ; create plot
    delete(res@cnLevels)
   res@cnLevels=levels2
   res@tiMainString="% of FFDI "+Fthresh+"p days where C-Haines >="+Cthresh
   plot(1) = gsn_csm_contour_map(wks,CHcor(1,:,:), res)     ; create plot

  resP                  = True                   ; modify the panel plot
  resP@gsnPanelFigureStrings= (/"a)","c)","b)","d)"/) ; add strings to panel
  resP@amJust   = "TopLeft"
  resP@gsnPanelBottom = 0.05
  gsn_panel(wks,plot,(/1,2/),resP)
