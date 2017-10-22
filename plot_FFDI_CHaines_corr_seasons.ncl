load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

levels=(/0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95/)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon
CH=ifile->CHaines

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI
FFDI2=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDI@_FillValue,FFDI)

CHcor=new((/4,dimsizes(lat),dimsizes(lon)/),"float")
CHcor!1="lat"
CHcor&lat=lat
CHcor!2="lon"
CHcor&lon=lon
skill=CHcor

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do m=1,3
  I=ind(time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CHcor(m,:,:)=spcorr_n(FFDI2(I,:,:),CH(I,:,:),0)
  delete(I)
end do

  I=ind(time2(:,1).eq.12.or.time2(:,1).le.2)
  CHcor(0,:,:)=spcorr_n(FFDI2(I,:,:),CH(I,:,:),0)
  delete(I)

  wks = gsn_open_wks("pdf","FFDI_CH_corr_seasons_19792016")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels
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

  plot = new(4,graphic)

  do p=0,3
   res@tiMainString=snames(p)+" correlation between daily FFDI and C-Haines"
   plot(p) = gsn_csm_contour_map(wks,CHcor(p,:,:), res)     ; create plot

  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  gsn_panel(wks,plot,(/4,1/),resP)
