load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

thresh=8
levels=ispan(0,50,5)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)

CHp=ifile->CHaines

years=ispan(1979,2016,1)

a=dimsizes(CHp)
CH95p=new((/dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!2="lat"
CH95p&lat=CHp&lat
CH95p!3="lon"
CH95p&lon=CHp&lon

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(y,m,:,:)=dim_num_n(CHp(I,:,:).ge.thresh,0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(y,0,:,:)=dim_num_n(CHp(I,:,:).ge.thresh,0)
  delete(I)
end do

CH95pave=dim_avg_n_Wrap(CH95p,0)

  wks = gsn_open_wks("pdf","chaines_erai_ge"+thresh+"_seasons_19792016_raster")             ; send graphics to PNG file

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
  res@cnFillMode = "RasterFill"
  res@cnLinesOn = False
 
  plot = new(4,graphic)

  do p=0,3
   res@tiMainString=snames(p)+" days with C-Haines >= "+thresh
   plot(p) = gsn_csm_contour_map(wks,CH95pave(p,:,:), res)     ; create plot
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  gsn_panel(wks,plot,(/4,1/),resP)
