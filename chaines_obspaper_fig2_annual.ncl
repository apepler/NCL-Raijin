load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

pthresh=95
levels=ispan(0,20,2)

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

years=ispan(1979,2016,1)

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


CH95p=new((/3,dimsizes(years),dimsizes(lat),dimsizes(lon)/),"integer")
CH95p!2="lat"
CH95p&lat=lat
CH95p!3="lon"
CH95p&lon=lon

seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
  I=ind(time2(:,0).eq.years(y))
  do x=0,1
    CH95p(x,y,:,:)=dim_num_n(data(x,I,:,:).ge.conform_dims(dimsizes(data(x,I,:,:)),CHpc(x,:,:),(/1,2/)),0)
  end do
  CH95p(2,y,:,:)=dim_num_n(data(0,I,:,:).ge.conform_dims(dimsizes(data(0,I,:,:)),CHpc(0,:,:),(/1,2/)).and.data(1,I,:,:).ge.conform_dims(dimsizes(data(1,I,:,:)),CHpc(1,:,:),(/1,2/)),0)
  delete(I)
end do

  wks = gsn_open_wks("pdf","Fig2_CHainesge"+pthresh+"p_FFDIge"+pthresh+"p_annual_recalc")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels
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
;  res@cnFillMode = "RasterFill"
;  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = True
;  res@gsnMaximize=True
  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "Draw"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"
  res@gsnStringFontHeightF = 0.035
  res@gsnRightString="Days"
  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02    

  plots = new(3,graphic)

  snames=(/"FFDI >= "+pthresh+"p","CH >="+pthresh+"p","Both"/)
  do x=0,2

;  Make the array as we go
  tmp=dim_avg_n_Wrap(CH95p(x,:,:,:),0)
  tmp!0="lat"
  tmp&lat=lat
  tmp!1="lon"
  tmp&lon=lon

  res@gsnLeftString=snames(x)
  plots(x) = gsn_csm_contour_map(wks,tmp, res)     ; create plot

  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnPanelLabelBar    = True
;  resP@gsnPanelYWhiteSpacePercent = 3
  resP@amJust   = "TopLeft"
  gsn_panel(wks,plots,(/1,3/),resP)

frame(wks)


