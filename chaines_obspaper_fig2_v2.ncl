load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Cthresh=8
pthresh=95
levels1=ispan(0,80,10)
levels2=ispan(0,20,3)

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

CH95p=new((/3,dimsizes(years),4,dimsizes(lat),dimsizes(lon)/),"integer")
CH95p!3="lat"
CH95p&lat=lat
CH95p!4="lon"
CH95p&lon=lon

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(0,y,m,:,:)=dim_num_n(data(0,I,:,:).ge.Cthresh,0)
  CH95p(1,y,m,:,:)=dim_num_n(data(0,I,:,:).ge.conform_dims(dimsizes(data(0,I,:,:)),CHpc(0,:,:),(/1,2/)),0)
  CH95p(2,y,m,:,:)=dim_num_n(data(0,I,:,:).ge.conform_dims(dimsizes(data(0,I,:,:)),CHpc(0,:,:),(/1,2/)).and.data(1,I,:,:).ge.conform_dims(dimsizes(data(1,I,:,:)),CHpc(1,:,:),(/1,2/)),0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(0,y,0,:,:)=dim_num_n(data(0,I,:,:).ge.Cthresh,0)
  CH95p(1,y,0,:,:)=dim_num_n(data(0,I,:,:).ge.conform_dims(dimsizes(data(0,I,:,:)),CHpc(0,:,:),(/1,2/)),0)
  CH95p(2,y,0,:,:)=dim_num_n(data(0,I,:,:).ge.conform_dims(dimsizes(data(0,I,:,:)),CHpc(0,:,:),(/1,2/)).and.data(1,I,:,:).ge.conform_dims(dimsizes(data(1,I,:,:)),CHpc(1,:,:),(/1,2/)),0)

  delete(I)
end do

CH95pave=dim_avg_n_Wrap(CH95p,0)

  wks = gsn_open_wks("pdf","Fig2_CHainesge"+Cthresh+"_"+pthresh+"p_FFDIge"+pthresh+"p_recalc")             ; send graphics to PNG file

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
;  res@cnFillMode = "RasterFill"
;  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False
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

  plots1 = new(4,graphic)
  plots2 = new(4,graphic)
  plots3 = new(4,graphic)

res2=res
res@cnLevels             = levels1
res2@cnLevels             = levels2
res2@tmYLLabelFontColor = "white"

labs1=(/"a","b","c","d"/)
labs2=(/"e","f","g","h"/)
labs3=(/"i","j","k","l"/)

  do m=0,3

;  Make the array as we go
  tmp=dim_avg_n_Wrap(CH95p(0,:,m,:,:),0)
  tmp!0="lat"
  tmp&lat=lat
  tmp!1="lon"
  tmp&lon=lon

  if(m.lt.3) then
   res@tmXBLabelFontColor = "white"
   res2@tmXBLabelFontColor = "white"
  else
   res@tmXBLabelFontColor = "black"
   res2@tmXBLabelFontColor = "black"
  end if

  res@gsnLeftString=labs1(m)+") "+snames(m)
  plots1(m) = gsn_csm_contour_map(wks,tmp, res)     ; create plot

  tmp=(/dim_avg_n_Wrap(CH95p(1,:,m,:,:),0)/)
  res2@gsnLeftString=labs2(m)+") "+snames(m)
  plots2(m) = gsn_csm_contour_map(wks,tmp, res2)     ; create plot

  tmp=(/dim_avg_n_Wrap(CH95p(2,:,m,:,:),0)/)
  res2@gsnLeftString=labs3(m)+") "+snames(m)
  plots3(m) = gsn_csm_contour_map(wks,tmp, res2)     ; create plot

  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnPanelMainString = "CH >= "+Cthresh
  resP@gsnPanelLabelBar    = True
;  resP@gsnPanelYWhiteSpacePercent = 3
  resP@amJust   = "TopLeft"
  resP@gsnPanelLeft=0.15
  resP@gsnPanelRight=0.4
  resP@gsnFrame = False ; don't advance frame yet
  resP@gsnPanelMainPosXF  = 0.28
  resP@gsnPanelFigureStringsFontHeightF = 0.012
  resP@gsnPanelMainFontHeightF = 0.018
  resP@pmLabelBarWidthF = 0.2
  resP@lbLabelFontHeightF = 0.012
  gsn_panel(wks,plots1,(/4,1/),resP)

  resP2=resP
  resP2@gsnPanelMainString = "CH >="+pthresh+"p"
  resP2@gsnPanelLeft=0.4
  resP2@gsnPanelRight=0.65
  resP2@gsnPanelMainPosXF  = 0.52
 gsn_panel(wks,plots2,(/4,1/),resP2)

  resP2=resP
  resP2@gsnPanelMainString = "CH >="+pthresh+"p & FFDI >= "+pthresh+"p"
  resP2@gsnPanelLeft=0.65
  resP2@gsnPanelRight=0.9
  resP2@gsnPanelMainPosXF  = 0.77
 gsn_panel(wks,plots3,(/4,1/),resP2)

frame(wks)


