load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Cthresh=6
Fthresh=25
levels=(/1,5,10,15,20,30,40,50,60,70,80/)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)

CHp=ifile->CHaines

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI
FFDI=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDI@_FillValue,FFDI)

years=ispan(1979,2016,1)

a=dimsizes(CHp)
CH95p=new((/3,dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!3="lat"
CH95p&lat=CHp&lat
CH95p!4="lon"
CH95p&lon=CHp&lon

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(0,y,m,:,:)=dim_num_n(CHp(I,:,:).ge.Cthresh,0)
  CH95p(1,y,m,:,:)=dim_num_n(FFDI(I,:,:).ge.Fthresh,0)
  CH95p(2,y,m,:,:)=dim_num_n(FFDI(I,:,:).ge.Fthresh.and.CHp(I,:,:).ge.Cthresh,0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(0,y,0,:,:)=dim_num_n(CHp(I,:,:).ge.Cthresh,0)
  CH95p(1,y,0,:,:)=dim_num_n(FFDI(I,:,:).ge.Fthresh,0)
  CH95p(2,y,0,:,:)=dim_num_n(FFDI(I,:,:).ge.Fthresh.and.CHp(I,:,:).ge.Cthresh,0)

  delete(I)
end do

CH95pave=dim_avg_n_Wrap(CH95p,0)

  wks = gsn_open_wks("pdf","seasonfreq_CHainesge"+Cthresh+"_FFDIge"+Fthresh+"_19792016_separated")             ; send graphics to PNG file

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
;  res@cnFillMode = "RasterFill"
;  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False


  plot = new(12,graphic)
  type=(/"C-Haines >="+Cthresh,"FFDI >="+Fthresh,"Both"/)
  plot = new(12,graphic)

  do m=0,3
   do x=0,2
    p=m*3+x

;  Make the array as we go
  tmp=dim_avg_n_Wrap(CH95p(x,:,m,:,:),0)
  tmp!0="lat"
  tmp&lat=CHp&lat
  tmp!1="lon"
  tmp&lon=CHp&lon

  res@gsnLeftString=snames(m)+": "+type(x)
  plot(p) = gsn_csm_contour_map(wks,tmp, res)     ; create plot
   end do
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPanelMainString = "Average frequency between 1979 and 2016"
  resP@gsnPanelLabelBar    = True
 resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/4,3/),resP)


