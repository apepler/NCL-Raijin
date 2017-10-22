load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Pthresh=90
levels=ispan(-100,100,10)
years=ispan(1979,2016,1)
snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)

CH=ifile->CHaines
lat=ifile->lat
lon=ifile->lon

CH2=new((/dimsizes(time),dimsizes(lat),dimsizes(lon),3/),typeof(CH))
CH2(:,:,:,0)=CH
CH2(:,:,:,1)=ifile->CA
CH2(:,:,:,2)=ifile->CB

; Need to recalculate the 90th percentile for each

pctile=new((/dimsizes(lat),dimsizes(lon),3/),"float")
do jlat=0,dimsizes(lat)-1
do ilon=0,dimsizes(lon)-1
do k=0,2
  I1=ind(.not.ismissing(CH2(:,jlat,ilon,k)))

  if(dimsizes(I1).gt.(0.9*dimsizes(time)).and.max(CH2(:,jlat,ilon,k)).gt.0) then
    CHsort=CH2(I1,jlat,ilon,k)
    qsort(CHsort)
    pctile(jlat,ilon,k)=(/CHsort(round(Pthresh*dimsizes(CHsort)/100,3))/)
  end if
  delete([/I1,CHsort/])
end do
end do
end do


CH3=new((/dimsizes(years),dimsizes(seasons),dimsizes(lat),dimsizes(lon),3/),"integer")

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH3(y,m,:,:,:)=dim_num_n(CH2(I,:,:,:).ge.conform_dims(dimsizes(CH2(I,:,:,:)),pctile,(/1,2,3/)),0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH3(y,0,:,:,:)=dim_num_n(CH2(I,:,:,:).ge.conform_dims(dimsizes(CH2(I,:,:,:)),pctile,(/1,2,3/)),0)
  delete(I)
end do

CHave=new((/3,4,dimsizes(lat),dimsizes(lon),3/),"float")

I1=ind(years.le.1997)
I2=ind(years.ge.1998)
CHave(0,:,:,:,:)=dim_avg_n_Wrap(CH3(I1,:,:,:,:),0)
CHave(1,:,:,:,:)=dim_avg_n_Wrap(CH3(I2,:,:,:,:),0)
CHave(2,:,:,:,:)=100*((CHave(1,:,:,:,:)/where(CHave(0,:,:,:,:).eq.0,CHave@_FillValue,CHave(0,:,:,:,:)))-1)

CHave!2="lat"
CHave&lat=lat
CHave!3="lon"
CHave&lon=lon


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
  res@mpLeftCornerLatF   = -45.
  res@mpLeftCornerLonF   = 110.
  res@mpRightCornerLatF  = -10.
  res@mpRightCornerLonF  = 155.
  res@mpGridLineDashPattern = 2
  res@mpGridLineThicknessF = 0.5
  res@mpProjection = "CylindricalEquidistant"
  res@cnLevels             = levels 
  res@cnLinesOn=False

  wks = gsn_open_wks("pdf","seasonfreq_CHaines"+Pthresh+"p_19791997_19982016_change_components")             ; send graphics to PNG file
  gsn_define_colormap(wks,"MPL_RdBu")
  gsn_reverse_colormap(wks)

  type=(/"C-Haines","Stability term","Moisture term"/)
  plot = new(12,graphic)

  do p1=0,3
   do p2=0,2
    p=p1*3+p2
    res@tiMainString=snames(p1)+" days "+type(p2)+">="+Pthresh+"p"
   plot(p) = gsn_csm_contour_map(wks,CHave(2,p1,:,:,p2), res)     ; create plot

  end do
  end do
  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  gsn_panel(wks,plot,(/4,3/),resP)


