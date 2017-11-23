load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

levels=fspan(-0.7,0.7,15)
years=ispan(1979,2016,1)

; Get the index

iname="SOI"
iyears=ispan(1876,2016,1)
;ifile="/short/eg3/asp561/Timeseries/SAM_NOAA_19792016.csv"
;ifile="/short/eg3/asp561/Timeseries/NOAA.DMI.monthly.csv"
ifile="/short/eg3/asp561/Timeseries/SOI_18762016.csv"
tmp=asciiread(ifile,(/dimsizes(iyears),13/),"float")

I=ind(tmp(:,0).ge.1979.and.tmp(:,0).le.2016)
index=tmp(I,0:12)
delete([/tmp,I,ifile/])

;read_data= asciiread(ifile,-1,"string")
;header    = read_data(0)        ; Header. Use for variable names.
;data      = read_data(1:)       ; Get rid of first line which is a header.

; Get the data

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon
CHp=ifile->CHaines

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI_percentile
tmp=ffile->FFDI

FFDI2=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDI@_FillValue,tmp)
FFDI2=where(tmp.eq.0..and.FFDI.gt.90,FFDI2@_FillValue,FFDI2)
delete(tmp)

a=dimsizes(CHp)
CH95p=new((/2,dimsizes(years),4,a(1),a(2)/),"float")
CH95p!3="lat"
CH95p&lat=lat
CH95p!4="lon"
CH95p&lon=lon

; Create seasonal count of days

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)
index2=new((/dimsizes(years),4/),"float")

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(0,y,m,:,:)=dim_avg_n(CHp(I,:,:),0)
  CH95p(1,y,m,:,:)=dim_avg_n(FFDI2(I,:,:),0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(0,y,0,:,:)=dim_avg_n(CHp(I,:,:),0)
  CH95p(1,y,0,:,:)=dim_avg_n(FFDI2(I,:,:),0)
  delete(I)
end do

do m=1,3
   index2(:,m)=dim_avg_n(index(:,seasons(m):(seasons(m)+2)),1)
end do
   index2(1:(dimsizes(years)-1),0)=(index(0:(dimsizes(years)-2),12)+dim_sum_n(index(1:(dimsizes(years)-1),1:2),1))/3

; Do correlations

CHcor=new((/2,4,a(1),a(2)/),"float")
CHcor!2="lat"
CHcor&lat=lat
CHcor!3="lon"
CHcor&lon=lon

do x=0,1
do m=0,3
CHcor(x,m,:,:)=spcorr_n(CH95p(x,:,m,:,:),conform(CH95p(x,:,m,:,:),index2(:,m),0),0)
end do
end do

; Make figures

  wks = gsn_open_wks("pdf","seasoncorr_"+iname+"_CHaines_FFDI_19792016_raster")             ; send graphics to PNG file

;  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  gsn_define_colormap(wks,"MPL_RdBu")

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
  res@gsnLeftString=""
  res@gsnRightString=""
  res@cnFillMode = "RasterFill"
  res@cnLinesOn = False

  plot = new(8,graphic)
  type=(/"Mean C-Haines","Mean FFDI"/)
  plot = new(8,graphic)

  do m=0,3
   do x=0,1
    p=m*2+x
    res@gsnLeftString=snames(m)+": "+type(x)
   plot(p) = gsn_csm_contour_map(wks,CHcor(x,m,:,:), res)     ; create plot

  end do
  end do
  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPanelMainString = "Correlation with seasonal "+iname
  gsn_panel(wks,plot,(/4,2/),resP)



