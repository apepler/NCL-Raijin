load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Cthresh=90
Fthresh=90
levels=fspan(-0.7,0.7,15)
years=ispan(1979,2016,1)
snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

; Get the index

inames=(/"SOI","DMI","SAM"/)
index2=new((/dimsizes(years),4,3/),"float")
istart=(/1876,1982,1979/)
idir="/short/eg3/asp561/Timeseries/"
ifiles=(/"SOI_18762016.csv","NOAA.DMI.monthly.csv","SAM_NOAA_19792016.csv"/)

do i=0,2
  iyears=ispan(istart(i),2016,1)
  tmp=asciiread(idir+ifiles(i),(/dimsizes(iyears),13/),"float")

  if(istart(i).le.1979) then 
    index=tmp(ind(tmp(:,0).ge.1979),:) 
  else
    index=new((/dimsizes(years),13/),"float")
    index(ind(years.ge.istart(i)),:)=tmp 
  end if

  do m=1,3
   index2(:,m,i)=dim_avg_n(index(:,seasons(m):(seasons(m)+2)),1)
   if(istart(i).gt.1979) then
    index2(ind(years.lt.istart(i)),m,i)=index2@_FillValue
   end if
  end do
  index2(1:(dimsizes(years)-1),0,i)=(index(0:(dimsizes(years)-2),12)+dim_sum_n(index(1:(dimsizes(years)-1),1:2),1))/3
  if(istart(i).gt.1979) then
   index2(ind(years.le.istart(i)),0,i)=index2@_FillValue
  end if

 delete([/tmp,index,iyears/])
end do

; Get the data

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon
CHp=ifile->CH_percentile

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI_percentile
tmp=ffile->FFDI

FFDI2=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDI@_FillValue,FFDI)
FFDI2=where(tmp.eq.0..and.FFDI2.gt.90,FFDI2@_FillValue,FFDI2)
delete(tmp)

a=dimsizes(CHp)
CH95p=new((/dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!2="lat"
CH95p&lat=lat
CH95p!3="lon"
CH95p&lon=lon

; Create seasonal count of days

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(y,m,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh.and.CHp(I,:,:).ge.Cthresh,0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(y,0,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh.and.CHp(I,:,:).ge.Cthresh,0)
  delete(I)
end do

; Do correlations

CHcor=new((/3,4,a(1),a(2)/),"float")
CHcor!2="lat"
CHcor&lat=lat
CHcor!3="lon"
CHcor&lon=lon

do x=0,2
do m=0,3
CHcor(x,m,:,:)=spcorr_n(CH95p(:,m,:,:),conform(CH95p(:,m,:,:),index2(:,m,x),0),0)
CHcor(x,m,:,:)=where(dim_avg_n(CH95p(:,m,:,:),0).lt.0.5,CHcor@_FillValue,CHcor(x,m,:,:))

end do
end do

; Make figures

  wks = gsn_open_wks("pdf","seasoncorr_many_CHaines"+Cthresh+"p_FFDI"+Fthresh+"p_19792016_raster")             ; send graphics to PNG file

;  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  gsn_define_colormap(wks,"MPL_RdBu")
  gsn_reverse_colormap(wks)

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
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

  opts_psl = True
     opts_psl@cnLineColor       = "Black"
     opts_psl@cnInfoLabelOn     = False
     opts_psl@cnLineLabelFontHeightF = 0.01
     opts_psl@cnLineLabelPerimOn = False
     opts_psl@gsnContourLineThicknessesScale = 2
     opts_psl@gsnDraw         = False
     opts_psl@gsnFrame        = False
     opts_psl@cnLevelSelectionMode = "ExplicitLevels"
     opts_psl@cnLevels             = (/ -100.,0.05,100./)
;     opts_psl@gsnContourPosLineDashPattern = 2


  plot = new(12,graphic)

  do m=0,3
   do x=0,2
    p=m*3+x
    res@gsnLeftString=snames(m)+": Correlation with "+inames(x)
   plot(p) = gsn_csm_contour_map(wks,CHcor(x,m,:,:), res)     ; create plot
   
   pval=rtest(CHcor(x,m,:,:),dimsizes(index2(:,m,x)),0)
   pval=where(dim_avg_n(CH95p(:,m,:,:),0).lt.0.5,pval@_FillValue,pval)
   pval!0="lat"
   pval&lat=lat
   pval!1="lon"
   pval&lon=lon
   contour_psl = gsn_csm_contour(wks,pval,opts_psl)
   overlay(plot(p), contour_psl)   
   delete([/pval,contour_psl/])
  end do
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPanelMainString = "Correlation of days with C-Haines>="+Cthresh+"p and FFDI>="+Fthresh+"p"
  resP@gsnPanelLabelBar    = True
  resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/4,3/),resP)



