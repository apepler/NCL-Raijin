load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Cthresh=90
Fthresh=90
levels=ispan(-100,100,10)
years=ispan(1979,2016,1)
I1=ind(years.le.1997)
I2=ind(years.ge.1998)

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
CH95p=new((/3,dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!0="variable"
CH95p!1="year"
CH95p!2="season"
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
  CH95p(0,y,m,:,:)=dim_num_n(CHp(I,:,:).ge.Cthresh,0)
  CH95p(1,y,m,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh,0)
  CH95p(2,y,m,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh.and.CHp(I,:,:).ge.Cthresh,0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(0,y,0,:,:)=dim_num_n(CHp(I,:,:).ge.Cthresh,0)
  CH95p(1,y,0,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh,0)
  CH95p(2,y,0,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh.and.CHp(I,:,:).ge.Cthresh,0)
  delete(I)
end do

; Make figures

wks = gsn_open_wks("pdf","seasonfreq_CHaines"+Cthresh+"p_FFDI"+Fthresh+"p_19791997_19982016_percentchange_sig")             ; send graphics to PNG file
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
  res@gsnRightString="%"
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
  type=(/"C-Haines >="+Cthresh+"p","FFDI >="+Fthresh+"p","Both"/)
  plot = new(12,graphic)

  do m=0,3
   do x=0,2
    p=m*3+x

;  Make the array as we go
  tmp1=dim_avg_n_Wrap(CH95p(x,I1,m,:,:),0)
  tmp2=dim_avg_n_Wrap(CH95p(x,I2,m,:,:),0)
  tmp=100*((tmp2/where(tmp1.eq.0,CH95p@_FillValue,tmp1))-1) 
  tmp=where(dim_avg_n_Wrap(CH95p(x,:,m,:,:),0).lt.0.5,CH95p@_FillValue,tmp)
  tmp!0="lat"
  tmp&lat=lat
  tmp!1="lon"
  tmp&lon=lon

  res@gsnLeftString=snames(m)+": "+type(x)
  plot(p) = gsn_csm_contour_map(wks,tmp, res)     ; create plot

; And also the significance
  r1a=CH95p(x,I1,m,:,:)
  r2a=CH95p(x,I2,m,:,:)
  xAve = dim_avg_n(r1a,0)              ; calculate means at each grid point 
  yAve = dim_avg_n(r2a,0)
  xVar = dim_variance_n(r1a,0)         ; calculate variances
  yVar = dim_variance_n(r2a,0)
  sigr = 0.05                        ; critical sig lvl for r
  xEqv = equiv_sample_size (r1a(lat|:,lon|:,year|:), sigr,0)
  yEqv = equiv_sample_size (r2a(lat|:,lon|:,year|:), sigr,0)                               
;  xN   = wgt_areaave (xEqv, 1., 1., 0)    ; wgty could be gaussian weights 
;  yN   = wgt_areaave (yEqv, 1., 1., 0) 
  iflag= False                        ; population variance similar
  pval = ttest(xAve,xVar,xEqv, yAve,yVar,yEqv, iflag, False) 
  pval=where(dim_avg_n_Wrap(CH95p(x,:,m,:,:),0).lt.0.5,pval@_FillValue,pval)
  pval!0="lat"
  pval!1="lon"
  pval&lat = lat
  pval&lon = lon

  contour_psl = gsn_csm_contour(wks,pval,opts_psl)    
  overlay(plot(p), contour_psl)

   delete([/tmp,tmp1,tmp2,r1a,r2a/])
  delete([/xAve,yAve,xVar,yVar,sigr,xEqv,yEqv,pval,contour_psl/])

   end do
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPanelMainString = "Percentage change between 1979-1997 and 1998-2016"
  resP@gsnPanelLabelBar    = True
 resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/4,3/),resP)



