load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

levels=ispan(-30,30,5)
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


CH3=new((/dimsizes(years),dimsizes(seasons),dimsizes(lat),dimsizes(lon),3/),"float")

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH3(y,m,:,:,:)=dim_avg_n(CH2(I,:,:,:),0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH3(y,0,:,:,:)=dim_avg_n(CH2(I,:,:,:),0)
  delete(I)
end do

I1=ind(years.le.1997)
I2=ind(years.ge.1998)

CH3!0="year"
CH3!2="lat"
CH3&lat=lat
CH3!3="lon"
CH3&lon=lon


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
  res@gsnRightString="%"
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


  wks = gsn_open_wks("pdf","seasonfreq_CHaines_mean_19791997_19982016_change_components_sig")             ; send graphics to PNG file
  gsn_define_colormap(wks,"MPL_RdBu")
  gsn_reverse_colormap(wks)

  type=(/"C-Haines","Stability term","Moisture term"/)
  plot = new(12,graphic)

  do m=0,3
   do x=0,2
    p=m*3+x

;  Make the array as we go
  tmp1=dim_avg_n_Wrap(CH3(I1,m,:,:,x),0)
  tmp2=dim_avg_n_Wrap(CH3(I2,m,:,:,x),0)
  tmp=100*((tmp2/where(tmp1.eq.0,CH3@_FillValue,tmp1))-1)
  tmp!0="lat"
  tmp&lat=lat
  tmp!1="lon"
  tmp&lon=lon

  res@gsnLeftString=snames(m)+": "+type(x)
  plot(p) = gsn_csm_contour_map(wks,tmp, res)     ; create plot

; And also the significance
  r1a=CH3(I1,m,:,:,x)
  r2a=CH3(I2,m,:,:,x)
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


