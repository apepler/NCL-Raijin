load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; This is complex - use same threshold (95 p) for everything
; Calculate the threshold again for each dataset, it's the safest way

pthresh=95
levels=ispan(-140,140,20)
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

type=(/"Stability Score (CA)","Humidity Score (CB)","C-Haines","FFDI","Pyroconvection Risk"/)
data=new((/4,dimsizes(time),dimsizes(lat),dimsizes(lon)/),"float")
data(0,:,:,:)=(/ifile->CA/)
data(1,:,:,:)=(/ifile->CB/)
data(2,:,:,:)=(/ifile->CHaines/)

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
data(3,:,:,:)=(/ffile->FFDI/)
data(3,:,:,:)=where(data(3,:,:,:).ge.999,data@_FillValue,data(3,:,:,:))

CHpc=new((/5,dimsizes(lat),dimsizes(lon)/),"float")

   do x=0,3
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


CH95p=new((/5,dimsizes(years),4,dimsizes(lat),dimsizes(lon)/),"integer")
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
  do x=0,3
    CH95p(x,y,m,:,:)=dim_num_n(data(x,I,:,:).gt.conform_dims(dimsizes(data(x,I,:,:)),CHpc(x,:,:),(/1,2/)),0)
  end do
  CH95p(4,y,m,:,:)=dim_num_n(data(3,I,:,:).gt.conform_dims(dimsizes(data(3,I,:,:)),CHpc(3,:,:),(/1,2/)).and.data(2,I,:,:).gt.conform_dims(dimsizes(data(2,I,:,:)),CHpc(2,:,:),(/1,2/)),0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  do x=0,3
    CH95p(x,y,0,:,:)=dim_num_n(data(x,I,:,:).gt.conform_dims(dimsizes(data(x,I,:,:)),CHpc(x,:,:),(/1,2/)),0)
  end do
  CH95p(4,y,0,:,:)=dim_num_n(data(3,I,:,:).gt.conform_dims(dimsizes(data(3,I,:,:)),CHpc(3,:,:),(/1,2/)).and.data(2,I,:,:).gt.conform_dims(dimsizes(data(2,I,:,:)),CHpc(2,:,:),(/1,2/)),0)

  delete(I)
end do

; Make figures

wks_type = "eps"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 9
wks = gsn_open_wks(wks_type,"Fig4_chaines_change_p"+pthresh+"_sig_forsub")             ; send graphics to PNG file
  gsn_define_colormap(wks,"MPL_RdBu")
  gsn_reverse_colormap(wks)

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
  res@gsnLeftString=""
  res@gsnRightString=""
  res@cnFillMode = "RasterFill"
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

  res@gsnStringFontHeightF = 0.025
  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02

  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "PostDraw"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"


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
     opts_psl@cnLineDrawOrder="Draw"
;     opts_psl@gsnContourPosLineDashPattern = 2

lablets=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)

  plot = new(20,graphic)

  do m=0,3
   do x=0,4
    p=m*5+x

;  Make the array as we go
  tmp1=dim_avg_n_Wrap(CH95p(x,I1,m,:,:),0)
  tmp2=dim_avg_n_Wrap(CH95p(x,I2,m,:,:),0)
  tmp=100*((tmp2/where(tmp1.eq.0,CH95p@_FillValue,tmp1))-1) 
  tmp=where(dim_avg_n_Wrap(CH95p(x,:,m,:,:),0).lt.0.5,CH95p@_FillValue,tmp)
  tmp!0="lat"
  tmp&lat=lat
  tmp!1="lon"
  tmp&lon=lon

  res@gsnLeftString=lablets(x*4+m)+") "+snames(m)+" "+type(x)

  if(m.gt.2) then
    res@tmXBLabelsOn   = True
  else
    res@tmXBLabelsOn   = False
  end if

  if(x.gt.0) then
;   res@tmYLLabelsOn   = False
   res@tmYLLabelFontColor = "white" 
  else
;   res@tmYLLabelsOn   = True
   res@tmYLLabelFontColor = "black" 
  end if
 

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
resP@gsnPaperOrientation = "portrait"
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="Long-term change (%)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
;  resP@gsnPanelMainString = "Percentage change in days >= "+pthresh+"p"
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/4,5/),resP)



