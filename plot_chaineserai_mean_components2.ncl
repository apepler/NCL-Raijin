load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; This is complex - use same threshold (95 p) for everything
; Calculate the threshold again for each dataset, it's the safest way

levels=fspan(-10,25,15)
years=ispan(1979,2016,1)
snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

; Get the data

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon

type=(/"700hPa Ta","850hPa Ta","850hPa Td","Ta diff (850 - 700)","850hPa Dewpoint Depression"/)
data=new((/dimsizes(time),dimsizes(lat),dimsizes(lon),5/),"float")
ffile=addfile(dir+"TaTd_19792016_0600UTC.nc","r")

data(:,:,:,0)=(/ffile->ta(:,1,:,:)/)-273.15
data(:,:,:,1)=(/ffile->ta(:,0,:,:)/)-273.15
data(:,:,:,2)=(/ffile->td(:,0,:,:)/)-273.15
data(:,:,:,3)=data(:,:,:,1)-data(:,:,:,0)
data(:,:,:,4)=data(:,:,:,1)-data(:,:,:,2)

CH3=new((/dimsizes(years),dimsizes(seasons),dimsizes(lat),dimsizes(lon),5/),"float")

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH3(y,m,:,:,:)=dim_avg_n(data(I,:,:,:),0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH3(y,0,:,:,:)=dim_avg_n(data(I,:,:,:),0)
  delete(I)
end do

CHave=dim_avg_n_Wrap(CH3,0)

CHave!1="lat"
CHave&lat=lat
CHave!2="lon"
CHave&lon=lon

; Make figures

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 9
wks = gsn_open_wks(wks_type,"chaines_mean_tatd_season")             ; send graphics to PNG file
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
 
  plot(p) = gsn_csm_contour_map(wks,CHave(m,:,:,x), res)     ; create plot

   end do
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
resP@gsnPaperOrientation = "portrait"
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="Seasonal mean (C)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
;  resP@gsnPanelMainString = "Percentage change in days >= "+pthresh+"p"
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/4,5/),resP)



