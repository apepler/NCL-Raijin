load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

Cthresh=95
Fthresh=95
levels=ispan(-100,100,10)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)

CHp=ifile->CH_percentile

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI_percentile

FFDI2=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDI@_FillValue,FFDI)


years=ispan(1979,2016,1)

a=dimsizes(CHp)
CH95p=new((/dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!2="lat"
CH95p&lat=CHp&lat
CH95p!3="lon"
CH95p&lon=CHp&lon

FFDI95p=CH95p

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(y,m,:,:)=dim_num_n(CHp(I,:,:).ge.Cthresh,0)
  FFDI95p(y,m,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh,0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(y,0,:,:)=dim_num_n(CHp(I,:,:).ge.Cthresh,0)
  FFDI95p(y,0,:,:)=dim_num_n(FFDI2(I,:,:).ge.Fthresh,0)
  delete(I)
end do

CHave=new((/3,4,a(1),a(2)/),"float")
CHave!2="lat"
CHave&lat=CHp&lat
CHave!3="lon"
CHave&lon=CHp&lon
FFDIave=CHave


I1=ind(years.le.1997)
I2=ind(years.ge.1998)
CHave(0,:,:,:)=dim_avg_n_Wrap(CH95p(I1,:,:,:),0)
CHave(1,:,:,:)=dim_avg_n_Wrap(CH95p(I2,:,:,:),0)
CHave(2,:,:,:)=100*((CHave(1,:,:,:)/where(CHave(0,:,:,:).eq.0,CHave@_FillValue,CHave(0,:,:,:)))-1)

FFDIave(0,:,:,:)=dim_avg_n_Wrap(FFDI95p(I1,:,:,:),0)
FFDIave(1,:,:,:)=dim_avg_n_Wrap(FFDI95p(I2,:,:,:),0)
FFDIave(2,:,:,:)=100*((FFDIave(1,:,:,:)/where(FFDIave(0,:,:,:).eq.0,FFDIave@_FillValue,FFDIave(0,:,:,:)))-1)


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

  wks = gsn_open_wks("pdf","seasonfreq_CHaines"+Cthresh+"p_FFDI"+Fthresh+"p_19791997_19982016_change_separate")             ; send graphics to PNG file
  gsn_define_colormap(wks,"MPL_RdBu")
  gsn_reverse_colormap(wks)

  plot = new(8,graphic)

  do p=0,3
   res@tiMainString="% change in "+snames(p)+" days with C-Haines>="+Cthresh+"p"
   plot(p*2) = gsn_csm_contour_map(wks,CHave(2,p,:,:), res)     ; create plot

   res@tiMainString="% change in "+snames(p)+" days with FFDI>="+Fthresh+"p"

   plot(p*2+1) = gsn_csm_contour_map(wks,FFDIave(2,p,:,:), res)     ; create plot   
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  gsn_panel(wks,plot,(/4,2/),resP)


