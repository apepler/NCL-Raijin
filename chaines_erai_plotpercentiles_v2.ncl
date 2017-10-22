load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

date=20030118
years=ispan(1979,2016,1)


dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,2)
T=ind(time2.eq.date)
lat=ifile->lat
lon=ifile->lon
CH=ifile->CHaines


  wks = gsn_open_wks("pdf","chaines_erai_percentiles_19792016_v2")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@tiMainString         = "C-Haines index: "+date    ; main title
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels             = fspan(0.,18.,10)
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

  plot = new(6,graphic)

  pvals=(/50,99,90,99.7,95,99.95/)
  counts=round(pvals*dimsizes(CH&time)/100.,3)
  counts(3)=dimsizes(CH&time)-dimsizes(years)
  counts(5)=dimsizes(CH&time)-round(dimsizes(years)/5.,3)
  pnames=(/"Median","99th percentile","90th percentile","1 year ARI","95th percentile","5 year ARI"/)

  cmin=(/1,6,3,6,4,6/)

  do p=0,5
   res@tiMainString=pnames(p)+" of ERAI C-Haines index"
   res@cnLevels             = ispan(cmin(p),cmin(p)+9,1)

    a=dimsizes(CH)-1
    CHpc=new((/dimsizes(lat),dimsizes(lon)/),typeof(CH))
    CHpc!0="lat"
    CHpc&lat=lat
    CHpc!1="lon"
    CHpc&lon=lon

    do jlat=0,dimsizes(lat)-1
    do ilon=0,dimsizes(lon)-1
      CHsort=CH(ind(.not.ismissing(CH(:,jlat,ilon))),jlat,ilon)
      if(dimsizes(CHsort).ge.0.95*dimsizes(CH&time)) then
        qsort(CHsort)
        tmp=round(counts(p)*1.0*dimsizes(CHsort)/dimsizes(CH&time),3)
        CHpc(jlat,ilon)=(/CHsort(tmp)/)
      end if
      delete(CHsort)
    end do
    end do

   plot(p) = gsn_csm_contour_map(wks,CHpc, res)     ; create plot
  end do

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  gsn_panel(wks,plot,(/3,2/),resP)
