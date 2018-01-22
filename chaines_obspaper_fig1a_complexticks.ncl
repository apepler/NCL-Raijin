load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("add_map_tickmarks2")
function add_map_tickmarks2(wks,plot,res)
local res2, bres, vpx, vpy, vpw, vph, xndc, yndc, npts, n, j, nlat, \
nlon, delta, bot_lon, top_lon, lft_lat, rgt_lat, xblabels, xbvalues, \
xtlabels, xtvalues, yllabels, ylvalues, yrlabels, yrvalues, xfix, \
xlat, xlon, yfix, annoid, anno_str
begin
;---Make a copy of the original resource list.
  res2 = res

;---Retrieve edges of plot in NDC space.
  getvalues plot
     "vpXF"      :  vpx
     "vpYF"      :  vpy
     "vpWidthF"  :  vpw
     "vpHeightF" :  vph
  end getvalues

;---Turn off tickmarks associated with map. We want to add our own.
  setvalues plot
    "pmTickMarkDisplayMode" : "Never"
  end setvalues

;---Initialize resources for tickmark plot. User shouldn't change these.
  bres                          = True
  bres@vpXF                     = vpx
  bres@vpYF                     = vpy
  bres@vpWidthF                 = vpw
  bres@vpHeightF                = vph
  bres@trXMinF                  = vpx
  bres@trXMaxF                  = vpx + vpw
  bres@trYMinF                  = vpy - vph
  bres@trYMaxF                  = vpy
  bres@tmEqualizeXYSizes        = True

;---This resource the user can change in main code if desired.
  bres@gsnTickMarksPointOutward = get_res_value(res2,"gsnTickMarksPointOutward",True)

;
; NDC Points to scan on X and Y axes. These arrays will be used to
; find the closest NDC pair that gets us close to the location where
; we want a tickmark.
;
  npts = 100000   ; Increase to get closer match for tickmarks
  xndc = fspan(vpx,vpx+vpw,npts)
  yndc = fspan(vpy-vph,vpy,npts)

  n    = dimsizes(yndc) 
  xfix = new(n,float)
  yfix = new(n,float)
  xlon = new(n,float)
  xlat = new(n,float)
  delta = 0.001

;---Left axis tickmarks
  if(isatt(res2,"tmYLValues")) then
    lft_lat    = get_res_value(res2,"tmYLValues",-1)
    nlat       = dimsizes(lft_lat)
    ylvalues = new(nlat,float)
    yllabels = new(nlat,string)
    xfix  = vpx + 0.0001 ; Just a smidge into the plot to make sure we don't
                         ; get missing values returned.
;
; Loop across each left latitude value that we want a tickmark for,
; and try to find the closest X,Y NDC coordinate pair along this axis.
;
    NhlNDCToData(plot,xfix,yndc,xlon,xlat)
    do j=0,dimsizes(lft_lat)-1
      NhlNDCToData(plot,xfix,yndc,xlon,xlat)
      ii = minind(fabs(xlat-lft_lat(j)))
      if(.not.any(ismissing(ii)).and.fabs(xlat(ii)-lft_lat(j)).le.delta)
        yllabels(j) = fabs(lft_lat(j)) + ""
        ylvalues(j) = yndc(ii(0))
        if(lft_lat(j).lt.0) then
          yllabels(j) = yllabels(j) + "~S~o~N~S"
        end if
        if(lft_lat(j).gt.0) then
          yllabels(j) = yllabels(j) + "~S~o~N~N"
        end if
      end if
      delete(ii)
    end do
    bres@tmYLMode   = "Explicit"
    bres@tmYLValues = ylvalues
    bres@tmYLLabels = get_res_value(res2,"tmYLLabels",yllabels)
    bres@tmYLLabelAngleF = 90
  else
    bres@tmYLOn       = False
    bres@tmYLLabelsOn = False
  end if

;---Right axis tickmarks
  if(isatt(res2,"tmYRValues")) then
    rgt_lat    = get_res_value(res2,"tmYRValues",-1)
    nlat       = dimsizes(rgt_lat)
    yrvalues = new(nlat,float)
    yrlabels = new(nlat,string)

    xfix  = vpx + vpw - 0.0001 ; Just a smidge into the plot to make sure we don't
                               ; get missing values returned.
;
; Loop across each right latitude value that we want a tickmark for,
; and try to find the closest X,Y NDC coordinate pair along this axis.
;
    do j=0,dimsizes(rgt_lat)-1
      NhlNDCToData(plot,xfix,yndc,xlon,xlat)
      ii = minind(fabs(xlat-rgt_lat(j)))
      if(.not.any(ismissing(ii)).and.fabs(xlat(ii)-rgt_lat(j)).le.delta)
        yrlabels(j) = fabs(rgt_lat(j)) + ""
        yrvalues(j) = yndc(ii(0))
        if(rgt_lat(j).lt.0) then
          yrlabels(j) = yrlabels(j) + "~S~o~N~S"
        end if
        if(rgt_lat(j).gt.0) then
          yrlabels(j) = yrlabels(j) + "~S~o~N~N"
        end if
      end if
      delete(ii)
    end do
    bres@tmYROn       = True
    bres@tmYRLabelsOn = True
    bres@tmYUseLeft   = False
    bres@tmYRMode     = "Explicit"
    bres@tmYRValues   = yrvalues
    bres@tmYRLabels   = get_res_value(res2,"tmYRLabels",yrlabels)
  else
    bres@tmYUseLeft   = False
    bres@tmYROn       = False
    bres@tmYRLabelsOn = False
  end if

;---Top axis tickmarks
  if(isatt(res2,"tmXTValues")) then
    top_lon    = get_res_value(res2,"tmXTValues",-1)
    nlon       = dimsizes(top_lon)
    xtvalues = new(nlon,float)
    xtlabels = new(nlon,string)

    yfix  = vpy - 0.0001 ; Just a smidge into the plot to make sure we don't
                         ; get missing values returned.
;
; Loop across each top longitude value that we want a tickmark for,
; and try to find the closest X,Y NDC coordinate pair along this axis.
;
    do j=0,dimsizes(top_lon)-1
      NhlNDCToData(plot,xndc,yfix,xlon,xlat)
      ii = minind(fabs(xlon-top_lon(j)))
      if(.not.any(ismissing(ii)).and.fabs(xlon(ii)-top_lon(j)).le.delta)
        xtlabels(j) = fabs(top_lon(j)) + ""
        xtvalues(j) = xndc(ii(0))
        if(top_lon(j).lt.0) then
          xtlabels(j) = xtlabels(j) + "~S~o~N~W"
        end if
        if(top_lon(j).gt.0) then
          xtlabels(j) = xtlabels(j) + "~S~o~N~E"
        end if
      end if
      delete(ii)
    end do
    bres@tmXTOn       = True
    bres@tmXTLabelsOn = True
    bres@tmXUseBottom = False
    bres@tmXTMode     = "Explicit"
    bres@tmXTValues   = xtvalues
    bres@tmXTLabels   = get_res_value(res2,"tmXTLabels",xtlabels)
  else
    bres@tmXUseBottom = False
    bres@tmXTOn       = False
    bres@tmXTLabelsOn = False
  end if

;---Bottom axis tickmarks
  if(isatt(res2,"tmXBValues")) then
    bot_lon    = get_res_value(res2,"tmXBValues",-1)
    nlon       = dimsizes(bot_lon)
    xbvalues = new(nlon,float)
    xblabels = new(nlon,string)

    yfix  = vpy-vph + 0.0001 ; Just a smidge into the plot to make sure 
                             ; we don't get missing values returned.
;
; Loop across each bottom longitude value that we want a tickmark for,
; and try to find the closest X,Y NDC coordinate pair along this axis.
;
    do j=0,dimsizes(bot_lon)-1
      NhlNDCToData(plot,xndc,yfix,xlon,xlat)
      ii = minind(fabs(xlon-bot_lon(j)))
      if(.not.any(ismissing(ii)).and.fabs(xlon(ii)-bot_lon(j)).le.delta)
        xblabels(j) = fabs(bot_lon(j)) + ""
        xbvalues(j) = xndc(ii(0))
        if(bot_lon(j).lt.0) then
          xblabels(j) = xblabels(j) + "~S~o~N~W"
        end if
        if(bot_lon(j).gt.0) then
          xblabels(j) = xblabels(j) + "~S~o~N~E"
        end if
      end if
      delete(ii)
    end do
    bres@tmXBMode   = "Explicit"
    bres@tmXBValues = xbvalues
    bres@tmXBLabels = get_res_value(res2,"tmXBLabels",xblabels)
  else
    bres@tmXBOn       = False
    bres@tmXBLabelsOn = False
  end if

;
; Now that we are done figuring out where to put tickmarks, and
; what labels to use, get any "tm" resources that might have been
; set by the user, and create a blank plot with thes new tickmarks.
;

;---Get rest of user resources that were set with "tm".
  bres = get_res_eq(res2,"tm")

  bres = True   ; Above call will set bres to True if no "tm" resources, so
                ; make sure it is True still.

  bres@gsnDraw  = False
  bres@gsnFrame = False

;
; Create blank plot with new tickmarks (don't use gsn_csm_blank_plot, 
; because it wants to scale the size of your X and Y axes.)
;
  blank  = gsn_blank_plot(wks,bres)

;
; Attach new tickmarks to original plot. This will allow resizing
; if desired. The default is to attach one plot to the center of
; the other one. These two plots are already the same size.
;
  annoid = gsn_add_annotation(plot,blank,False)

;
; Be sure to return the annotation id, otherwise the
; tickmarks will disappear.
;
  anno_str = unique_string("annoid")
  plot@$anno_str$ = annoid

  return(plot)
end

begin

date=20030118
years=ispan(1979,2016,1)
pthresh=95
levels=ispan(1,13,1)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,2)
T=ind(time2.eq.date)
lat=ifile->lat
lon=ifile->lon
CH=ifile->CHaines

mfile=addfile("/g/data1/ub4/erai/netcdf/fx/ei_invariant_075x075_90N0E90S35925E.nc","r")
lat2=mfile->latitude
lon2=mfile->longitude
I=ind(lon2.ge.min(lon).and.lon2.le.max(lon))
J=ind(lat2.ge.min(lat).and.lat2.le.max(lat))
landmask=short2flt(mfile->lsm(0,J,I))

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 5.1
wks_type@wkPaperHeightF = 4.6

  wks = gsn_open_wks(wks_type,"Fig1a_chaines_erai_p"+pthresh+"_ForRevision")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@gsnMaximize=True
  res@gsnPaperMargin = 0.2
  res@gsnPaperOrientation = "portrait" 
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels
  res@gsnAddCyclic = False
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpFillOn              = True
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "Draw"           ; Draw map outlines last
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
  res@gsnRightString=""
;  res@lbLabelBarOn        = False
  res@cnLineLabelPlacementMode = "constant"
  res@cnLineLabelsOn = True
  res@cnLineDashSegLenF = 0.3
  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@lbTitleOn        =  True                ; turn on title
  res@lbTitlePosition  = "Bottom"              ; title position
;  res@lbTitleFontHeightF= .015                ; make title smaller
  res@lbTitleDirection = "Across"             ; title direction
;  res@lbLabelPosition  = "Center"             ; label position
  res@cnLabelBarEndStyle="IncludeMinMaxLabels"
  res@lbLabelStrings=ispan(0,14,1)
  res@tmXBLabelFontHeightF = 0.015                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.015
  res@lbTitleFontHeightF = 0.015
  res@tmYLLabelAngleF = 90.
  res@gsnMajorLatSpacing = 10  

;  res@tmYLMode   = "Explicit"
;  res@tmYLValues = ispan(-40,-10,10)
;  res@tmYLLabels = (/"40~S~o~N~S","30~S~o~N~S",\
;                     "20~S~o~N~S","10~S~o~N~S"/)

res@pmLabelBarOrthogonalPosF = .15
res@pmLabelBarWidthF = 0.55
res@pmLabelBarHeightF = 0.05

  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "Draw"

;  res@mpAreaMaskingOn       = True
;  res@mpGridMaskMode        = "MaskLand"
;  res@mpMaskAreaSpecifiers = "water"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"
;
;  res@mpFillAreaSpecifiers        = (/"water",       "land" /)
;  res@mpSpecifiedFillColors       = (/"white","white"/)


                                               ; last.
  nlevels = dimsizes(levels)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(levels(i).lt.6.or.mod(levels(i),2).ne.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

  counts=round(pthresh*dimsizes(CH&time)/100.,3)

;   res@gsnLeftString="a) "+pthresh+"th percentile of CH based on ERAI reanalysis"
    res@lbTitleString=pthresh+"~S~th~N~ percentile of CH based on reanalysis"
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
        tmp=round(counts*1.0*dimsizes(CHsort)/dimsizes(CH&time),3)
        CHpc(jlat,ilon)=(/CHsort(tmp)/)
      end if
      delete(CHsort)
    end do
    end do

    g6=where(CHpc.ge.6,1,0)
    print(sum(g6*landmask)/sum(landmask))
 
;   CHpc=mask(CHpc,landmask,0)

   plot=gsn_csm_contour_map_ce(wks,CHpc, res)     ; create plot

tmres          = True

;---Set the values you want for the map tickmarks.
;  tmres@tmXBLabelFontHeightF = 0.008         ; Default is a bit large

;---We're putting tickmarks on all four axes.
  tmres@tmYLValues = ispan(-40,-20,10)
  tmres@tmXBValues = ispan(110,150,10)
  res@tmYLLabelAngleF = 90.

;---Attach the new map tickmarks.
  map = add_map_tickmarks2(wks,plot,tmres)

   ptnames=(/"Alice Springs","Townsville","Perth","Wagga Wagga"/)
   ptlats=(/-23.7951,-19.2483,-31.9275,-35.1583/)
   ptlons=(/133.8890,146.7661,115.9764,147.4575/)

   gsres = True
   gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
   gsres@gsMarkerSizeF = 0.015

   txres               = True
   txres@txFontHeightF = 0.025
   txres@txJust="CenterRight"
   txres@txFontThicknessF = 3.   

  dum1=gsn_add_polymarker(wks,plot,ptlons,ptlats,gsres)
  dum2=gsn_add_text(wks,plot,ptnames,ptlons-1,ptlats,txres)

   txres2               = True
   txres2@txFontHeightF = 0.035
;   dum3=gsn_add_text(wks,plot,"a)",113,-13,txres2)
   gsn_text_ndc(wks,"a)",0.03,0.95,txres2)
   draw(plot)
   maximize_output(wks,True)
end
