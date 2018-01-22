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

levels1=(/0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
levels2=(/0.5,1.,1.5,2.,2.5/)

levels1a=(/0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1/)
levels2a=(/0,0.5,1.,1.5,2.,2.5,3./)

pthresh=95

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)
lat=ifile->lat
lon=ifile->lon

data=new((/2,dimsizes(time),dimsizes(lat),dimsizes(lon)/),"float")
data(0,:,:,:)=(/ifile->CHaines/)

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
data(1,:,:,:)=(/ffile->FFDI/)
data(1,:,:,:)=where(data(1,:,:,:).ge.999,data@_FillValue,data(1,:,:,:))

CHcor=new((/2,dimsizes(lat),dimsizes(lon)/),"float")
CHcor!1="lat"
CHcor&lat=lat
CHcor!2="lon"
CHcor&lon=lon

  CHcor(0,:,:)=spcorr_n(data(0,:,:,:),data(1,:,:,:),0)

CHpc=new((/2,dimsizes(lat),dimsizes(lon)/),"float")

   do x=0,1
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

atmp=dim_num_n(data(0,:,:,:).gt.conform_dims(dimsizes(data(0,:,:,:)),CHpc(0,:,:),(/1,2/)).and.data(1,:,:,:).gt.conform_dims(dimsizes(data(1,:,:,:)),CHpc(1,:,:),(/1,2/)),0)
btmp=dim_num_n(data(1,:,:,:).gt.conform_dims(dimsizes(data(1,:,:,:)),CHpc(1,:,:),(/1,2/)),0)

;CHcor(1,:,:)=(/100*atmp/where(btmp.eq.0,CHcor@_FillValue,btmp)/)
CHcor(1,:,:)=(/100.0*atmp/conform_dims(dimsizes(atmp),dimsizes(time),-1)/)
  CHcor(1,:,:)=where(CHcor(0,:,:).eq.CHcor@_FillValue,CHcor@_FillValue,CHcor(1,:,:)) 

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10.2
wks_type@wkPaperHeightF = 4.6
  wks = gsn_open_wks(wks_type,"Fig2_FFDICHcorr_propF"+pthresh+"pC"+pthresh+"p_ForRevision")             ; send graphics to PNG file

  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels             = levels
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
;  res@gsnLeftStringFontHeightF = 0.03
  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02
res@pmLabelBarOrthogonalPosF = .15
  res@cnLineLabelPlacementMode = "constant"
  res@cnLineLabelsOn = True
  res@cnLineDashSegLenF = 0.3

  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "PostDraw"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"
 res@pmLabelBarWidthF = 0.6

; Label bar title

  res@lbTitleOn        =  True                ; turn on title
  res@lbTitlePosition  = "Bottom"              ; title position
;  res@lbTitleFontHeightF= .015                ; make title smaller
  res@lbTitleDirection = "Across"             ; title direction
;  res@lbLabelPosition  = "Center"             ; label position
  res@cnLabelBarEndStyle="IncludeMinMaxLabels"


  plot = new(2,graphic)

   res@cnLevels=levels1
  nlevels = dimsizes(levels1)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(mod(i,2).ne.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

;   res@gsnLeftString="a) Correlation coefficient (r)"
  res@lbTitleString="Correlation coefficient (r)"
  res@lbLabelStrings=levels1a
   plot(0) = gsn_csm_contour_map(wks,CHcor(0,:,:), res)     ; create plot

tmres          = True
  tmres@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  tmres@tmYLLabelFontHeightF = 0.02
  tmres@tmYLValues = ispan(-40,-20,10)
  tmres@tmXBValues = ispan(110,150,10)
  plot(0) = add_map_tickmarks2(wks,plot(0),tmres)


    delete([/res@cnLevels,flags,res@cnLevelFlags,res@lbLabelStrings/])

  nlevels = dimsizes(levels2)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(mod(i,2).eq.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

;   res@tmYLLabelFontColor = "white"   
   res@cnLevels=levels2
  res@lbLabelStrings=levels2a
  res@lbTitleString="Compound event days (%)"

   plot(1) = gsn_csm_contour_map(wks,CHcor(1,:,:), res)     ; create plot
  delete(tmres@tmYLValues) 
  plot(1) = add_map_tickmarks2(wks,plot(1),tmres)

txres2               = True
   txres2@txFontHeightF = 0.025
gsn_text_ndc(wks,"a)",0.015,0.75,txres2)
gsn_text_ndc(wks,"b)",0.515,0.75,txres2)

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize=True
  resP@gsnPaperOrientation = "portrait"
  resP@amJust   = "TopLeft"
;  resP@gsnPanelFigureStrings = (/"a)","b)"/)
  resP@gsnPanelFigureStringsPerimOn = False
;  resP@gsnPanelBottom = 0.05
  gsn_panel(wks,plot,(/1,2/),resP)

end
