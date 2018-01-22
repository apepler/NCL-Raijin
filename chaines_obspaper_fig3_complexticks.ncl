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

Cthresh=8
levels1=ispan(0,80,10)
levels1a=ispan(0,90,10)

CthreshP=95
levels2=ispan(0,20,2)
levels2a=ispan(0,22,2)
Fthresh=95

levels3=ispan(0,10,1)
levels3a=ispan(0,11,1)

levels=(/1,5,10,15,20,30,40,50,60,70,80/)

dir="/g/data/eg3/asp561/CHaines_ERAI/"
ifile=addfile(dir+"CHaines_19792016_06UTC.nc","r")
time=ifile->time
time2=cd_calendar(time,0)

CH=ifile->CHaines
CHp=ifile->CH_percentile

ffile=addfile(dir+"FFDI_AWAP_19792016_regrid.nc","r")
FFDI=ffile->FFDI
FFDIp=ffile->FFDI_percentile

FFDIp=where(conform_dims(dimsizes(FFDI),dim_max_n(FFDI,0).eq.dim_min_n(FFDI,0),(/1,2/)),FFDIp@_FillValue,FFDIp)

years=ispan(1979,2016,1)

a=dimsizes(CHp)
CH95p=new((/3,dimsizes(years),4,a(1),a(2)/),"integer")
CH95p!3="lat"
CH95p&lat=CHp&lat
CH95p!4="lon"
CH95p&lon=CHp&lon

snames=(/"DJF","MAM","JJA","SON"/)
seasons=(/12,3,6,9/)

do y=0,dimsizes(years)-1
do m=1,3
  I=ind(time2(:,0).eq.years(y).and.time2(:,1).ge.seasons(m).and.time2(:,1).le.seasons(m)+2)
  CH95p(0,y,m,:,:)=dim_num_n(CH(I,:,:).gt.Cthresh,0)
  CH95p(1,y,m,:,:)=dim_num_n(CHp(I,:,:).gt.CthreshP,0)
  CH95p(2,y,m,:,:)=dim_num_n(FFDIp(I,:,:).gt.Fthresh.and.CHp(I,:,:).gt.CthreshP,0)
  delete(I)
end do
end do

do y=1,dimsizes(years)-1
  I=ind((time2(:,0).eq.years(y-1).and.time2(:,1).eq.12).or.(time2(:,0).eq.years(y).and.time2(:,1).le.2))
  CH95p(0,y,0,:,:)=dim_num_n(CH(I,:,:).gt.Cthresh,0)
  CH95p(1,y,0,:,:)=dim_num_n(CHp(I,:,:).gt.CthreshP,0)
  CH95p(2,y,0,:,:)=dim_num_n(FFDIp(I,:,:).gt.Fthresh.and.CHp(I,:,:).gt.CthreshP,0)

  delete(I)
end do

CH95pave=dim_avg_n_Wrap(CH95p,0)
wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 6
wks_type@wkPaperHeightF = 5

  wks = gsn_open_wks(wks_type,"Fig3_CHainesge"+Cthresh+"_"+CthreshP+"p_FFDIge"+Fthresh+"_ForRevision")             ; send graphics to PNG file

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
;  res@cnFillMode = "RasterFill"
;  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False
;  res@gsnMaximize=True
  res@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
  res@cnLineDrawOrder="PreDraw"
  res@mpDataSetName               = "Earth..4"
  res@mpFillDrawOrder             = "Draw"
  res@mpLandFillColor             = "transparent"
  res@mpOceanFillColor            = "white"
  res@gsnStringFontHeightF = 0.035
  res@gsnRightString=""
  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02    
;  res@cnLabelBarEndStyle="IncludeMinMaxLabels"
  plots1 = new(4,graphic)
  plots2 = new(4,graphic)
  plots3 = new(4,graphic)

res2=res
res@cnLevels             = levels1
res2@cnLevels             = levels2
;res2@tmYLLabelFontColor = "white"

labs1=(/"a","b","c","d"/)
labs2=(/"e","f","g","h"/)
labs3=(/"i","j","k","l"/)

  do m=0,3

;  Make the array as we go
  tmp=dim_avg_n_Wrap(CH95p(0,:,m,:,:),0)
  tmp!0="lat"
  tmp&lat=CHp&lat
  tmp!1="lon"
  tmp&lon=CHp&lon

;  if(m.lt.3) then
;   res@tmXBLabelFontColor = "white"
;   res2@tmXBLabelFontColor = "white"
;  else
;   res@tmXBLabelFontColor = "black"
;   res2@tmXBLabelFontColor = "black"
;  end if

  res@gsnLeftString=labs1(m)+") "+snames(m)
  plots1(m) = gsn_csm_contour_map(wks,tmp, res)     ; create plot

  tmres          = True
  tmres@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  tmres@tmYLLabelFontHeightF = 0.02
  tmres@tmYLValues = ispan(-40,-10,10)
  tmres@tmXBValues = ispan(110,150,10)
  if(m.lt.3) then
   tmres@tmXBLabelFontColor = "white"
  else
   tmres@tmXBLabelFontColor = "black"
  end if

  plots1(m) = add_map_tickmarks2(wks,plots1(m),tmres)
  delete(tmres@tmYLValues)

  tmp=(/dim_avg_n_Wrap(CH95p(1,:,m,:,:),0)/)
  res2@gsnLeftString=labs2(m)+") "+snames(m)
  plots2(m) = gsn_csm_contour_map(wks,tmp, res2)     ; create plot
  plots2(m) = add_map_tickmarks2(wks,plots2(m),tmres)
  tmp=(/dim_avg_n_Wrap(CH95p(2,:,m,:,:),0)/)

  res2@cnLevels             = levels3
  res2@gsnLeftString=labs3(m)+") "+snames(m)
  plots3(m) = gsn_csm_contour_map(wks,tmp, res2)     ; create plot
  plots3(m) = add_map_tickmarks2(wks,plots3(m),tmres)
  res2@cnLevels             = levels2
  delete(tmres)
  end do

  resP                  = True                   ; modify the panel plot
;  resP@gsnPanelMainString = "CH >= "+Cthresh
  resP@gsnMaximize=True
  resP@gsnPaperOrientation = "portrait"
  
  resP@lbTitleString="CH > "+Cthresh+" (days)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
;  resP@cnLabelBarEndStyle="IncludeMinMaxLabels"
;  resP@lbLabelStrings=levels1a
  resP@pmLabelBarOrthogonalPosF = -0.02

  resP@gsnPanelLabelBar    = True
;  resP@gsnPanelYWhiteSpacePercent = 1
  resP@amJust   = "TopLeft"
  resP@gsnPanelLeft=0.15
  resP@gsnPanelRight=0.4
 resP@gsnPanelBottom=0.05
 resP@gsnPanelTop=1.
  resP@gsnFrame = False ; don't advance frame yet
  resP@gsnPanelMainPosXF  = 0.28
  resP@gsnPanelFigureStringsFontHeightF = 0.012
  resP@gsnPanelMainFontHeightF = 0.018
  resP@pmLabelBarWidthF = 0.2
  resP@lbLabelFontHeightF = 0.012
  resP@lbTitleFontHeightF = 0.012
  gsn_panel(wks,plots1,(/4,1/),resP)

  resP2=resP
  resP2@lbTitleString="CH > "+CthreshP+"th percentile (days)"
;  delete(resP2@lbLabelStrings)
;  resP2@lbLabelStrings=levels2a
  resP2@gsnPanelLeft=0.4
  resP2@gsnPanelRight=0.65
  resP2@gsnPanelMainPosXF  = 0.52
 gsn_panel(wks,plots2,(/4,1/),resP2)

  resP2=resP

  resP2@lbTitleString= "Compound events (days)"
  resP2@gsnPanelLeft=0.65
  resP2@gsnPanelRight=0.9
  resP2@gsnPanelMainPosXF  = 0.77
 gsn_panel(wks,plots3,(/4,1/),resP2)

frame(wks)
;maximize_output(wks,True)

end
