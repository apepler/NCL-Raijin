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

type=(/"Stability Score (CA)","Humidity Score (CB)","C-Haines","FFDI","Compound Events"/)
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

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 8
wks = gsn_open_wks(wks_type,"Fig4_chaines_change_p"+pthresh+"_sig_ForRevision")             ; send graphics to PNG file
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

  tmres          = True
  tmres@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  tmres@tmYLLabelFontHeightF = 0.02
  if(x.eq.0) then
    tmres@tmYLValues = ispan(-40,-10,10)
  end if
  if(m.eq.3) then
    tmres@tmXBValues = ispan(110,150,10)
  end if
  
  plot(p) = gsn_csm_contour_map(wks,tmp, res)     ; create plot
  plot(p) = add_map_tickmarks2(wks,plot(p),tmres)
  delete(tmres)

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


end
