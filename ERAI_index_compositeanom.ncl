load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

rdir="/g/data/eg3/asp561/ERAI/"
years=ispan(1979,2016,1)
sname="MJJA"
mlist=(/5,8/)
q=0.9

iname="STRI"
ifile="/short/eg3/asp561/Timeseries/STRI_18902016.csv"
iyears=ispan(1890,2016,1)

ofile="ERAIcomposite_"+iname+"_Q"+q*100+"_"+sname+"_v2"

; Step 1 - load the index, and restrict to 1979-2016

  tmp=asciiread(ifile,(/dimsizes(iyears),13/),"float")

  if(min(iyears).le.1979) then
    index=tmp(ind(tmp(:,0).ge.1979.and.tmp(:,0).le.2016),:)
  else
    index=new((/dimsizes(years),13/),"float")
    index(ind(years.ge.min(iyears)),:)=tmp(ind(tmp(:,0).le.2016),:)
  end if

  delete(tmp)

; Step 2 - identify when index is above and below quantile, for the seasons of interest

 index2=new(dimsizes(index),typeof(index))
 index2=0
 index2(:,0)=index(:,0)

 do m=1,12
  tmp=index(:,m)
  qsort(tmp)
  q1=tmp(round(q*dimsizes(tmp),3))
  q2=tmp(round((1-q)*dimsizes(tmp),3))
  
  index2(:,m)=where(index(:,m).ge.q1,1,where(index(:,m).le.q2,-1,0))
 end do 

 index3=reshape(index2(:,1:12),(/12*dimsizes(years),1/))

; Step two - load the ERAI u/v/slp anomalies

a=addfile(rdir+"ERAI.slp.monanom.19792016.nc","r")
lon=tofloat(a->lon)
lat=tofloat(a->lat)
time=cd_calendar(a->time,0)
slp=a->psl
slp=slp/100

a=addfile(rdir+"ERAI.uv.850hPa.monanom.19792016.nc","r")
u=a->ua(:,0,:,:)
v=a->va(:,0,:,:)

a=addfile(rdir+"ERAI.tas.monanom.19792016.nc","r")
T2=a->tas

; Set up all the plotting options

  res = True
  res@mpFillOn = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is Lowres
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@gsnDraw         = False
  res@gsnFrame        = False
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res@mpGridLineDashPattern = 2 ; lat/lon lines as dashed
  res@mpGridAndLimbOn = False
  res@mpGridSpacingF = 10.
  res@mpProjection = "CylindricalEquidistant"
  res@mpGeophysicalLineThicknessF = 2.0
  res@mpMinLonF = 90
  res@lbOrientation     = "Vertical"
  res@mpMaxLonF = 180
  res@mpMinLatF = -50
  res@mpMaxLatF = 0
  res@mpGridLineThicknessF = 0.5
  res@cnFillOn        = True
  res@cnLinesOn       = False
;  res@lbLabelBarOn        = False
  res@gsnStringFontHeightF = 0.02
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -0.5                ; set min contour level
  res@cnMaxLevelValF       =  0.5              ; set max contour level
  res@cnLevelSpacingF      =  0.1              ; set contour spacing
  res@gsnLeftString = ""
  res@gsnRightString=sname
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

  vres                         = True             ; plot mods desired
  vres@vcRefLengthF            = 0.05            ; define length of vec ref
  vres@vcMinDistanceF          = 0.017
  vres@gsnDraw         = False
  vres@gsnFrame        = False
  vres@vcLineArrowThicknessF = 1.1
  vres@vcRefMagnitudeF           = 3.0             ; define vector ref mag
  vres@tiMainFontHeightF = 0.02
  vres@gsnLeftString = ""
  vres@gsnRightString=""
  vres@vcLineArrowColor="black"
  vres@vcVectorDrawOrder         = "PostDraw"
  vres@vcRefAnnoOrthogonalPosF = -1.1   ; move ref vector into plot
  vres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  vres@vcRefAnnoString1 = "$VMG$ m/s"
  vres@vcRefAnnoString2On = False

opts_psl = True
     opts_psl@cnLineColor       = "Black"
     opts_psl@cnInfoLabelOn     = False
     opts_psl@cnLineLabelFontHeightF = 0.01
     opts_psl@cnLineLabelPerimOn = False
     opts_psl@gsnContourLineThicknessesScale = 3
     opts_psl@gsnDraw         = False
     opts_psl@gsnFrame        = False
     opts_psl@cnLevelSelectionMode = "ExplicitLevels"
     opts_psl@cnLevels             = fspan(-10,10,21)
     opts_psl@cnLineDrawOrder="Draw"
     opts_psl@cnLineLabelPlacementMode = "constant"
;     opts_psl@cnLineDashSegLenF = 0.4
     opts_psl@cnLineLabelFontThicknessF = 2.
  opts_psl@gsnLeftString = ""
  opts_psl@gsnRightString=""
  opts_psl@gsnContourNegLineDashPattern = 1

; Step three - for both cases, make mean anomaly for months of choice. 
; For slp maybe also show where sig diff to zero?

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF =10
wks_type@wkPaperHeightF = 5
wks = gsn_open_wks(wks_type,ofile)
  gsn_define_colormap(wks,"MPL_RdBu")
  gsn_reverse_colormap(wks)
  plots=new(2,graphic)
;pattern = "$_$_$_$_$_"
;   opts_psl@gsnContourNegLineDashPattern =NhlNewDashPattern(wks,pattern)

I=ind(time(:,1).ge.mlist(0).and.time(:,1).le.mlist(1).and.index3(:,0).eq.-1)
res@gsnLeftString = iname+" <= "+100*(1-q)+"th percentile"
plots(0)=gsn_csm_contour_map(wks,dim_avg_n_Wrap(T2(I,:,:),0),res)

contour_psl=gsn_csm_contour(wks,dim_avg_n_Wrap(slp(I,:,:),0),opts_psl)
overlay(plots(0),contour_psl)

wind = gsn_csm_vector(wks,dim_avg_n_Wrap(u(I,:,:),0),dim_avg_n_Wrap(v(I,:,:),0),vres)
overlay(plots(0),wind)
delete(I)

I=ind(time(:,1).ge.mlist(0).and.time(:,1).le.mlist(1).and.index3(:,0).eq.1)
res@gsnLeftString = iname+" >= "+100*q+"th percentile"
plots(1)=gsn_csm_contour_map(wks,dim_avg_n_Wrap(T2(I,:,:),0),res)

contour_psl=gsn_csm_contour(wks,dim_avg_n_Wrap(slp(I,:,:),0),opts_psl)
overlay(plots(1),contour_psl)

wind = gsn_csm_vector(wks,dim_avg_n_Wrap(u(I,:,:),0),dim_avg_n_Wrap(v(I,:,:),0),vres)
overlay(plots(1),wind)

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPaperOrientation = "portrait"
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelFigureStrings= (/"a)","b)"/) ; add strings to panel
  resP@amJust   = "TopLeft"
  resP@lbOrientation       = "vertical"
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="(~S~o~N~C)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@gsnPanelLabelBar    = True
  resP@lbBottomMarginF = 0.3
  resP@lbLabelFontHeightF =.01
  gsn_panel(wks,plots,(/1,2/),resP)

