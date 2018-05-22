load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

edir="/g/data1/ub4/erai/netcdf/6hr/atmos/"
odir=""
date2=(/2012,6,5,0/)
date=date2(0)+sprinti("%2.2i",date2(1))+sprinti("%2.2i",date2(2))+sprinti("%2.2i",date2(3))

;hdir="/g/data/eg3/asp561/CycloneTracking/ERAI/proj100_highs_rad10cv0.075/"
hdir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/"
hname="austhighs_10deg_proj100_rad10cv0.075_"
;cdir="/g/data/eg3/asp561/CycloneTracking/ERAI/proj100_lows_rad5cv0.15/"
cdir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/"
cname="austlows_8deg_proj100_rad5cv0.15_"
fdir="/g/data/eg3/asp561/Fronts/"
fname="austfronts_len3_3deg_"
jdir="/g/data/eg3/asp561/CattoData_20052015/"
jname="storm_types_Aus333_"
; Think about cold fronts later

; Getting file with spaces is a pain
var="Storm types for use in collaboration with A. Dowdy & J. Catto:"
a=addfile(jdir+jname+date2(0)+".nc","r")
catto=a->$var$
catto&lat@units="degrees north"
catto&lon@units="degrees east"

time=(a->time)*6
time@units="hours since "+date2(0)+"-1-1 00:00:00"
time2=cd_calendar(time,0)
I1=ind(time2(:,1).eq.date2(1).and.time2(:,2).eq.date2(2).and.time2(:,3).eq.date2(3))

; Set up code to make a plot of Catto fronts for selected day
; In the lambert conformal projection

  wks_type = "pdf"
  wks_type@wkPaperWidthF = 6
  wks_type@wkPaperHeightF = 4

   res                 = True                    ; plot mods desired
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@gsnPaperOrientation = "portrait"
  res@gsnMaximize      = True
   res@cnLinesOn       = True
   res@gsnAddCyclic = False
   res@cnFillOn        = False               ; color plot desired
  res@lbLabelBarOn        = False
  res@mpFillOn              = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = True                ; Turn on lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks   
  res@tmXTOn = False
  res@tmYROn = False
  res@mpGridSpacingF = 10.
  res@mpProjection = "CylindricalEquidistant"
  res@mpGeophysicalLineThicknessF = 2.0
  res@mpMinLonF = 110
  res@lbOrientation     = "Vertical"
  res@mpMaxLonF = 160
  res@mpMinLatF = -50
  res@mpMaxLatF = -10
  res@mpGridLineThicknessF = 0.5
  res@lbLabelBarOn = False

levels=ispan(990,1030,4)
  res@cnLineColor       = "Black"
  res@cnInfoLabelOn     = False
  res@cnLineLabelFontHeightF = 0.01
  res@cnLineLabelPerimOn = False
  res@gsnContourLineThicknessesScale = 2.5
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = levels
  res@cnLineDrawOrder="Draw"
  res@cnLineLabelPlacementMode = "constant"
  res@cnLineDashSegLenF = 0.4
  res@cnLineLabelFontThicknessF = 2.
  res@gsnLeftString = "Simmonds method"

  nlevels = dimsizes(levels)
  flags   = new(nlevels,string)

;---Loop through each level and set the corresponding "flags" value.
  do i=0,nlevels-1
    if(mod(levels(i),4).ne.0) then
      flags(i) = "LineOnly"
    else
      flags(i) = "LineAndLabel"
    end if
  end do
  res@cnLevelFlags=flags

  res2=True
  res2@gsnDraw  = False                          ; don't draw
  res2@gsnFrame = False   
  res2@cnLevels             = 0.9
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 1., 0., 0., 0.25 /)/)
  res2@cnFillOn        = True
  res2@cnLinesOn  = False
  res2@gsnLeftString=""
  res2@gsnRightString=""
  res2@cnInfoLabelOn     = False
  res2@lbLabelBarOn = False
  res2@cnLineLabelsOn  = False 

   date3=date2(0)+"-"+date2(1)+"-"+date2(2)+" "+date2(3)+" UTC"
   wks = gsn_open_wks(wks_type,"WeatherSystems_"+date+"_v2")          ; send graphics to PNG file
    plots = new(2,graphic)

   res@gsnRightString=date3
  Data=addfile(systemfunc("ls -1 "+edir+"oper_an_sfc/v01/psl/psl_6hrs_ERAI_historical_an-sfc_"+date2(0)+sprinti("%2.2i",date2(1))+"*"),"r")
   time3=cd_calendar(Data->time,0)
   I2=ind(time3(:,2).eq.date2(2).and.time3(:,3).eq.date2(3))
  p=Data->psl(I2,:,:)
  p=p/100.
  plots(0) = gsn_csm_contour_map(wks,p,res)
  res@gsnLeftString = "Catto method"
  plots(1) = gsn_csm_contour_map(wks,p,res)

; Now, need to add the various systems
; To test, just start with my cyclones

  FILES1=addfile(cdir+cname+date2(0)+".nc","r")
  systems=FILES1->systems(I1,:,:)
  contour2 = gsn_csm_contour(wks,systems,res2)
  overlay(plots(0),contour2)

; Add my fronts
  delete(systems)
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 0., 0., 1., 0.25 /)/)
  FILES1=addfile(fdir+fname+date2(0)+".nc","r")
  systems=FILES1->systems(I1-1,:,:)
  contour3 = gsn_csm_contour(wks,systems,res2)
  overlay(plots(0),contour3)

; And my highs

  delete(systems)
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 1., 0.25, 0., 0.5 /)/)
  FILES1=addfile(hdir+hname+date2(0)+".nc","r")
  systems=FILES1->systems(I1,:,:)
  contour4 = gsn_csm_contour(wks,systems,res2)
  overlay(plots(0),contour4)

; Finally, Jen data
; I already have it loaded

;  cattof=where(catto(I1,:,:).eq.2.or.catto(I1,:,:).eq.4.or.catto(I1,:,:).eq.6.or.catto(I1,:,:).eq.7,1,0)
;  copy_VarMeta(catto(I1,:,:),cattof)
  FILES1=addfile(jdir+"cold_fronts_Aus"+date2(0)+"_v2.nc","r")
  cattof=FILES1->systems(I1,:,:)

  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 0., 0., 1., 0.25 /)/)
  contour5 = gsn_csm_contour(wks,cattof,res2)
  overlay(plots(1),contour5)

  cattoc=where(catto(I1,:,:).eq.1.or.catto(I1,:,:).eq.4.or.catto(I1,:,:).eq.5.or.catto(I1,:,:).eq.7,1,0)
  copy_VarMeta(catto(I1,:,:),cattoc)
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 1., 0., 0., 0.25 /)/)
  contour6 = gsn_csm_contour(wks,cattoc,res2)
  overlay(plots(1),contour6)

  cattot=where(catto(I1,:,:).eq.3.or.catto(I1,:,:).ge.5.,1,0)
  copy_VarMeta(catto(I1,:,:),cattot)
  res2@cnFillPalette = (/(/ 0., 0., 0., 0. /),(/ 1., 0.75, 0, 0.5 /)/)
  contour7 = gsn_csm_contour(wks,cattot,res2)
  overlay(plots(1),contour7)

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
resP@gsnPaperOrientation = "portrait"
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plots,(/1,2/),resP)


