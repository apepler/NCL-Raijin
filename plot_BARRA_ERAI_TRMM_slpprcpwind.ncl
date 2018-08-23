load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;load "./WRFUserARW.ncl"

bdir="/g/data/ma05/BARRA_R/v1/"
edir="/g/data1/ub4/erai/netcdf/6hr/atmos/"
tdir="/g/data/eg3/asp561/TRMM/"
odir=""
date2=(/2013,1,29,0/)
date=date2(0)+sprinti("%2.2i",date2(1))+sprinti("%2.2i",date2(2))+sprinti("%2.2i",date2(3))

; Set up plotting resources

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
  res@mpMinLonF = 145
  res@lbOrientation     = "Vertical"
  res@mpMaxLonF = 160
  res@mpMinLatF = -40
  res@mpMaxLatF = -25
  res@mpGridLineThicknessF = 0.5
  res@cnFillOn        = True
  res@cnLinesOn       = False
;  res@lbLabelBarOn        = False
  res@gsnStringFontHeightF = 0.02
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/ 0.1, 0.5, 1., 2., 5., 10., 20., 30., 40., 50. /)


  res@gsnLeftString = ""
  res@gsnRightString=""
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

  vres                         = True             ; plot mods desired
  vres@vcRefLengthF            = 0.05            ; define length of vec ref
  vres@vcMinDistanceF          = 0.017
  vres@gsnDraw         = False
  vres@gsnFrame        = False
  vres@vcLineArrowThicknessF = 1.1
  vres@vcRefMagnitudeF           = 10.0             ; define vector ref mag
  vres@tiMainFontHeightF = 0.02
  vres@gsnLeftString = ""
  vres@gsnRightString=""
  vres@vcLineArrowColor="black"
  vres@vcVectorDrawOrder         = "PostDraw"
  vres@vcRefAnnoOrthogonalPosF = -1   ; move ref vector into plot
  vres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  vres@vcRefAnnoString1 = "$VMG$ m/s"
  vres@vcRefAnnoString2On = False

levels=ispan(990,1030,2)
opts_psl = True
     opts_psl@cnLineColor       = "Black"
     opts_psl@cnInfoLabelOn     = False
     opts_psl@cnLineLabelFontHeightF = 0.01
     opts_psl@cnLineLabelPerimOn = False
     opts_psl@gsnContourLineThicknessesScale = 2.5
     opts_psl@gsnDraw         = False
     opts_psl@gsnFrame        = False
     opts_psl@cnLevelSelectionMode = "ExplicitLevels"
     opts_psl@cnLevels             = levels
     opts_psl@cnLineDrawOrder="Draw"
     opts_psl@cnLineLabelPlacementMode = "constant"
     opts_psl@cnLineDashSegLenF = 0.4
     opts_psl@cnLineLabelFontThicknessF = 2.
  opts_psl@gsnLeftString = ""
  opts_psl@gsnRightString=""

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

; Set up panel

wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 5

wks = gsn_open_wks(wks_type,odir+"ERAI_BARRA_TRMM_850wind_prcp_mslp_"+date)
  gsn_define_colormap(wks,"precip_11lev")
  plots = new(3,graphic)

; ERAI plotting

  Data=addfile(systemfunc("ls -1 "+edir+"oper_an_pl/v01/ua/ua_6hrs_ERAI_historical_an-pl_"+date2(0)+sprinti("%2.2i",date2(1))+"*"),"r")
  time=Data->time
  time2=cd_calendar(time,0)
  I=ind(time2(:,2).eq.date2(2).and.time2(:,3).eq.date2(3))
  lev=Data->lev
  J=ind(lev.eq.85000.)
  u=Data->ua(I,J,:,:)
  delete(Data)

  Data=addfile(systemfunc("ls -1 "+edir+"oper_an_pl/v01/va/va_6hrs_ERAI_historical_an-pl_"+date2(0)+sprinti("%2.2i",date2(1))+"*"),"r")
  v=Data->va(I,J,:,:)
  delete(Data)

  Data=addfile(systemfunc("ls -1 "+edir+"oper_an_sfc/v01/psl/psl_6hrs_ERAI_historical_an-sfc_"+date2(0)+sprinti("%2.2i",date2(1))+"*"),"r")
  p=Data->psl(I,:,:)
  p=p/100.
  delete([/Data,I,J/])

  ; What file is precip? And call it T!
  delete([/time,time2/])
  Data=addfile(systemfunc("ls -1 /g/data1/ub4/erai/netcdf/3hr/atmos/oper_fc_sfc/v01/tp/tp_3hrs_ERAI_historical_fc-sfc_"+date2(0)+sprinti("%2.2i",date2(1))+"*"),"r")
  time=Data->time
  time2=cd_calendar(time,0)
  I=ind(time2(:,2).eq.date2(2).and.time2(:,3).eq.date2(3))
  prcp=dim_sum_n_Wrap(Data->tp(I:(I+1),:,:),0)*1000.

  copy_VarMeta(p,prcp)
  delete([/Data,time,time2,I/])

  res@gsnLeftString = "ERAI"
  plots(0) = gsn_csm_contour_map(wks,prcp,res)
  contour_psl = gsn_csm_contour(wks,p,opts_psl)
  wind = gsn_csm_vector(wks,u,v,vres)

  overlay(plots(0),wind)
  overlay(plots(0),contour_psl)

  delete([/prcp,p,u,v,lev/])

; BARRA plotting
  res@gsnAddCyclic=False
  vres@vcMinDistanceF = 0.03

  date3=date2(0)+sprinti("%2.2i",date2(1))+sprinti("%2.2i",date2(2))+"T"+sprinti("%2.2i",date2(3))+"00Z"

  Data=addfile(bdir+"analysis/prs/wnd_ucmp/"+date2(0)+"/"+sprinti("%2.2i",date2(1))+"/wnd_ucmp-an-prs-PT0H-BARRA_R-v1-"+date3+".nc","r")
  lev=Data->pressure
  J=closest_val(850.,lev)
  u=Data->wnd_ucmp(J,:,:)
  delete(Data)

  Data=addfile(bdir+"analysis/prs/wnd_vcmp/"+date2(0)+"/"+sprinti("%2.2i",date2(1))+"/wnd_vcmp-an-prs-PT0H-BARRA_R-v1-"+date3+".nc","r")
  v=Data->wnd_vcmp(J,:,:)
  delete(Data)

  Data=addfile(bdir+"analysis/spec/mslp/"+date2(0)+"/"+sprinti("%2.2i",date2(1))+"/mslp-an-spec-PT0H-BARRA_R-v1-"+date3+".nc","r")
  p=Data->mslp
  p=p/100.
  delete(Data)

  Data=addfile(bdir+"/forecast/spec_proc/accum_prcp/"+date2(0)+"/"+sprinti("%2.2i",date2(1))+"/accum_prcp-fc-spec_proc-PT1H-BARRA_R-v1-"+date3+".nc","r")
  prcp=Data->accum_prcp(2,:,:)
  delete(Data)

  res@gsnLeftString = "BARRA"
  plots(1) = gsn_csm_contour_map(wks,prcp,res)

  contour_psl = gsn_csm_contour(wks,p,opts_psl)
  wind = gsn_csm_vector(wks,u,v,vres)

  overlay(plots(1),wind)
  overlay(plots(1),contour_psl)

; TRMM plotting
  delete(prcp)
  date4=date2(0)+sprinti("%2.2i",date2(1))+sprinti("%2.2i",date2(2))
  if(date2(0).ge.2000.and.date2(0).le.2010) then
    dfile=tdir+date2(0)+"/3B42."+date4+".7A.nc"
    if (fileexists(dfile)) then
      FILE1=addfile(dfile,"r")
      prcp=FILE1->pcp(date2(3)/3,:,:)
    end if
  else if(date2(0).gt.2013) then
      dfile=tdir+date2(0)+"/3B42."+date4+"."+sprinti("%2.2i",date2(3))+".7.nc"
      if (fileexists(dfile)) then
       FILE1=addfile(dfile,"r")
       prcp=FILE1->pcp(0,:,:)
      end if
  else
    dfile=tdir+date2(0)+"/3B42."+date4+".7.nc"
    if (fileexists(dfile)) then
     FILE1=addfile(dfile,"r")
     prcp=FILE1->pcp(date2(3)/3,:,:)
    end if
  end if
  end if
; print(dfile)
  prcp=prcp*6
  res@gsnLeftString = "TRMM (instantaneous)"
  plots(2) = gsn_csm_contour_map(wks,prcp,res)

; Plot panel

  resP                  = True                   ; modify the panel plot
  resP@gsnMaximize      = True                   ; use full page
resP@gsnPaperOrientation = "portrait"
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="6h rainfall accumulation (mm)"
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plots,(/1,3/),resP)

