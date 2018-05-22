load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

edir="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_sfc/v01/psl/"
dir="/short/eg3/asp561/cts.dir/gcyc_out/"

reanal="ERAI"
proj="proj100"
thresh="rad10cv0.075"
type="highs"
hdir=dir+reanal+"/"+proj+"_"+type+"_"+thresh+"/"

year=2017
daylist=(/"0601","0607","0617","0902","0915"/) ; Dates I want
;numlist=ispan(100,1400,100)
;numlist=180

; Load in highs
filename=hdir+"tracks_"+year+".dat"
  read_data= asciiread(filename,-1,"string")
  header    = read_data(0)        ; Header. Use for variable names.
  data      = read_data(1:)       ; Get rid of first line which is a header.
  quote = inttochar(34)
  delim=",: "+quote
  date1=tointeger(str_get_field(data,3,delim))
  time1= tointeger(str_get_field(data,4,delim))
  I=ind(date1/10000.eq.max(date1/10000).and.time1.eq.0)
  tmp = str_split_by_length(str_get_field(data(I),3,delim),(/4,4/))
  hdate=tmp(:,1)
  hlat = tofloat(str_get_field(data(I),8,delim))
  hlon= tofloat(str_get_field(data(I),7,delim))
  delete([/date1,time1,tmp,I/])

; ERAI slp
  efiles=addfiles(systemfunc (" ls -1 " + edir + "psl_6hrs_ERAI_historical_an-sfc_" + year + "*.nc"),"r")
  etime=cd_calendar(efiles[:]->time,0)
  edate=etime(:,1)*100+etime(:,2)

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
   res@cnLineLabelsOn  = True              ; turn off contour lines
;  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels             = (/0.9,1.2/)
;  res@cnFillPalette = (/"white","blue","white"/)
  res@mpFillOn              = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = True                ; Turn on lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks   
  res@tmXTOn = False
  res@tmYROn = False
   res@mpLimitMode       = "Corners"            ; choose range of map
   res@mpLeftCornerLatF  = -45.
   res@mpLeftCornerLonF  = 85.
   res@mpRightCornerLatF = -5.
   res@mpRightCornerLonF = 170.
;  res@cnFillMode = "RasterFill"   
   res@mpProjection        = "LambertConformal"
   res@mpLambertParallel1F = -18
   res@mpLambertParallel2F = -36
   res@mpLambertMeridianF  = 136
   res@gsnLeftString="Highs"
  res@cnLevelSelectionMode           = "ExplicitLevels"
  res@cnLevels                       = ispan(960,1040,4)
  res@cnLineLabelBackgroundColor     = -1    ; transparent
  res@cnLineThicknessF               = 2.5
  res@cnInfoLabelOn              = False
  res@cnLineLabelPlacementMode = "Constant"
  res@cnLineDashSegLenF = 0.25

  slp_res                                = True
  slp_res@gsnFrame       = False
  slp_res@gsnDraw        = False
  slp_res@gsnLeftString  = ""
  slp_res@gsnRightString = ""
  slp_res@cnLineColor                    = "blue"
  slp_res@cnLevelSelectionMode           = "ExplicitLevels"
  slp_res@cnLevels                       = (/0.9,1.2/)
  slp_res@cnLineLabelBackgroundColor     = -1    ; transparent
  slp_res@cnLineThicknessF               = 5
  slp_res@cnInfoLabelOn              = False
  slp_res@cnLineLabelsOn  = False

   pres=True
   pres@gsMarkerColor = "blue"
   pres@gsMarkerSizeF=0.02
   pres@gsMarkerThicknessF=5
   pres@gsMarkerIndex="x"

do p=0,dimsizes(daylist)-1
   date=tointeger(year+daylist(p))
   n=ind(edate.eq.tointeger(daylist(p)).and.etime(:,3).eq.0)
   wks = gsn_open_wks(wks_type,"ERAI_highs_"+date+"")          ; send graphics to PNG file
   res@gsnRightString=date

; Add ERAI slp to plot
   psl=efiles[:]->psl(n,:,:)
   slp=psl/100
   copy_VarMeta(psl,slp)

   plot = gsn_csm_contour_map(wks,slp,res)

; Add Irina lines to plot
   I=ind(hdate.eq.daylist(p))
   dum=gsn_add_polymarker(wks,plot,hlon(I),hlat(I),pres)
 draw(plot)

; Add legend to plot somehow?

  frame(wks)
  delete([/wks,plot,dum,I/])
end do
