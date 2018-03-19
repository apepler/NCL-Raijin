load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

ofile="ERAIcorr_MUCAPE_control_Australia_methodpanel_nonzero"

years=ispan(2000,2016,1)
lonlim=(/110,160/)
latlim=(/-45,-10/)
levels=(/0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95/)
hours=(/0,6,12,18/)

snames=(/"DJF","MAM","JJA","SON"/)
smonths=new((/4,2/),"integer")
smonths(:,0)=(/12,3,6,9/)-1
smonths(:,1)=(/2,5,8,11/)-1

; Setup - need to have the latitude and longitude
; Get CHaines to use for correlations

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

I=ind(lon.ge.min(lonlim).and.lon.le.max(lonlim))
J=ind(lat.ge.min(latlim).and.lat.le.max(latlim))

delete(ifile)

ifile=addfiles( systemfunc (" ls -1 " + cdir+ types(4) + "/mucape_2?????.nc"), "r")
    control=ifile[:]->mucape(:,J,I)
    time=cd_calendar(ifile[:]->time,0)


; Set up plot first

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
  res@mpLeftCornerLatF   = latlim(0)
  res@mpLeftCornerLonF   = lonlim(0)
  res@mpRightCornerLatF  = latlim(1)
  res@mpRightCornerLonF  = lonlim(1)
  res@mpGridLineDashPattern = 2
  res@mpGridLineThicknessF = 0.5
  res@mpProjection = "CylindricalEquidistant"
  res@cnLevels             = levels
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

  resP                  = True                   ; modify the panel plot
  resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="MUCAPE corr with control:"+tnames(4)
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
  resP@gsnPanelYWhiteSpacePercent = 5


  lablets=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)

  wks_type = "pdf"
  ;wks_type@wkOrientation = "landscape"
  wks_type@wkPaperWidthF = 10
  wks_type@wkPaperHeightF = 8

  hcorr=new((/4,4,dimsizes(J),dimsizes(I)/),"float")
  hcorr!0="Hour"
  hcorr!1="Type"
  hcorr!2="lat"
  hcorr&lat=lat(J)
  hcorr!3="lon"
  hcorr&lon=lon(I)
  
  scorr=hcorr
  scorr!0="Season"

  ; Calculate all the correlations

  do p=0,3
   print(types(p))
   ifile=addfiles( systemfunc (" ls -1 " + cdir+ types(p) + "/mucape_2?????.nc"), "r")
    mucape=ifile[:]->mucape(:,J,I)



    do h=0,3
     T=ind(time(:,3).eq.hours(h))
     hcorr(h,p,:,:)=(/spcorr_n(mask(control(T,:,:),control(T,:,:).gt.0,True),mask(mucape(T,:,:),mucape(T,:,:).gt.0,True),0)/)

;     hcorr(h,p,:,:)=(/spcorr_n(control(T,:,:),mucape(T,:,:),0)/)
     delete(T)
    end do

    do s=0,3
     if(smonths(s,1).ge.smonths(s,0)) then
       T=ind(time(:,1).ge.smonths(s,0).and.time(:,1).le.smonths(s,1))
     else
       T=ind(time(:,1).ge.smonths(s,0).or.time(:,1).le.smonths(s,1))
     end if
     scorr(s,p,:,:)=(/spcorr_n(mask(control(T,:,:),control(T,:,:).gt.0,True),mask(mucape(T,:,:),mucape(T,:,:).gt.0,True),0)/)
;     scorr(s,p,:,:)=(/spcorr_n(control(T,:,:),mucape(T,:,:),0)/)
     delete(T)
    end do

    delete(mucape)
  end do

  ; Make the hourly plots

  wks = gsn_open_wks(wks_type,ofile+"_vhour")             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  plot=new(16,graphic)
 
  do h=0,3
   do x=0,3
    p=h*4+x
    res@gsnLeftString=lablets(x*4+h)+") "+hours(h)+" UTC"
    res@gsnRightString=tnames(x)
    plot(p)=gsn_csm_contour_map_ce(wks,hcorr(h,x,:,:), res)     ; create plot
   end do
 end do

 gsn_panel(wks,plot,(/4,4/),resP)
 delete([/wks,plot/])

; And season plot

  wks = gsn_open_wks(wks_type,ofile+"_vseason")             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  plot=new(16,graphic)

  do s=0,3
   do x=0,3
    p=s*4+x
    res@gsnLeftString=lablets(x*4+s)+") "+snames(s)
    res@gsnRightString=tnames(x)
    plot(p)=gsn_csm_contour_map_ce(wks,scorr(s,x,:,:), res)     ; create plot
   end do
 end do

 gsn_panel(wks,plot,(/4,4/),resP)
 delete([/wks,plot/])


