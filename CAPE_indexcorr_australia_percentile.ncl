load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

cdir="/g/data/eg3/ajd548/NESP_EWP/CAPE/"
types=(/"erainterim_adiabatic_largemix","erainterim_adiabatic_smallmix","erainterim_adiabatic_nomix","erainterim_reversible_largemix","erainterim_reversible_smallmix","erainterim_reversible_nomix"/)
tnames=(/"Adiabatic, Large Mixing","Adiabatic, Small Mixing","Adiabatic, No Mixing","Reversible, Large Mixing","Reversible,  Small Mixing","Reversible, No Mixing"/)

levels=fspan(-0.7,0.7,15)
years=ispan(1979,2016,1)
snames=(/"MJJASO","NDJFMA","DJF","MAM","JJA","SON"/)
seasons1=(/5,11,12,3,6,9/)
seasons2=(/10,4,2,5,8,11/)

; Get the index

inames=(/"SOI","DMI","SAM"/)
istart=(/1876,1982,1979/)
idir="/short/eg3/asp561/Timeseries/"
ifiles=(/"SOI_18762016.csv","NOAA.DMI.monthly.csv","SAM_NOAA_19792016.csv"/)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

lonlim=(/110,160/)
latlim=(/-45,-10/)

; Setup - need to have the latitude and longitude
; Just do mean CAPE for each month to start

ifile=addfile(cdir+types(0)+"/mucape_200001.nc","r")
lat=ifile->latitude
lon=ifile->longitude

ifile=addfile("mucape_95pc_19792016_monthly.nc","r")
mucape=ifile->mucape

do i=0,dimsizes(inames)-1
  iyears=ispan(istart(i),2016,1)
  tmp=asciiread(idir+ifiles(i),(/dimsizes(iyears),13/),"float")

  if(istart(i).le.1979) then
    index=tmp(ind(tmp(:,0).ge.1979),:)
  else
    index=new((/dimsizes(years),13/),"float")
    index(ind(years.ge.istart(i)),:)=tmp
  end if

 delete([/tmp,iyears/])

do s=0,dimsizes(snames)-1

sname=snames(s)
seasons=(/seasons1(s),seasons2(s)/)

ofile="MUCAPE_ERAI_Australia_methodpanel_95pc_"+inames(i)+"corr_"+sname

Cmean=new((/dimsizes(types),dimsizes(years),dimsizes(lat),dimsizes(lon)/),"integer")
Cmean!1="year"
Cmean!2="lat"
Cmean&lat=lat
Cmean!3="lon"
Cmean&lon=lon

index2=new(dimsizes(years),"float")

if(seasons(1).ge.seasons(0)) then
 mlist=ispan(seasons(0),seasons(1),1)-1
 Cmean=(/dim_sum_n_Wrap(mucape(:,:,mlist,:,:),2)/)
 index2=dim_avg_n(index(:,mlist+1),1)
 delete(mlist)
else
  list1=ispan(seasons(0)-1,11,1)
  list2=ispan(0,seasons(1)-1,1)
 do y=1,dimsizes(years)-1
  if(dimsizes(list1).eq.1) then
    Cmean(:,y,:,:)=(mucape(:,y-1,list1,:,:)+dim_sum_n_Wrap(mucape(:,y,list2,:,:),1))
  else
    Cmean(:,y,:,:)=(dim_sum_n_Wrap(mucape(:,y-1,list1,:,:),1)+dim_sum_n_Wrap(mucape(:,y,list2,:,:),1))
 end if
    index2(y)=(sum(index(y-1,list1+1))+sum(index(y,list2+1)))/(dimsizes(list1)+dimsizes(list2))
 end do
 delete([/list1,list2/])
end if 

; Plotting

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
  res@gsnLeftString="ERAI"
  res@cnLinesOn = False
  res@lbLabelBarOn        = False
  res@cnLineLabelsOn = False

opts_psl = True
     opts_psl@cnLineColor       = "Black"
     opts_psl@cnInfoLabelOn     = False
     opts_psl@cnLineLabelFontHeightF = 0.01
     opts_psl@cnLineLabelPerimOn = False
     opts_psl@gsnContourLineThicknessesScale = 4
     opts_psl@gsnDraw         = False
     opts_psl@gsnFrame        = False
     opts_psl@cnLevelSelectionMode = "ExplicitLevels"
     opts_psl@cnLevels             = (/ -100.,0.05,100./)
     opts_psl@cnLineDrawOrder="Draw"
;     opts_psl@gsnContourPosLineDashPattern = 2


wks_type = "pdf"
;wks_type@wkOrientation = "landscape"
wks_type@wkPaperWidthF = 10
wks_type@wkPaperHeightF = 6
  wks = gsn_open_wks(wks_type,ofile)             ; send graphics to PNG file
  gsn_define_colormap(wks,"MPL_RdBu")
  gsn_reverse_colormap(wks)

  plot = new(dimsizes(types),graphic)

  do p=0,dimsizes(types)-1
   res@gsnRightString=tnames(p)

   ; Do correlations

  CHcor=spcorr_n(Cmean(p,:,:,:),conform(Cmean(p,:,:,:),index2,0),0)
  CHcor!0="lat"
  CHcor&lat=lat
  CHcor!1="lon"
  CHcor&lon=lon

  plot(p) = gsn_csm_contour_map(wks,CHcor, res)     ; create plot

; And also the significance

   pval=rtest(CHcor,dimsizes(index2),0)
   pval!0="lat"
   pval&lat=lat
   pval!1="lon"
   pval&lon=lon
   contour_psl = gsn_csm_contour(wks,pval,opts_psl)

  overlay(plot(p), contour_psl)

  end do

  resP                  = True                   ; modify the panel plot
resP@gsnPaperOrientation = "portrait"
  resP@gsnMaximize      = True                   ; use full page
  resP@lbTitleFontHeightF=0.012
  resP@lbTitleString="Correlation between days MUCAPE >= 95th percentile and "+inames(i)+" in "+sname
  resP@lbTitleOn        =  True                ; turn on title
  resP@lbTitlePosition  = "Bottom"              ; title position
  resP@lbTitleDirection = "Across"             ; title direction
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@gsnPanelLabelBar    = True
  resP@lbLabelFontHeightF =.01
 resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelRowSpec = True
  gsn_panel(wks,plot,(/3,3/),resP)

delete([/Cmean,index2/])
end do
end do

