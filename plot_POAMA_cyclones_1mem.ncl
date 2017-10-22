load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

startdate=20130101
keydate1=20130112
keydate2=20130131
;emn="e01"
; Load POAMA data

f1=addfile("/g/data2/rr8/POAMA2/hindcast/e24/e24a/"+emn+"/dac_"+startdate+".nc","r")
lat1=f1->lat
lon1=f1->lon
time1=f1->time
time1a=cd_calendar(time1,2)

x=ind(time1a.ge.keydate1.and.time1a.le.keydate2)
print(x)
prcp1=dim_sum_n_Wrap(f1->hr24_prcp(x,:,:),0)

prcp1@lat=lat1
prcp1@lon=lon1


; Plotting info 

 res                 = True
; Set some contouring resources.
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn   = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/  1., 5., 10., 25., 50., 100., 200., 300. /)
  res@cnFillColors         = (/2,11,10,6,7,8,4,5,17/)
  res@gsnDraw         = False
  res@gsnFrame        = False
; Add map resources
  res@mpFillOn              = False
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res@mpLimitMode        = "Corners"               ; Portion of map to zoom
  res@mpLeftCornerLatF   = -50.
  res@mpLeftCornerLonF   = 110.
  res@mpRightCornerLatF  = 0.
  res@mpRightCornerLonF  = 180.
  res@mpGridLineDashPattern = 2
  res@mpGridLineThicknessF = 0.5
  res@mpProjection = "CylindricalEquidistant"
; Add label bar resources
  res@lbLabelAutoStride = True
  res@lbBoxMinorExtentF        = 0.13
  res@lbLabelFontHeightF       = 0.012
  res@gsnLeftStringFontHeightF = 0.01
  res@gsnRightStringFontHeightF = 0.008
  res@gsnMaximize     = False    ; Maximize plot in frame
  res@gsnLeftString="24 hour rainfall (mm)"
  res@gsnRightString=""

; Plotting

wks = gsn_open_wks("pdf","POAMA_"+emn+"_prcyclones_"+startdate+"_"+keydate1+"_"+keydate2)
gsn_define_colormap(wks,"prcp_1")                ; choose colormap

res@tiMainString="POAMA rain: "+keydate1+" to "+keydate2 
plot = gsn_csm_contour_map(wks, prcp1, res)


; For start point - 
pmres=True
pmres@gsMarkerColor = "black"
pmres@gsMarkerIndex = 5
pmres@gsMarkerThicknessF = 2
pmres@gsMarkerSizeF = 0.01

; For line

  lres                  = True           ; polyline mods desired
  lres@gsLineThicknessF = 4.0            ; line thickness   
  lres@gsLineColor      = "black"         ; color of lines

; Get all cyclones
f2="/short/eg3/asp561/cts.dir/gcyc_out/POAMA/proj240_rad5cv0.15/e24a/"+emn+"/tracks_20130101_4month.dat"
data=readAsciiTable(f2,14,"float",0)
I=ind(data(:,6).ge.135.and.data(:,6).le.155.and.data(:, 7).ge.-40.and.data(:,7).le.-10.and.data(:,2).ge.mod(keydate1,1000000).and.data(:,2).le.mod(keydate2,1000000))
ids=get_unique_values(data(I,0))

dum1=new(dimsizes(ids),graphic)
dum2=dum1

do ip=0,dimsizes(ids)-1
  I2=ind(data(:,0).eq.ids(ip))
  dum1(ip)=gsn_add_polyline(wks,plot,data(I2,6),data(I2,7),lres)
  dum2(ip)=gsn_add_polymarker(wks,plot,data(I2(0),6),data(I2(0),7),pmres)
  delete(I2)
end do


draw(plot)

