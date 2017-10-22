load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

year=2009
month=2
day=7
hour=6 ; This is in GMT

Edir="/g/data1/ub4/erai/netcdf/6hr/atmos/oper_an_pl/v01/"

tfile=addfile(Edir+"ta/ta_6hrs_ERAI_historical_an-pl_"+year+sprinti("%2.2i",month)+"01_"+year+sprinti("%2.2i",month)+days_in_month(year,month)+".nc","r")

hfile=addfile(Edir+"hur/hur_6hrs_ERAI_historical_an-pl_"+year+sprinti("%2.2i",month)+"01_"+year+sprinti("%2.2i",month)+days_in_month(year,month)+".nc","r")

time=tfile->time
lev=tfile->lev
lat=tfile->lat
lon=tfile->lon

I=ind(lon.ge.110..and.lon.le.155.)
J=ind(lat.ge.-45..and.lat.le.-10.)
L=ind(lev.eq.85000..or.lev.eq.70000.)

; COULD do time properly, but for now will do it lazily
; 124 hours -> data is six-hourly. 06Z will be best
T=4*(day-1)+(hour/6)

; Now, need our three variables

AirTemp=tfile->ta(T,L,J,I)
RH=hfile->hur(T,L,J,I)
DewPoint=dewtemp_trh(AirTemp,RH)

;
CA=0.5*(AirTemp(0,:,:)-AirTemp(1,:,:))
CB=0.3333*(AirTemp(0,:,:)-DewPoint(0,:,:))
CB=where(CB.ge.9.,9.,CB)
CB=where(CB.ge.5,5+0.5*CB,CB)
CHaines=CB+CA

CHaines!0="lat"
CHaines&lat=lat(J)
CHaines&lat@units="degrees_north"
CHaines!1="lon"
CHaines&lon=lon(I)
CHaines&lon@units="degrees_east"

; Play with plotting, since this is a test

date=year+sprinti("%2.2i",month)+sprinti("%2.2i",day)+"_"+sprinti("%2.2i",hour)+"UTC"

  wks = gsn_open_wks("pdf","chaines_"+date)             ; send graphics to PNG file
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  res                      = True               ; plot mods desired
  res@tiMainString         = "C-Haines index: "+date    ; main title
  res@cnFillOn             = True               ; turn on color fill
  res@gsnAddCyclic = False
  res@mpFillOn              = False
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

  plot = gsn_csm_contour_map(wks,CHaines, res)     ; create plot

