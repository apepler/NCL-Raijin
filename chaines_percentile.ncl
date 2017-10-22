load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

dir="/g/data/eg3/asp561/CHaines_ERAI/"

ifile=addfile(dir+"CHaines_19792016_06UTC.nc","w")
CH=ifile->CHaines

CH_percentile=new(dimsizes(CH),"integer")
copy_VarMeta(CH,CH_percentile)

a=dimsizes(CH)-1
CHpc=new((/101,a(1)+1,a(2)+1/),typeof(CH))
CHpc!0="percentile"
CHpc&percentile=ispan(0,100,1)
CHpc!1="lat"
CHpc&lat=CH&lat
CHpc!2="lon"
CHpc&lon=CH&lon

do jlat=0,a(1)
do ilon=0,a(2)
  CHsort=CH(ind(.not.ismissing(CH(:,jlat,ilon))),jlat,ilon)
  qsort(CHsort)
  do p=0,99
    CHpc(p,jlat,ilon)=(/CHsort(round(p*dimsizes(CHsort)/100,3))/)
  end do
  CHpc(100,jlat,ilon)=max(CHsort) ;

  do p=0,98
   I=ind(CH(:,jlat,ilon).ge.CHpc(p,jlat,ilon).and.CH(:,jlat,ilon).lt.CHpc(p+1,jlat,ilon))
   CH_percentile(I,jlat,ilon)=p+1
   delete(I)
  end do
  I=ind(CH(:,jlat,ilon).ge.CHpc(99,jlat,ilon)) ; making sure max included in p100
  CH_percentile(I,jlat,ilon)=100
  delete(I)
  delete(CHsort)
end do
end do

CH_percentile@standard_name="C_Haines_percentile"
CH_percentile@long_name = "C-Haines Index converted to percentile ranges, 1-100"

CHpc@standard_name="C_Haines_percentile_thresholds"
CHpc@long_name = "Value of 0-100 percentiles of C-Haines index"
  
    ifile->CH_percentile=CH_percentile
    ifile->percentile_values=CHpc
