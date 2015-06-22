#!/bin/csh

set date=`date +%y%m%d` ; # Set date for dynamic run

# This creates the files like
# MODEL.OI2.ice.mnly.yyyymm-YYYYMM.unf.nc
# MODEL.OI2.sst.mnly.yyyymm-YYYYMM.unf.nc
# ==============================NCL====================================

### setenv NCARG_ROOT /usr/local/bin

cat >! main.ncl << "END_MAIN_NCL"
;; Not needed from 6.2 onwarded
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

procedure linInterpFlip( X, nPass, NX, NY, endPt)
; source: coast_land_NearNbor_bilin.ncl
;         modified to take 0.5 to 359.5 as input
;         use loop to minimize memory
local np, y, nMSG, nt, nx, ny, dimX, ntim
begin
  nx = NX
  ny = NY
  
  dimX = dimsizes(X)
  ntim = dimX(0)
 do nt=0,ntim-1
   x = lonFlip( X(nt,:,:) )
  do np=1,nPass
     x      = linmsg(x, (/endPt, nx/))    ; interpolate in longitude
     if (ny.gt.0) then                    ; interpolate in y
         y  = x(lon|:,lat|:)       ; reorder for interpolation
         y  = linmsg(y, (/endPt,ny/))     

         x  = y(lat|:,lon|:)       ; reorder to original
     end if

     nx     = nx + 10                     ; expand the number of pts
     if (np.gt.4) then
         ny = ny + 3
     end if
  end do
  X(nt,:,:) = (/ lonFlip(x) /)
 end do

  nMSG = num(ismissing(X))
  print("****************************")
  print("linInterpFlip: nMSG="+nMSG)
  print("****************************")
  if (nMSG.ne.0) then
      printVarSummary(X)
      do nt=0,ntim-1
         nMSG = num(ismissing(X(nt,:,:)))
         if (nMSG.ne.0) then
             print("nt="+nt+"  nMSG="+nMSG )
         end if
      end do
      exit
  end if
end

external SSTICE  "./sstice.so" 
external JIMH "./consistent.so"
external NEAR_NEIGHBOR "./coast_land_NearNbor.so"

begin                       ; MAIN NCL DRIVER
  netCDF      = True 
  DEBUG       = False 
  PRNT_DEBUG  = True 
  PLOT_DEBUG  = False
  PLOT_TYPE   = "ps"
                            ; main input directory
 ;diri   = "/ptmp/shea/SSTICE/NEW_SST/"
 ;diro   = "/ptmp/shea/SSTICE/"
 ;dirm   = "/cgd/cas/shea/JHURRELL/"             ; mask ... not used
 ;diri   = "/project/cas/shea/SSTICE/SST_NEW/"
 ;diro   = "/project/cas/shea/SSTICE/"
 ;dirm   = "/project/cas/shea/SSTICE/"           ; mask ... not used
 ;diri   = "/work/durack1/Shared/150219_AMIPForcingData/SST_NEW/$date/"
  diri   = "/work/durack1/Shared/150219_AMIPForcingData/SST_NEW3/150610/"
  diro   = "/work/durack1/Shared/150219_AMIPForcingData/SST_NEW3/"
  dirm   = "/work/durack1/Shared/150219_AMIPForcingData/SST_NEW3/"
  film   = "lstags.onedeg.dat"  ; sst mask [not used]

  fils   = systemfunc ("cd "+diri+"; ls oiv2mon.* ")  ; the ftp's file names
  nfils  = dimsizes(fils)

;;filoi = "MODEL.OI2.ice.mnly.200501-200903.unf.nc"
;;filos = "MODEL.OI2.sst.mnly.200501-200903.unf.nc"
  sfxc  = stringtochar( get_file_suffix(fils(0),0) )
  YYYYMM_START = chartostring(sfxc(1:))
  sfxc  = stringtochar( get_file_suffix(fils(nfils-1),0) )
  YYYYMM_LAST  = chartostring(sfxc(1:))
  filoi = "MODEL.OI2.ice.mnly."+YYYYMM_START+"-"+YYYYMM_LAST+".unf.nc"
  filos = "MODEL.OI2.sst.mnly."+YYYYMM_START+"-"+YYYYMM_LAST+".unf.nc"
  print("filoi="+filoi)
  print("filos="+filos)

  ntim   = nfils

  nfStrt = 0                 ; nfStrt = nfils-1
  nfLast = nfils-1

  if (PLOT_DEBUG) then
     wks  = gsn_open_wks(PLOT_TYPE,"test")    ; open workstation (plot destination)
     gsn_define_colormap(wks,"BlGrYeOrReVi200") ; choose colormap  
     gray = NhlNewColor(wks,0.8,0.8,0.8)   ; add gray to colormap

     res                      = True
     res@cnFillOn             = True        ; turn on color
     res@gsnSpreadColors      = True        ; use full range of colormap
     res@gsnSpreadColorStart  = 2           ; start at color 2
     res@gsnSpreadColorEnd    = -3          ; don't use added gray
     res@cnLinesOn            = False       ; no contour lines
     res@cnLineLabelsOn       = False       ; no contour lines
     res@cnInfoLabelOn        = False       ; turn off cn info label
    ;res@cnFillDrawOrder      = "PreDraw"   ; draw contours before continents
     res@gsnMaximize          = True        ; maximize plot
    ;res@mpFillOn             = False
     res@mpCenterLonF         = 200.

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF       =   0.               ; set min contour level
     res@cnMaxLevelValF       =  32.               ; set max contour level
     res@cnLevelSpacingF      =   2.               ; set contour spacing

     RES = res
     RES@cnMinLevelValF       =   0.               ; set min contour level
     RES@cnMaxLevelValF       = 100.               ; set max contour level
     RES@cnLevelSpacingF      =   5.               ; set contour spacing

     nfStrt = nfils-1
     ntim   = nfLast-nfStrt+1
  end if
 ;print(fils(nfStrt:nfLast))

  year   = new (ntim, "integer") 
  month  = new (ntim, "integer") 

  print("ntim="+ntim+" nfStrt="+nfStrt+"  nfLast="+nfLast)

  nt = -1
  do nf=nfStrt,nfLast
     tmp_c     = stringtochar(fils(nf))
     nt = nt + 1
     year(nt)  = stringtointeger((/tmp_c( 8:11)/))
     month(nt) = stringtointeger((/tmp_c(12:13)/))
     delete(tmp_c)
  end do

print(year+"   "+month)
  ntStrt = 0
  ntLast = ntim

  day    = new (ntim, "integer")
  day    = 16                         ; default
  hour   = new (ntim, "integer")
  hour   = 12                         ; default
  minute = new (ntim, "integer")
  minute = 0
  sec    = new (ntim, "double")
  sec    = 0.0d0

  datesec= new (ntim, "integer")
  delete(datesec@_FillValue)
  datesec@units = "current seconds of current date"
  datesec!0     = "time" 
  datesec       = 43200               ; default

  date_frac= new (ntim, "double")
  delete(date_frac@_FillValue)
  date_frac@units = "yyyymmdd.fraction_of_day"
  date_frac!0     = "time" 
  date_frac       = 0.5d0

  i      = ind(month.eq.2)            ; ignore leap year
  if (.not.any(ismissing(i))) then
      day(i) = 15
  end if
  delete(i)

  i      = ind(month.eq.2 .or. month.eq.4 .or. month.eq.6 .or. \
               month.eq.9 .or. month.eq.11     ) 
  if (.not.any(ismissing(i))) then
      hour(i)      = 0  
      datesec(i)   = 0
      date_frac(i) = 0.d0             
  end if
  delete(i)

  units  = "days since 1800-01-01 00:00:00"
  time   = ut_inv_calendar(year,month,day,hour,minute,sec, units, 0)
  time!0 = "time"
  delete(time@_FillValue)
  time@information = "middle of month"
  printVarSummary(time)

  date       = ut_calendar(time, -2)     
  date!0     = "time"
  date@units = "yyyymmdd"

  date_frac  = (/ date + date_frac /)

  print("nfStrt="+nfStrt)
  print("nfLast="+nfLast)
  print(fils(nfStrt:nfLast) \
            +"  "+year+"  "+month+"  "+day+"  "+hour \
            +"  "+date+"  "+datesec+"  "+date_frac)

  nlat = 180
  mlon = 360

  lat  = latGlobeFo(nlat, "lat", "latitude", "degrees_north")
  lon  = lonGlobeFo(mlon, "lon", "longitude", "degrees_east")

  smsg = 1.e20
  sice = -1.8
                         ; create various arrays to be used
  sst           = new ( (/ntim,nlat,mlon/) , float, smsg)
  sst!0         = "time"
  sst!1         = "lat"
  sst!2         = "lon"
  sst&time      = time
  sst&lat       = lat
  sst&lon       = lon
                                          ; not really needed
                                          ; done in fortran
  lsMask_OI     = new ( (/nlat,mlon/) , float) ; 0=land , 1=ocean
  lsMask_OI@long_name = "land-sea mask"
  delete(lsMask_OI@_FillValue)
  lsMask_OI!0   = "lat"
  lsMask_OI!1   = "lon"
  lsMask_OI&lat =  lat
  lsMask_OI&lon =  lon

  sst@longName  = "Sea Surface Temperature"
  sst@units     = "degC"
  sst@info      = "sst-ice consistency enforced" 

  ice = sst          ; copy meta data
  ice@longName  = "sea-ice concentration"
  ice@units     = "%"

  SST = sst          ; SST will be the data written to grid
  ICE = ice          ; ICE will be the data written to grid
  
  SST@ice       = sice    

;********************************************************
; Loop over NCEP-OI files
; For memory reasons do manipulations on a per/file basis
;********************************************************

  YYYYMM    = date/100

;********************************************************
; Empirical relationship
;********************************************************
                                        ; Emp=> Empirical
  iceEmp = fspan(0,0.90,91)
  sstEmp = 9.328*(0.729-iceEmp^3) - 1.8  ; sst = f(ice) 
  sstEmp@long_name = "SST MAX: EMPIRICAL"
  sstEmp@units     = "C"
  sstCrit = 5.0     ; 9.328*0.729 - 1.8

  print(fils)

  do nf=nfStrt,nfLast    ; loop over all files
     FNAME = diri+fils(nf)                          ; 19 NOV 2007
     sfx = get_file_suffix(FNAME,0) 

     print("nf="+nf+"  sfx="+sfx+"   FNAME="+FNAME)
     if (sfx.eq.".gz") then
         print("SUFFIX .gz encountered: Better to manually gzip -d")
	 exit
        ;system("gzip -d "+FNAME)
     end if

     filc = stringtochar(fils(nf) )
     yyyy = stringtointeger((/filc( 8:11)/))
     mm   = stringtointeger((/filc(12:13)/))
     mm_s = (/filc(12:13)/) ; string

     yyyymm = yyyy*100 + mm
     nt     = ind(yyyymm.eq.YYYYMM)
     print("nf="+nf+"  nt="+nt+"  yyyymm="+yyyymm)

     if (ismissing(nt(0))) then
         print("nt is missing:  nf="+nf+"  yyyymm="+yyyymm)
         exit
     end if
print("DEBUG A:+++++++++++++++++++++++")
print("fils(nf)="+fils(nf))
                                   ; read NCEP data 
                                   ; lsMask_OI added 3/27/2006
     SSTICE::sstice  (dirm, diri, film,  fils(nf) \
                     ,yyyy, mm, nlat, mlon, ice(nt,:,:) \
                     ,sst(nt,:,:),lsMask_OI, sst@_FillValue  )
print("DEBUG B:+++++++++++++++++++++++")

;********************************************************
; gross error checks  [basically, protect from round off error]
;********************************************************
     ice(nt,:,:) = ice(nt,:,:) > 0.0     
     ice(nt,:,:) = ice(nt,:,:) < 100.0     ; %

     dim2d = dimsizes(sst(nt,:,:))

;********************************************************
; Jim H's consistency fortran code [minor additions by DJS]
; This code assume ice@units = %
;********************************************************
     JIMH::ssticejh (nlat, mlon, ice(nt,:,:), sst(nt,:,:), sst@_FillValue )
    ;if (DEBUG) then
         print("date="+date(nt)+"  min(ice)="+min(ice(nt,:,:)) \
                               +"  max(ice)="+max(ice(nt,:,:)) \
                               +"  min(sst)="+min(sst(nt,:,:)) \
                               +"  max(sst)="+max(sst(nt,:,:)) )
    ;end if

;========================================================
; For *clarity*, explicitly do each  of the
; following data operations separately. 
;========================================================

;********************************************************
; Set all sst < -1.8 to -1.8   [JimH did this in fortran]
; Safety check to make sure this happened
;********************************************************
     sst(nt,:,:)   = sst(nt,:,:) > sice

;********************************************************
; where ice > 90 [%], set corresponding sst to -1.8
; f90:  where(ice > 90.) sst = -1.8
; Safety check to make sure this happened
;********************************************************
;    sst(nt,:,:) = where (ice(nt,:,:).ge.90, sice, sst(nt,:,:))
;********************************************************

     sst1d = ndtooned( sst(nt,:,:) )
     ice1d = ndtooned( ice(nt,:,:) )
     i     = ind(ice1d.ge.90)
     if (.not.any(ismissing(i))) then
         sst1d(i) = sice
     end if
     sst(nt,:,:)  = onedtond(sst1d, dim2d)
     delete(i)
     delete(sst1d)
     delete(ice1d)

;********************************************************
; where ice < 15 [%], reset ice to 0.0  
; Safety check to make sure this happened
;********************************************************
;    ice(nt,:,:) = where (ice(nt,:,:).lt.15,  0.0 , ice(nt,:,:))
;********************************************************

     sst1d = ndtooned( sst(nt,:,:) )
     ice1d = ndtooned( ice(nt,:,:) )
     i     = ind(ice1d.lt.15)
     if (.not.any(ismissing(i))) then
         ice1d(i) = 0.0
     end if
     ice(nt,:,:)  = onedtond(ice1d, dim2d)
     delete(i)
     delete(sst1d)
     delete(ice1d)

;********************************************************
; Any place where sst>sstCrit set the ice=0
; Safety check to make sure this happened
;********************************************************
;    ice(nt,:,:) = where (ice(nt,:,:).gt.0 .and. sst(nt,:,:).ge.sstCrit \
;                        , 0.0 , ice(nt,:,:))
;********************************************************

     ice1d = ndtooned( ice(nt,:,:) )
     sst1d = ndtooned( sst(nt,:,:) )
     i     = ind(ice1d.gt.0 .and. sst1d.ge.sstCrit)
     if (.not.any(ismissing(i))) then
         ice1d(i) = 0.0
     end if
     ice(nt,:,:) = onedtond(ice1d, dim2d)
     delete(i)
     delete(sst1d)
     delete(ice1d)

;********************************************************
; Under the assumption that SST is more reliable than 
; sea-ice concentration, where (15<=ice<=90)  and the
; sst exceed the empirically determined max .. adjust the sea-ice.
; Make sure tenths are used. Then put into ice_pc
;********************************************************
                                              ; empirical formula expects tenths
     ice1d = ndtooned( ice(nt,:,:)*0.01 )     ; tenths  
     sst1d = ndtooned( sst(nt,:,:) )  

     n1590 = num(ice1d.ge.0.15 .and. ice1d.lt.0.90\
                               .and. .not.ismissing(sst1d) )
     i     = ind(ice1d.ge.0.15 .and. ice1d.lt.0.90\
                               .and. .not.ismissing(sst1d) )

     if (.not.any(ismissing(i))) then
                                                      ; alter only following
         ni = dimsizes(i)
         do n=0,ni-1
            sstmx = 9.328*(0.729-ice1d(i(n))^3) - 1.8
            if (sst1d(i(n)).gt.sstmx) then
                ice1d(i(n))  = exp( log(0.729-((sstmx+1.8)/9.328))/3.0 ) ; tenths
            end if
         end do

    
         sstmx1d    = sst1d                           ; exact copy
         sstmx1d(i) = 9.328*(0.729-ice1d(i)^3) - 1.8  ; empirical max sst
         k          = ind(ice1d.ge.0.15 .and. ice1d.lt.0.90 \
                                        .and. sst1d.gt.sstmx1d)
         dimk       = dimsizes(k)
         if (.not.any(ismissing(k))) then
                                   ; calculate the empirical ice max
             ice1d(k)  = exp( log(0.729-((sstmx1d(k)+1.8)/9.328))/3.0 ) ; tenths
         end if

         delete(k)
         delete(sst1d)
         delete(sstmx1d)
     end if

     ice(nt,:,:) = onedtond(ice1d*100, dim2d)      ; return to %
     ice(nt,:,:) = ice(nt,:,:) > 0.0

     delete(i)
     delete(ice1d)

;********************************************************
; Any place where ice < 15% set to 0.0
;********************************************************
;    ice(nt,:,:) = where (ice(nt,:,:).lt.15,  0.0 , ice(nt,:,:))
;********************************************************

     ice1d = ndtooned( ice(nt,:,:) )
     i     = ind(ice1d.lt.15)
     if (.not.any(ismissing(i))) then
         ice1d(i) = 0.0                    ; %
     end if
     ice(nt,:,:) = onedtond(ice1d, dim2d)
     delete(i)
     delete(ice1d)

  end do

  N18   = num(sst.lt.sice)                       
  N15   = num(ice.gt.0.   .and. ice.lt.15)                       
  N15C  = num(ice.gt.0.   .and. ice.lt.15 .and. sst.gt.sstCrit)                       
  N1590 = num(ice.ge.15 .and. ice.lt.90 .and. sst.gt.sstEmp(1) )
  N90   = num(ice.ge.90 .and. sst.gt.sice)
  print("AFTER:  N15="+N15+"   N15C="+N15C+"  N18="+N18+"  N90="+N90+"  N1590="+N1590)

  sst@info = "sst-ice consistency enforced"
  ice@info = "sst-ice consistency enforced"

  if (PLOT_DEBUG) then
      nt   = ntim-1
      res@gsnCenterString = "SST Grid: "+date(nt)
      plot = gsn_csm_contour_map_ce(wks,sst(nt,:,:), res)  
      RES@gsnCenterString = "ICE Grid: "+date(nt)
      plot = gsn_csm_contour_map_ce(wks,ice(nt,:,:), RES)  

      res_OI = True
      res_OI@gsnSpreadColors      =  True
      res_OI@gsnSpreadColorEnd    = -3          ; don't use added gray
      res_OI@cnFillOn             =  True
      res_OI@cnFillMode           = "CellFill"
     ;res_OI@cnMissingValFillColor= "yellow"
      res_OI@mpFillOn             =  False
      res_OI@mpFillDrawOrder      =  "PostDraw"     
      res_OI@gsnCenterString      = "SST Missing"
      res_OI@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
      res_OI@cnMinLevelValF       =   0                ; set min contour level
      res_OI@cnMaxLevelValF       =   2                ; one less than max
      res_OI@cnLevelSpacingF      =   1                ; set contour spacing
     ;res_OI@lbLabelStrings       = ispan(1,3,1)   
      if (isatt(res,"mpCenterLonF")) then
          res_OI@mpCenterLonF     = res@mpCenterLonF
      end if
      plot = gsn_csm_contour_map_ce(wks,lsMask_OI, res_OI)  
     ;print(lat+"  "+lsMask_OI(:,{0.5})+"    "+lsMask_OI(:,{179.5}) )
  end if
                                      ; Antarctic mask  
                                      ; NCEP: 0=land , ocean=1
  latMask1d  = ndtooned(conform(lsMask_OI, lat, 0))
  lsMask1d   = ndtooned(lsMask_OI)
  iMask      = ind(latMask1d.le.-60 .and. lsMask1d.eq.0)
                                      ; force certain values 
  ice(:,{-45:35} ,:) =   0.0

 ;wcStrt = systemfunc("date")
  do nt=0,ntim-1
     x1d = ndtooned(sst(nt,:,:))
     x1d(iMask)   = sice              ; Antarctic land
     sst(nt,:,:)  = onedtond(x1d, (/nlat,mlon/) )

     x1d = ndtooned(ice(nt,:,:))
     x1d(iMask)   = 100.              ; Antarctic land
     ice(nt,:,:)  = onedtond(x1d, (/nlat,mlon/) )
  end do
 ;wallClockElapseTime(wcStrt, "MASK", 0)

  if (PLOT_DEBUG) then
      nt   = ntim-1
      res@gsnCenterString = "set -1.8"
      plot = gsn_csm_contour_map_ce(wks,sst(nt,:,:), res)  
      RES@gsnCenterString = "set 100"
      plot = gsn_csm_contour_map_ce(wks,ice(nt,:,:), RES)  
  end if

                               ; will be inserted over Siberia
  sst_zon= dim_avg_Wrap(sst(:,:,{0:180}))
  sst_zon@long_name = "Zonal Mean SST"
  ice_zon= dim_avg_Wrap(ice(:,:,{0:180}))
  ice_zon@long_name = "Zonal Mean SEAICE"
                               ; completly bogus [again]
                               ; set China/Siberia to the  zonal average
                               ; makes code 'converge' faster
  if (PRNT_DEBUG) then
      nMsgS  = num(ismissing(sst))         ; total number of _FillValue
      print("SST: start: nMsgS ="+nMsgS)
      nMsgI   = num(ismissing(ice))        ; total number of _FillValue
      print("ICE: start: nMsgI ="+nMsgI)
  end if
  wcStrt = systemfunc("date")

;***************************************************
; Begin Interpolation over land
;***************************************************
  nDx = 2                      ; longitude only 
  sst = linmsg (sst, (/0,nDx/)); linearly interpolate over small distances
  ice = linmsg (ice, (/0,nDx/)) 

  if (PRNT_DEBUG) then
      nMsgS  = num(ismissing(sst))         ; total number of _FillValue
      print("SST: linmsg: nDX="+nDx+" : nMsgS ="+nMsgS)
      nMsgI  = num(ismissing(ice))         ; total number of _FillValue
      print("ICE: linmsg: nDX="+nDx+" : nMsgI ="+nMsgI)
  end if

  if (PLOT_DEBUG) then
      nt   = ntim-1
      res@gsnCenterString = "linmsg="+nDx
      plot = gsn_csm_contour_map_ce(wks,sst(nt,:,:), res)  
      plot = gsn_csm_contour_map_ce(wks,ice(nt,:,:), RES)  
  end if
                               ; assumes -180 to +180
  nPtx  = 3                    ; small inland nearest neighbor   
  nPty  = 1  
                               ; coastnn works better when 179.5W -> 179.5E
 ;NEAR_NEIGHBOR::coastnn(sst,mlon,nlat,ntim,sst@_FillValue,nPtx,nPty)
 ;NEAR_NEIGHBOR::coastnn(ice,mlon,nlat,ntim,ice@_FillValue,nPtx,nPty)
                               
  do nt=0,ntim-1               ; use loop to minimize memory
     tmp = lonFlip(sst(nt:nt,:,:))         ; make NCEP 179.5W -> 179.5E
     NEAR_NEIGHBOR::coastnn(tmp,mlon,nlat, 1  ,tmp@_FillValue,nPtx,nPty)
     sst(nt:nt,:,:) = (/ lonFlip(tmp) /)

     tmp = lonFlip(ice(nt:nt,:,:))
     NEAR_NEIGHBOR::coastnn(tmp,mlon,nlat, 1  ,tmp@_FillValue,nPtx,nPty)
     ice(nt:nt,:,:) = (/ lonFlip(tmp) /)
  end do

  if (PRNT_DEBUG) then
      nMsgS  = num(ismissing(sst))         ; total number of _FillValue
      print("SST: NEAR NEIGHBOR: nPtx="+nPtx+"  nPty="+nPty+"  nMsgS ="+nMsgS)
      nMsgI  = num(ismissing(ice))         ; total number of _FillValue
      print("ICE: NEAR NEIGHBOR: nPtx="+nPtx+"  nPty="+nPty+"  nMsgI ="+nMsgI)
  end if

  if (PLOT_DEBUG) then
      nt   = ntim-1
      res@gsnCenterString = "NearNeighbor: nPtx="+nPtx+"  nPty="+nPty
      plot = gsn_csm_contour_map_ce(wks,sst(nt,:,:), res)  
      plot = gsn_csm_contour_map_ce(wks,ice(nt,:,:), RES)  
  end if
                                                    ; hasten fill area
  sst(:,{25:70},{100}) = (/ sst_zon(:,{25:70}) /)   ; arbitrary
  sst(:,{50:60},{ 55}) = (/ sst_zon(:,{50:60}) /)   ; arbitrary
  ice(:,{50:60},{ 55}) = (/ ice_zon(:,{50:60}) /)   ; arbitrary
                     ; large scale iterative bilinear interpolation
  nx     = 10        ; each is one deg  
  ny     = 3
  endPt  = 0
  nPass  = 8
  linInterpFlip( sst, nPass, nx, ny, endPt)
  if (PRNT_DEBUG) then
      nMsgS  = num(ismissing(sst))         ; total number of _FillValue
      print("SST: nPass="+nPass+" : nMsgS ="+nMsgS)
      print(" ")
  end if
  nPass  = 8
  linInterpFlip( ice, nPass, nx, ny, endPt)

  wallClockElapseTime(wcStrt, "NearNeighbor-LinInterp", 0)

  if (PRNT_DEBUG) then
     ;nMsgS  = num(ismissing(sst))         ; total number of _FillValue
     ;print("SST: nPass="+nPass+" : nMsgS ="+nMsgS)
     ;print(" ")
      nMsgI  = num(ismissing(ice))         ; total number of _FillValue
      print("ICE: nPass="+nPass+" : nMsgI ="+nMsgI)
      print(" ")
  end if

  if (PLOT_DEBUG) then
      nt   = ntim-1
      res@gsnCenterString = "linInterp: nPass="+nPass
      plot = gsn_csm_contour_map_ce(wks,sst(nt,:,:), res)  
      plot = gsn_csm_contour_map_ce(wks,ice(nt,:,:), RES)  
  end if

  printMinMax(sst , True)     
  print("nMsg(sst)="+num(ismissing(sst)))
  printMinMax(ice , True)         
  print("nMsg(ice)="+num(ismissing(ice)))

  sst@info = "sst-ice consistency enforced"
  ice@info = "sst-ice consistency enforced"

;********************************************************
; Force consistency: make sure land interp does not screw up
;********************************************************
  do nt=0,ntim-1

;********************************************************
; Set all sst < -1.8 to -1.8   [JimH did this in fortran]
; Safety check to make sure this happened
;********************************************************
     sst(nt,:,:)   = sst(nt,:,:) > sice

;********************************************************
; where ice > 90 [%], set corresponding sst to -1.8
; f90:  where(ice > 90.) sst = -1.8
; Safety check to make sure this happened
;********************************************************
;    sst(nt,:,:) = where (ice(nt,:,:).ge.90, sice, sst(nt,:,:))
;********************************************************

     sst1d = ndtooned( sst(nt,:,:) )
     ice1d = ndtooned( ice(nt,:,:) )
     i     = ind(ice1d.ge.90)
     if (.not.any(ismissing(i))) then
         sst1d(i) = sice
     end if
     sst(nt,:,:)  = onedtond(sst1d, dim2d)
     delete(i)
     delete(sst1d)
     delete(ice1d)

;********************************************************
; where ice < 15 [%], reset ice to 0.0  
; Safety check to make sure this happened
;********************************************************
;    ice(nt,:,:) = where (ice(nt,:,:).lt.15,  0.0 , ice(nt,:,:))
;********************************************************

     sst1d = ndtooned( sst(nt,:,:) )
     ice1d = ndtooned( ice(nt,:,:) )
     i     = ind(ice1d.lt.15)
     if (.not.any(ismissing(i))) then
         ice1d(i) = 0.0
     end if
     ice(nt,:,:)  = onedtond(ice1d, dim2d)
     delete(i)
     delete(sst1d)
     delete(ice1d)

;********************************************************
; Any place where sst>sstCrit set the ice=0
; Safety check to make sure this happened
;********************************************************
;    ice(nt,:,:) = where (ice(nt,:,:).gt.0 .and. sst(nt,:,:).ge.sstCrit \
;                        , 0.0 , ice(nt,:,:))
;********************************************************

     ice1d = ndtooned( ice(nt,:,:) )
     sst1d = ndtooned( sst(nt,:,:) )
     i     = ind(ice1d.gt.0 .and. sst1d.ge.sstCrit)
     if (.not.any(ismissing(i))) then
         ice1d(i) = 0.0
     end if
     ice(nt,:,:) = onedtond(ice1d, dim2d)
     delete(i)
     delete(sst1d)
     delete(ice1d)

;********************************************************
; Under the assumption that SST is more reliable than 
; sea-ice concentration, where (15<=ice<=90)  and the
; sst exceed the empirically determined max .. adjust the sea-ice.
; Make sure tenths are used. Then put into ice_pc
;********************************************************
                                              ; empirical formula expects tenths
     ice1d = ndtooned( ice(nt,:,:)*0.01 )  
     sst1d = ndtooned( sst(nt,:,:) )  

     n1590 = num(ice1d.ge.0.15 .and. ice1d.lt.0.90\
                               .and. .not.ismissing(sst1d) )
     i     = ind(ice1d.ge.0.15 .and. ice1d.lt.0.90\
                               .and. .not.ismissing(sst1d) )

     if (.not.any(ismissing(i))) then
                                                      ; alter only following
         ni = dimsizes(i)
         do n=0,ni-1
            sstmx = 9.328*(0.729-ice1d(i(n))^3) - 1.8
            if (sst1d(i(n)).gt.sstmx) then
                ice1d(i(n))  = exp( log(0.729-((sstmx+1.8)/9.328))/3.0 ) ; tenths
            end if
         end do

    
         sstmx1d    = sst1d                           ; exact copy
         sstmx1d(i) = 9.328*(0.729-ice1d(i)^3) - 1.8  ; empirical max sst
         k          = ind(ice1d.ge.0.15 .and. ice1d.lt.0.90 \
                                        .and. sst1d.gt.sstmx1d)
         dimk       = dimsizes(k)
         if (.not.any(ismissing(k))) then
                                   ; calculate the empirical ice max
             ice1d(k)  = exp( log(0.729-((sstmx1d(k)+1.8)/9.328))/3.0 ) ; tenths
         end if

         delete(k)
         delete(sst1d)
         delete(sstmx1d)
     end if

     ice(nt,:,:) = onedtond(ice1d*100, dim2d)      ; return to %
     ice(nt,:,:) = ice(nt,:,:) > 0.0

     delete(i)
     delete(ice1d)

;********************************************************
; Any place where ice < 15% set to 0.0
;********************************************************
;    ice(nt,:,:) = where (ice(nt,:,:).lt.15,  0.0 , ice(nt,:,:))
;********************************************************

     ice1d = ndtooned( ice(nt,:,:) )
     i     = ind(ice1d.lt.15)
     if (.not.any(ismissing(i))) then
         ice1d(i) = 0.0                    ; %
     end if
     ice(nt,:,:) = onedtond(ice1d, dim2d)
     delete(i)
     delete(ice1d)

  end do


  N18   = num(sst.lt.sice)                       
  N15   = num(ice.gt.0.   .and. ice.lt.15)                       
  N15C  = num(ice.gt.0.   .and. ice.lt.15 .and. sst.gt.sstCrit)                       
  N1590 = num(ice.ge.15 .and. ice.lt.90 .and. sst.gt.sstEmp(1) )
  N90   = num(ice.ge.90 .and. sst.gt.sice)
  print("AFTER:  N15="+N15+"   N15C="+N15C+"  N18="+N18+"  N90="+N90+"  N1590="+N1590)

;********************************************************
; Write netCDF
;********************************************************
if (netCDF) then

  nline         = inttochar(10)

  fAtt          = True
  fAtt@title    = "Hurrell Consistent SST with no missing values over land"
  fAtt@OI_clim  = "1971-2000"
  fAtt@NCEP_OI_clim  = "1971-2000"

  fAtt@story    = nline + \
"    NCL version of JimH sst-ice.consistency.f                        " + nline + \
"                                                                     " + nline + \
"[a] all sst < -1.8 to -1.8  [sea-ice value]                          " + nline + \
"[b] where(ice > 0.9) sst = -1.8                                      " + nline + \ 
"[c] adjust the sst-ice concentration data via Jim Hack formula       " + nline + \
"    sstm = 9.328 * (0.729-si**3) - 1.8  [where si is frac sea-ice]   " + nline + \
"[d] Use ad-hoc  method to interpolate values over land               " + nline + \
"    values near coasts set to nearest neighbor                       " + nline + \
"    function cssgrid used to interpolat [not nice]                   " 

  fAtt@creator  = "Dennis Shea, CGD"
  fAtt@creation_date = systemfunc( "date" )

  ncfile = diro+filos
  print (ncfile)
  system ("/bin/rm -f " + ncfile)         ; remove an pre-file              

  ncdf   = addfile(ncfile,"c")            ; "c"reate the netCDF file

  setfileoption(ncdf,"DefineMode",True)   ; EFFICIENCY

  fileattdef( ncdf, fAtt )

  dimNames = (/ "time", "lon", "lat"  /)
  dimSizes = (/   -1  ,  mlon, nlat   /)
  dimUnlim = (/  True , False, False  /)
  filedimdef( ncdf, dimNames, dimSizes, dimUnlim )

  filevardef   ( ncdf, "time", typeof(time), getvardims(time) )
  filevarattdef( ncdf, "time", time )
                                         ; Define 1D variables.
  filevardef   ( ncdf, "lon", typeof(lon), getvardims(lon) )
  filevarattdef( ncdf, "lon", lon )

  filevardef   ( ncdf, "lat", typeof(lat), getvardims(lat) )
  filevarattdef( ncdf, "lat", lat )

  filevardef   ( ncdf, "date", typeof(date), getvardims(date))
  filevarattdef( ncdf, "date", date )

  filevardef   ( ncdf, "datesec", typeof(datesec), getvardims(datesec) )
  filevarattdef( ncdf, "datesec", datesec )

  filevardef   ( ncdf, "date_frac", typeof(date_frac), getvardims(date_frac) )
  filevarattdef( ncdf, "date_frac", date_frac )

  filevardef   ( ncdf, "SST", typeof(sst), getvardims(sst))
  filevarattdef( ncdf, "SST", sst )

  setfileoption(ncdf,"DefineMode",False)       ; (not really necessary)
                                            
  ncdf->time     = (/time  /)  
  ncdf->lat      = (/lat   /)   
  ncdf->lon      = (/lon   /)   
  ncdf->date     = (/date  /)   
  ncdf->datesec  = (/datesec    /)   
  ncdf->date_frac= (/date_frac  /)   
  ncdf->SST      = (/sst   /) 


; ------------- sea-ice netCDF

  gAtt          = True
  gAtt@title    = "Hurrell Consistent ICE with no missing values over land"

  gAtt@story    = nline + \
"    NCL version of JimH sst-ice.consistency.f                        " + nline + \
"                                                                     " + nline + \
"[a] all sst < -1.8 to -1.8  [sea-ice value]                          " + nline + \
"[b] where(ice > 0.9) sst = -1.8                                      " + nline + \ 
"[c] adjust the sst-ice concentration data via Jim Hack formula       " + nline + \
"    sstm = 9.328 * (0.729-si**3) - 1.8  [where si is frac sea-ice]   " + nline + \
"[d] minor sea-ice adjustments                                        " 

  gAtt@creators = "Dennis Shea, CGD"
  gAtt@creation_date = systemfunc( "date" )

  NCFILE = diro+filoi
  print (NCFILE)
  system ("/bin/rm -f " + NCFILE)  ; remove an pre-file              

  NCDF   = addfile(NCFILE,"c")       ; "c"reate the netCDF file

  setfileoption(NCDF,"DefineMode",True)   ; EFFICIENCY

  fileattdef( NCDF, gAtt )

  dimNames = (/ "time", "lon", "lat"  /)
  dimSizes = (/   -1  ,  mlon, nlat   /)
  dimUnlim = (/  True , False, False  /)
  filedimdef( NCDF, dimNames, dimSizes, dimUnlim )

  filevardef   ( NCDF, "time", typeof(time), getvardims(time) )
  filevarattdef( NCDF, "time", time )
                                         ; Define 1D variables.
  filevardef   ( NCDF, "lon", typeof(lon), getvardims(lon) )
  filevarattdef( NCDF, "lon", lon )

  filevardef   ( NCDF, "lat", typeof(lat), getvardims(lat) )
  filevarattdef( NCDF, "lat", lat )

  filevardef   ( NCDF, "date", typeof(date), getvardims(date) )
  filevarattdef( NCDF, "date", date )

  filevardef   ( NCDF, "datesec", typeof(datesec), getvardims(time) )
  filevarattdef( NCDF, "datesec", datesec )

  filevardef   ( NCDF, "date_frac", typeof(date_frac), getvardims(time) )
  filevarattdef( NCDF, "date_frac", date_frac )

  filevardef   ( NCDF, "SEAICE", typeof(ice), getvardims(ice))
  filevarattdef( NCDF, "SEAICE", ice )

  setfileoption(NCDF,"DefineMode",False)       ; (not really necessary)
                                            
  NCDF->time     = (/time  /)  
  NCDF->lat      = (/lat   /)   
  NCDF->lon      = (/lon   /)   
  NCDF->date     = (/date  /)   
  NCDF->datesec  = (/datesec    /)   
  NCDF->date_frac= (/date_frac  /)   
  NCDF->SEAICE   = (/ice   /) 
end if

end
"END_MAIN_NCL"

# ===========================FORTRAN===================================
 
cat >! sstice.f << "END_SSTICE"
C NCLFORTSTART
      subroutine sstice (dirm, diri, film,  fildata
     +                  ,yyyy, mm, nlat, mlon, ice, sst, tagls, zmsg )
      implicit      none
      character*(*) diri, dirm, film, fildata 
      integer       yyyy, mm, nlat, mlon
      real          ice(mlon,nlat), sst(mlon,nlat), zmsg
      real          tagls(mlon,nlat)
C NCLEND
c *************************************************************
c  This section reads in the NCEP OI data [Jim Hurrell  read.f]
c *************************************************************

c
c  This subroutine reads a individual yyyymm reynolds OI SST fields
c
c  The geo-location of the SST array elements are:
c    SST(1,1)     = 0.5E, 89.5S
c    SST(1,2)     = 0.5E, 88.5S
c    SST(2,1)     = 1.5E, 89.5S
c    SST(360,180) = 359.5E, 89.5N
c
c    a land/sea mask should be used to mask out OI SST analyzed values
c    not located in the ocean, e.g. data file lstags.onedeg.dat
c    land=0  ocean=1
c
c  sst   - sea surface temperature array (deg C)
c  ice   - ice concentration array (%)  (0-100,   >100 = land or coast)
c  iyrst - year of start date of analysis
c  imst  - month of start date of analysis
c  idst  - day of start date of analysis
c  iyrnd - year of end date of analysis
c  imnd  - month of end date of analysis
c  idnd  - day of end date of analysis
c  ndays - number of days in analysis (start date thru enddate)
c  index - analysis version for reference
c  xlon  - longitude of center of grid square
c  xlat  - latitude of center of grid square
c  tagls - land/sea tag array (0=land, 1=water)

c  NOTES:  
c   - land values for sst do not necessarily coincide with land values 
c       from ice analysis

      integer    id, jd, krecl
      parameter (id=360,jd=180)
      parameter (krecl=id*jd*4)

      character*1  cice(id,jd)
      character*6  yyyymm

c c c real    sstnew(id,jd),icenew(id,jd), tagls(id,jd)
c c c real    sstnew(id,jd),icenew(id,jd)
      real*8  ix 
c c c real    sstmax(91),icemax(91)
      real    si, sstm

      integer iyrst,imst,idst,iyrnd,imnd,idnd,ndays,index
      integer i,j,ientry, ic
      data    ientry /0/
      save    ientry
c c c save    tagls 
c upon initial entry: open and read the land / sea mask

c *************************************************************
c  This section reads in the land-sea mask.
c  The OI analysis is done over all ocean areas and 
c  the Great Lakes. There is no analysis over land.  
c  The land values are filled by a Cressman interpolation 
c  to produce a complete grid for possible interpolation to 
c  other grids.  
c  ocean: lstag=1         land: lstag=0
c *************************************************************
c  For the ice fields, the value 122 represents land or coast.  
c  Note, the ice land mask is a function of the ice analysis 
c  and may change periodically.
c *************************************************************
      print *, "ENTER SSTICE"

      if (ientry.eq.0) then
          print *,"dirm//film=",trim(dirm)//trim(film)
          open(50,file=trim(dirm)//trim(film), convert="big_endian",
     *         form='unformatted',access='direct',recl=krecl)
          print *,"Past open statement for mask"
          read(50,rec=1) tagls
          print *,"Past read statement"
          close(50)
          ientry = 1
          print *, "=> LAND_SEA MASK READ OK  <="
      end if
      print *, "AFTER IENTRY:*************************"

c *************************************************************
c  This section reads in the data [Jim Hurrell  read.f]
c *************************************************************

c read sst and sea-ice concentration data

      print *,"diri=",diri
      print *,"fildata=",fildata
      open(10,file=diri//fildata, convert="big_endian",
     *        form='unformatted')
      print *, "AFTER OPEN"
      read(10) iyrst,imst,idst,iyrnd,imnd,idnd,ndays,index
      print *, "=>fortran: ",iyrst,imst,idst,iyrnd,imnd,idnd,ndays,index
      read(10) ((sst(i,j),i=1,id),j=1,jd)
      read(10) ((cice(i,j),i=1,id),j=1,jd)
      close (10)

      do j = 1, jd
        do i = 1, id
           ice(i,j) = float(ichar(cice(i,j)))
        enddo
      enddo

c set SST and ice to missing over land

      do j = 1, jd
        do i = 1, id
           if (tagls(i,j).eq.0) sst(i,j) = zmsg
           if (ice(i,j).gt.100.) ice(i,j) = zmsg
        enddo
      enddo

      return
      end
"END_SSTICE"

# consistent.f and coast_land_NearNbor.f same as 
# CGD: /fs/cgd/home0/shea/ncld/ncld2/ncld3/hadley 
 
## Portland Group compiler no longer available
##WRAPIT -d -pg -fPIC sstice.f 
##WRAPIT -d -pg -fPIC consistent.f 
## WRAPIT -d -pg -fPIC coast_land_NearNbor.f

WRAPIT -d sstice.f 
WRAPIT -d consistent.f 
WRAPIT -d coast_land_NearNbor.f

# ============================Execute===================================

 ncl main.ncl 

# ============================Clean UP==================================
 /bin/rm -f  main.ncl    # this is local
 /bin/rm -f  coast_land_NearNbor*o
 /bin/rm -f  consistent*o
 /bin/rm -f  objects
 /bin/rm -f  sstice.f    # this is local
 /bin/rm -f  sstice*o
 /bin/rm -f  WRAPIT* 
 /bin/rm -f  core

# ============================copy files================================

#scp /ptmp/shea/SSTICE/*nc tramhill.cgd.ucar.edu:/project/cas/shea/hadley/.  
exit 
