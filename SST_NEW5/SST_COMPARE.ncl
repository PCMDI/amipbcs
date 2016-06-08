; =================================================
; This script will check the overlay years on
; the update file to make sure the differencess are 0.0
; For some reason the 1st two months show differences
; over land. No reason why!!!!
; =================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

                                 ; This is the file to be updated
 ;dira = "/project/cas/shea/hadley/" 
 ;dira = "/project/cas/shea/SSTICE/"
  dira = "../"
 ;fila = "MODEL.SST.HAD187001-198110.OI198111-201103.nc"
 ;fila = "MODEL.SST.HAD187001-198110.OI198111-201112.nc"
 ;fila = "MODEL.SST.HAD187001-198110.OI198111-201203.nc"
  fila = "MODEL.SST.HAD187001-198110.OI198111-201403.nc"
  print("fila="+fila)
                                 ; This is the raw update file
                                 ; This includes a few overlap years
 ;dirb = "/project/cas/shea/SSTICE/"
  dirb = "./"
 ;filb = systemfunc("cd "+dirb+" ; ls MODEL.OI2.sst.mnly.*.unf.nc")
 ;filb = "MODEL.OI2.sst.mnly.20xxyy-20XXYY.unf.nc"
 ;filb = "MODEL.OI2.sst.mnly.201001-201108.unf.nc"
 ;filb = "MODEL.OI2.sst.mnly.201001-201203.unf.nc"
 ;filb = "MODEL.OI2.sst.mnly.201201-201303.unf.nc"
  filb = "MODEL.OI2.sst.mnly.201401-201503.unf.nc"
  print("filb="+filb)
     
  fa   = addfile(dira+fila, "r")
  fb   = addfile(dirb+filb, "r")

  datea = fa->date
  dateb = fb->date

  datea = datea/100       ; yyyymm
  dateb = dateb/100
  na    = dimsizes(datea)
  nb    = dimsizes(dateb)

  idatea = ind(datea.ge.dateb(0) .and. datea.le.datea(na-1))
  idateb = ind(dateb.ge.dateb(0) .and. dateb.le.datea(na-1))
  print(datea(idatea)+"   "+dateb(idateb))

  ssta   = fa->SST(idatea,:,:)
  sstb   = fb->SST(idateb,:,:)

  printVarSummary(ssta)
  printVarSummary(sstb)
 
  printMinMax(ssta, True)
  printMinMax(sstb, True)
  print("=================================")
  print("=================================")
  print("=================================")

  dims   = dimsizes(ssta)
  print(dims)
  ntim   = dims(0)

  diff   = sstb-ssta
  mxdiff = max(abs(diff))
  print("Max diff over ALL dates: ="+mxdiff)

  do nt=0,ntim-1
     print(dateb(idateb(nt)) +"   mxdiff="+ max(abs(diff(nt,:,:))))
  end do

  diff@long_name = "SST Diff"
  copy_VarCoords(ssta, diff)
  
;************************************************
; create plot
;************************************************
  ;wks = gsn_open_wks("x11","TEST")                 ; open a ps file
  wks = gsn_open_wks("pdf","TEST")                 ; open a ps file
  res = True                                    ; plot mods desired
  res@gsnMaximize          = True
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       =  -4.               ; set min contour level
  res@cnMaxLevelValF       =   4.               ; set max contour level
  res@cnLevelSpacingF      =   0.5              ; set contour spacing

 ;res@tiMainString         = "CCM2 T42 July"    ; plot title
 ;res@cnInfoLabelOrthogonalPosF = -0.07         ; move the label inside th plot

  do nt=0,ntim-1
     res@gsnCenterString   =  dateb(idateb(nt))
     plot = gsn_csm_contour_map_ce(wks,diff(nt,:,:), res)  ; create plot
  end do

