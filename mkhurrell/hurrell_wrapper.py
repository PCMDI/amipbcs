#!/bin/env python
# -*- coding: utf-8 -*-

"""

Stephen Po-Chedley 15 November 2018

Hurrell boundary conditions example. 

@author: pochedls
"""

import hurrellfx
import cdms2
import numpy as np

# set variables
# tos_file = '/p/user_pub/work/input4MIPs/CMIP6/CMIP/PCMDI/PCMDI-AMIP-1-1-3/ocean/mon/tos/gn/v20171031/tos_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc'
tos_file = '/p/user_pub/work/input4MIPs/CMIP6/CMIP/PCMDI/PCMDI-AMIP-1-1-3/seaIce/mon/siconc/gn/v20171031/siconc_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc'
fnOut = 'test.nc'
nyears = 5

# load data
f = cdms2.open(tos_file)
tos = f('siconc')
# tos = f('tos')
f.close()

# set target grid if desired
latn = np.arange(-89, 90, 2)
lonn = np.arange(1, 360, 2)
latAxis = cdms2.createAxis(latn)
latAxis.id = 'latitude'
latAxis.units = 'degrees_north'
latAxis.designateLatitude()
lonAxis = cdms2.createAxis(lonn)
lonAxis.id = 'longitude'
lonAxis.units = 'degrees_east'
lonAxis.designateLongitude()
lonAxis.isCircularAxis()
lonAxis.isCircular()
targetGrid = cdms2.createRectGrid(latAxis, lonAxis, order='yx', type='generic')

# option to add in mask
# fn_mask = '/work/cmip-dyn/CMIP5/fx/sftlf/CMIP5.CMIP.amip.NOAA-GFDL.GFDL-HIRAM-C180.r0i0p0.fx.sftlf.atmos.glb-2d-gu.v1.0000000.0.xml'
# f = cdms2.open(fn_mask)
# sftlf = f('sftlf')
# diag = {}
# sftlf = sftlf.regrid(targetGrid, regridTool='esmf', regridMethod = 'linear', missing=np.nan, coordSys='deg', diag = diag, periodicity = 1) 
# sftlf[sftlf > 75] = -999
# sftlf[sftlf >= 0] = 1
# f.close()

# create midpoint values
# tosp = hurrellfx.createMonthlyMidpoints(tos, 'sst', 'celcius', nyears)
tosp = hurrellfx.createMonthlyMidpoints(tos, 'ice', 'percent', nyears)

# save output file
f = cdms2.open(fnOut, 'w')
f.write(tosp)
f.close()


