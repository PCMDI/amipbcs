#!/bin/env python
# -*- coding: utf-8 -*-

"""
Stephen Po-Chedley 15 November 2018

Hurrell boundary conditions example.

PJD 15 Jul 2019     - Updated to use PCMDI-AMIP-1-1-5 as input
PJD 15 Jul 2019     - Updated to read grid and sftof from previous files
PJD 15 Jul 2019     - Added output file purge if exists
PJD 15 Jul 2019     - Added PCMDI-AMIP-1-2-0 as additional input
PJD 15 Jul 2019     - Added siconc in addition to tos
PJD  8 Aug 2019     - Update for latest hurrellfx.py

@author: pochedls and durack1
"""

import hurrellfx
import cdms2 as cdm
import os
#import pdb
#import numpy as np

# Set netcdf file criterion - turned on from default 0s
cdm.setCompressionWarnings(0) ; # Suppress warnings
cdm.setNetcdfShuffleFlag(0)
cdm.setNetcdfDeflateFlag(1)
cdm.setNetcdfDeflateLevelFlag(9)
# Hi compression: XXMb file ; # Single tosbcs variable
# No compression: 314.8Mb

# set variables
tosFile = '/data_crunchy_oceanonly/oceanonly_work/durack1/Shared/150219_AMIPForcingData/input4MIPs/CMIP6/CMIP/PCMDI/PCMDI-AMIP-1-2-0/ocean/mon/tos/gn/v20190715/tos_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-2-0_gn_187001-201812.nc'
tosOut = 'tosbcs-v1-2-0.nc'
sicFile = '/data_crunchy_oceanonly/oceanonly_work/durack1/Shared/150219_AMIPForcingData/input4MIPs/CMIP6/CMIP/PCMDI/PCMDI-AMIP-1-2-0/seaIce/mon/siconc/gn/v20190715/siconc_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-2-0_gn_187001-201812.nc'
sicOut = 'siconcbcs-v1-2-0.nc'
nyears = 5
print('filenames set..')

# load data
f = cdm.open(tosFile)
tos = f('tos')
f.close()
f = cdm.open(sicFile)
sic = f('siconc')
f.close()
print('tos/sic read..')

# set target grid if desired
'''
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
'''
# Get target grid
targetGrid = tos.getGrid()
print('grid generated..')

# option to add in mask
fn_mask = '/data_crunchy_oceanonly/oceanonly_work/durack1/Shared/150219_AMIPForcingData/input4MIPs/CMIP6/CMIP/PCMDI/PCMDI-AMIP-1-2-0/ocean/fx/sftof/gn/v20190715/sftof_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-2-0_gn.nc'
f = cdm.open(fn_mask)
sftof = f('sftof')
f.close()
print('sftof read..')

# create tos midpoint values
print('Entering function..')
tosBcs = hurrellfx.createMonthlyMidpoints(tos, 'sst', 'celcius', nyears, 'tosbcs') #, grid=targetGrid, mask=sftof)
print('Exiting function..')
# save output file
if os.path.exists(tosOut):
    print("".join(['** File exists.. removing: ',tosOut,' **']))
    os.remove(tosOut)
f = cdm.open(tosOut, 'w')
f.write(tosBcs)
f.close()

# create siconc midpoint values
print('Entering function..')
sicBcs = hurrellfx.createMonthlyMidpoints(sic, 'ice', '%', nyears, 'siconcbcs') #, grid=targetGrid, mask=sftof)
print('Exiting function..')
# save output file
if os.path.exists(sicOut):
    print("".join(['** File exists.. removing: ',sicOut,' **']))
    os.remove(sicOut)
f = cdm.open(sicOut, 'w')
f.write(sicBcs)
f.close()