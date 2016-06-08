#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:08:03 2015

Paul J. Durack 16th Jun 2015

This script sanitizes nc files written by the lats library

PJD 16 Jun 2015     - Started
PJD 16 Jun 2015     - Finalised outputs - organized by variables and subdir
PJD 19 Jun 2015     - Added tarfile creation
PJD 19 Jun 2015     - Added sanitize/package control options
PJD 22 Jun 2015     - Following zos_AVISO obs4MIPS data added and further clarified global_atts
PJD 22 Jun 2015     - Updated source_id = 'PCMDI-AMIPBCS-1.0' -> 'PCMDI-AMIP-1.0.0'
PJD 27 May 2016     - Renamed sanitizeAndPackage.py -> sanitize.py
PJD 31 May 2016     - Updated to follow CMIP6 data conventions
PJD  1 Jun 2016     - Created areacello variable
PJD  2 Jun 2016     - Further standardised variables for consistency across output files
PJD  2 Jun 2016     - Incremented version number v1.0.1 -> v1.1.0
PJD  3 Jun 2016     - Updated to use netcdf4 in CMOR3
PJD  3 Jun 2016     - Updated areacello to include landsea mask
PJD  7 Jun 2016     - Forced 'gregorian' input calendar type (check old vs new calendars)
PJD  7 Jun 2016     - Updated to CMOR 3.0.6 and using cmor.set_deflate for deflation
PJD  8 Jun 2016     - Updated to generate calendar using makeCalendar (rather than the input data time axis)
                      solves issue with hour offset between cdms and cmor written files
PJD  8 Jun 2016     - Update to write bcs vars with no time_bounds
PJD  8 Jun 2016     - Turned off landsea/masking in areacello to maintain consistency with the other variables
                    - TODO:

@author: durack1
"""

import cmor,datetime,gc,glob,os,pytz ; #pdb
import cdms2 as cdm
import cdutil as cdu
import numpy as np
from durolib import globalAttWrite,makeCalendar,mkDirNoOSErr

#%% Set cdms format criteria
cdm.setNetcdfClassicFlag(0) ; # 1 = amipbc files 437.2MB; amipobs files 437.3MB
cdm.setNetcdfDeflateFlag(1) ; # 1 = amipbc files sic 75.2MB, tos 239.2MB; amipobs files sic 29.3MB, tos 225.9MB

#%% Set version info
activity_id         = 'input4MIPs'
comment             = 'Based on Hurrell SST/sea ice consistency criteria applied to merged HadISST (1870-01 1981-10) & NCEP-0I2 (1981-11 to 2015-12)' ; # WILL REQUIRE UPDATING
contact             = 'pcmdi-cmip@lists.llnl.gov'
dataVer             = 'PCMDI-AMIP-1-1-0' ; # WILL REQUIRE UPDATING
dataVerSht          = 'v1.1.0' ; # WILL REQUIRE UPDATING
data_structure      = 'grid'
further_info_url    = 'http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php' ; # WILL REQUIRE UPDATING - point to GMD paper when available
institute_id        = 'PCMDI'
institution         = 'Program for Climate Model Diagnosis and Intercomparison (LLNL), Livermore, CA, USA'
last_year           = '2015' ; # WILL REQUIRE UPDATING
license_txt         = 'AMIP boundary condition data produced by PCMDI is licensed under a Creative Commons Attribution \"Share Alike\" 4.0 International License (http://creativecommons.org/licenses/by/4.0/). The data producers and data providers make no warranty, either express or implied, including but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law.'
mip_specs           = 'AMIP CMIP5 CMIP6'
project_id          = 'AMIP'
ref_obs             = 'Hurrell, J. W., J. J. Hack, D. Shea, J. M. Caron, and J. Rosinski (2008) A New Sea Surface Temperature and Sea Ice Boundary Dataset for the Community Atmosphere Model. J. Climate, 22 (19), pp 5145-5153. doi: 10.1175/2008JCLI2292.1'
ref_bcs             = 'Taylor, K.E., D. Williamson and F. Zwiers, 2000: The sea surface temperature and sea ice concentration boundary conditions for AMIP II simulations. PCMDI Report 60, Program for Climate Model Diagnosis and Intercomparison, Lawrence Livermore National Laboratory, 25 pp. Available online: http://www-pcmdi.llnl.gov/publications/pdf/60.pdf'
source              = 'PCMDI-AMIP 1.1.0: Merged SST based on UK MetOffice HadISST and NCEP OI2' ; # WILL REQUIRE UPDATING
time_period         = '187001-201512' ; # WILL REQUIRE UPDATING

#%% Set directories
homePath    = '/work/durack1/Shared/150219_AMIPForcingData/'
sanPath     = os.path.join(homePath,'_'.join(['360x180',dataVerSht,'san']))

#%% Get files
newList = sorted(glob.glob(os.path.join(homePath,'_'.join(['360x180',dataVerSht]),'*.nc')))

#%% Create memory objects for building composite field
fileCount = 0
for fileName in newList:
    if 'amipbc_sst' in fileName:
        fileCount = fileCount + 1
fileCount   = fileCount - 1
varComp     = np.ma.zeros([fileCount*12,180,360])
timeComp    = np.zeros([fileCount*12])

#%% Get variable into memory
count = 0
for filePath in newList:
    obsVsBC     = filePath.split('/')[-1].split('_')[0]
    varName     = filePath.split('/')[-1].split('_')[1]
    climCheck   = filePath.split('/')[-1].split('_')[-1]
    if obsVsBC == 'bcinfo' or obsVsBC == 'spinup' or climCheck == 'clim.nc':
        print 'Invalid file, skipping..'
        continue
    if 'bc' in obsVsBC:
        BC = 'bcs'
        cdm.setAutoBounds(2) ; # Turn off time_bounds - only latitude/longitude written
    else:
        BC = ''
        cdm.setAutoBounds(1) ; # Turn on time_bounds - time/latitude/longitude written
    if varName == 'sst':
        varLoad = ''.join(['tos',BC])
        varPath = 'tos'
    else:
        varLoad = ''.join([varName,BC])
        varPath = varName
    print filePath
    fH      = cdm.open(filePath)
    var     = fH(varLoad)
    varLen  = var.shape[0]

    #%% Cleanup coord atts
    # time
    time                    = var.getAxis(0)
    time.standard_name      = 'time'
    time.long_name          = 'time'
    time.calendar           = 'gregorian' ; # Force Gregorian
    time.axis               = 'T'
    time.toRelativeTime('days since 1870-1-1') ; # Fix negative values
    cdu.setTimeBoundsMonthly(time) ; # Resolve issues with bounds being mid-time values rather than month-end/start values
    if BC == 'bcs':
        time._bounds_ = None ; # Required to purge bounds created by cdu call above

    #%% Write timestep to composite variable
    countUp                     = count + 12
    varComp[count:countUp]      = var
    timeComp[count:countUp]     = time
    count                       = countUp

    #%% Cleanup coord atts and create areacello
    if last_year in filePath:
        # latitude
        latitude                = var.getAxis(1)
        latitude.id             = 'lat'
        latitude.standard_name  = 'latitude'
        latitude.long_name      = 'latitude'
        latitude.axis           = 'Y'
        delattr(latitude,'realtopology')
        # longitude
        longitude               = var.getAxis(2)
        longitude.id            = 'lon'
        longitude.standard_name = 'longitude'
        longitude.long_name     = 'longitude'
        longitude.axis          = 'X'
        delattr(longitude,'realtopology')

        # Create areacello variable
        if 'amipbc_sic' in filePath:
            areacello                   = cdu.area_weights(var[0,]) ; # areacello.sum() = 1.0
            earthSurfaceArea            = 510.1 ; # million km2
            earthSurfaceAreaM2          = earthSurfaceArea*1e9 ; # m2
            areacelloM2                 = areacello*earthSurfaceAreaM2
            areacelloM2.standard_name   = 'cell_area'
            areacelloM2.long_name       = 'Ocean Grid-Cell Area'
            areacelloM2.units           = 'm2'
            areacelloM2.id              = 'areacello'
            # pdb.set_trace()
            # Set landsea mask
            #areacelloM2.mask            = cdu.generateLandSeaMask(var) ; # Enforce landsea mask on areacello
            #areacelloM2.missing         = 1e20
            # Create output file
            outFileA = os.path.join(sanPath,filePath.split('/')[-1])
            outFileA = outFileA.replace('_'.join(['',last_year]),'')
            outFileA = outFileA.replace('sst','areacello')
            outFileA = outFileA.replace('sic','areacello')
            print 'Processing: ',outFileA
            if not os.path.exists(sanPath):
                mkDirNoOSErr(sanPath)
            if os.path.exists(outFileA):
                os.remove(outFileA)
            fO = cdm.open(outFileA,'w')
            # global atts
            globalAttWrite(fO,options='noid'); # 'noid' option prevents data_contact and institution being written
            fO.Conventions      = 'CF-1.6' ; fO.sync()
            fO.title            = areacelloM2.long_name
            fO.activity_id      = activity_id
            fO.further_info_url = further_info_url
            fO.comment          = comment ; fO.sync()
            fO.contact          = contact ; fO.sync() ; # Overwritten globalAttWrite
            local               = pytz.timezone("America/Los_Angeles")
            time_now            = datetime.datetime.now();
            local_time_now      = time_now.replace(tzinfo = local)
            utc_time_now        = local_time_now.astimezone(pytz.utc)
            time_format         = utc_time_now.strftime("%Y-%m-%dT%H:%M:%SZ")
            fO.creation_date    = time_format
            fO.data_structure   = data_structure
            fO.institute_id     = institute_id ; fO.sync()
            fO.institution      = institution ; fO.sync()
            fO.license          = license_txt
            fO.mip_specs        = mip_specs
            fO.project_id       = project_id
            fO.realm            = 'fx'
            fO.references       = ref_obs ; fO.sync()
            fO.source           = source
            fO.source_id        = dataVer
            fO.write(areacelloM2.astype('float32'));
            fO.close()
            del(areacello,areacelloM2,earthSurfaceArea,earthSurfaceAreaM2) ; gc.collect()
    del(var,time) ; gc.collect()
    fH.close()

    #%% Write to outfile
    if last_year in filePath:
        '''
        # Convert variables to cdms transient variables
        var     = cdm.createVariable(varComp)
        time    = cdm.createAxis(timeComp,id='time')
        # Assign new coord atts
        var.setAxis(0,time)
        var.setAxis(1,latitude)
        var.setAxis(2,longitude)
        time.standard_name      = 'time'
        time.long_name          = 'time'
        time.units              = 'days since 1870-1-1'
        time.calendar           = 'gregorian'
        time.axis               = 'T'
        '''
        time = makeCalendar('1870',last_year,calendarStep='months')
        var = cdm.createVariable(varComp)
        var.setAxis(0,time)
        var.setAxis(1,latitude)
        var.setAxis(2,longitude)

        if BC == 'bcs':
            var.cell_methods    = 'time: point'
            longTxt             = 'constructed mid-month'
            dataUsageTips       = 'The mid-month data should be linearly interpolated in time and then clipped for use as boundary conditions to drive AMIP simulations as described at: http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php'
            refTxt              = ref_bcs
        else:
            var.cell_methods    = 'time: mean'
            longTxt             = 'observed monthly mean'
            dataUsageTips       = 'The observed monthly-mean data should *NOT* be used to drive AMIP simulations. For further information see: http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php'
            refTxt              = ref_obs

        if varName == 'sst':
            var.id              = ''.join(['tos',BC])
            var.name            = ''.join(['tos',BC])
            var.long_name       = ' '.join(['AMIP',longTxt,'sea surface temperature'])
            var.standard_name   = 'sea_surface_temperature'
            var.units           = 'K'
            realmTxt            = 'ocean'
        elif varName == 'sic':
            var.id              = ''.join(['siconc',BC])
            var.name            = ''.join(['siconc',BC])
            var.long_name       = ' '.join(['AMIP',longTxt,'sea ice area fraction'])
            var.standard_name   = 'sea_ice_area_fraction'
            var.units           = '%'
            realmTxt            = 'seaIce'

        #%% Create output file
        outFile = os.path.join(sanPath,filePath.split('/')[-1])
        outFile = outFile.replace(last_year,time_period)
        outFile = outFile.replace('sst','tos')
        outFile = outFile.replace('sic','siconc')
        print 'Processing: ',outFile
        if not os.path.exists(sanPath):
            mkDirNoOSErr(sanPath)
        if os.path.exists(outFile):
            os.remove(outFile)
        fO = cdm.open(outFile,'w')
        # global atts
        globalAttWrite(fO,options='noid'); # 'noid' option prevents data_contact and institution being written
        fO.Conventions      = 'CF-1.6' ; fO.sync()
        fO.title            = var.long_name
        fO.activity_id      = activity_id
        fO.further_info_url = further_info_url
        fO.comment          = comment ; fO.sync()
        fO.contact          = contact ; fO.sync() ; # Overwritten globalAttWrite
        local               = pytz.timezone("America/Los_Angeles")
        time_now            = datetime.datetime.now();
        local_time_now      = time_now.replace(tzinfo = local)
        utc_time_now        = local_time_now.astimezone(pytz.utc)
        time_format         = utc_time_now.strftime("%Y-%m-%dT%H:%M:%SZ")
        fO.creation_date    = time_format
        fO.data_structure   = data_structure
        fO.data_usage_tips  = dataUsageTips ; fO.sync()
        fO.frequency        = 'mon' ; fO.sync()
        fO.institute_id     = institute_id ; fO.sync()
        fO.institution      = institution ; fO.sync()
        fO.license          = license_txt
        fO.mip_specs        = mip_specs
        if BC == '':
            fO.product      = 'observations'
        else:
            fO.product      = 'derived'
        fO.project_id       = project_id
        fO.realm            = realmTxt
        fO.references       = refTxt ; fO.sync()
        fO.source           = source
        fO.source_id        = dataVer
        fO.write(var.astype('float32'));
        del(time,latitude,longitude) ; gc.collect()
        fO.close()

        #%% CMORise
        print 'CMOR start'
        '''
        Compression: deflate=1,deflate_level=x,shuffle=1 (cdms2DefaultDeflate/CMOR3NoDeflate)
        areacello   1: 55.6KB, 3: xx.xKB, 9: xx.xKB ( 36.1KB/303.7KB)
        siconc      1: 44.6MB, 3: 43.2MB, 9: 39.9MB ( 29.3MB/433.1MB)
        siconcbcs   1: 73.0MB, 3: 72.3MB, 9: 69.6MB ( 75.2MB/433.1MB)
        tos         1:177.1MB, 3:174.8MB, 9:171.3MB (225.9MB/433.1MB)
        tosbcs      1:185.1MB, 3:182.6MB, 9:178.4MB (239.2MB/433.1MB)
        '''
        # Write bcs variables
        if 'amipbc' in filePath:
            cmor.setup(inpath='CMOR/input4MIPs-cmor-tables/Tables',netcdf_file_action=cmor.CMOR_REPLACE_4)
            cmor.dataset_json("CMOR/drive_input4MIPs_bcs.json")
            f       = cdm.open(outFile)
            d       = f[var.id]
            lat     = d.getLatitude()
            lon     = d.getLongitude()
            time    = d.getTime()
            cmor.set_cur_dataset_attribute('history',f.history) ; # Force local file attribute as history
            table   = 'input4MIPs.json'
            cmor.load_table(table)
            axes    = [ {'table_entry': 'time2',
                         'units': 'days since 1870-01-01',
                         },
                         {'table_entry': 'latitude',
                          'units': 'degrees_north',
                          'coord_vals': lat[:],
                          'cell_bounds': lat.getBounds()},
                         {'table_entry': 'longitude',
                          'units': 'degrees_east',
                          'coord_vals': lon[:],
                          'cell_bounds': lon.getBounds()},
                      ]
            axis_ids = list()
            for axis in axes:
                axis_id = cmor.axis(**axis)
                axis_ids.append(axis_id)
            varid   = cmor.variable(var.id,var.units,axis_ids)
            values  = np.array(d[:],np.float32)
            cmor.set_deflate(varid,1,1,1) ; # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
            cmor.write(varid,values,time_vals=time[:])
            f.close()
            cmor.close()
            # Cleanup
            del(outFile,var,f,d,lat,lon,time) ; gc.collect()

      # Write obs variables
        if 'amipobs' in filePath:
            cmor.setup(inpath='CMOR/input4MIPs-cmor-tables/Tables',netcdf_file_action=cmor.CMOR_REPLACE_4)
            cmor.dataset_json("CMOR/drive_input4MIPs_obs.json")
            f       = cdm.open(outFile)
            d       = f[var.id]
            lat     = d.getLatitude()
            lon     = d.getLongitude()
            time    = d.getTime()
            cmor.set_cur_dataset_attribute('history',f.history) ; # Force local file attribute as history
            table   = 'input4MIPs.json'
            cmor.load_table(table)
            axes    = [ {'table_entry': 'time',
                         'units': 'days since 1870-01-01',
                         },
                         {'table_entry': 'latitude',
                          'units': 'degrees_north',
                          'coord_vals': lat[:],
                          'cell_bounds': lat.getBounds()},
                         {'table_entry': 'longitude',
                          'units': 'degrees_east',
                          'coord_vals': lon[:],
                          'cell_bounds': lon.getBounds()},
                      ]
            axis_ids = list()
            for axis in axes:
                axis_id = cmor.axis(**axis)
                axis_ids.append(axis_id)
            varid   = cmor.variable(var.id,var.units,axis_ids)
            values  = np.array(d[:],np.float32)
            cmor.set_deflate(varid,1,1,1) ; # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
            #cmor.write(varid,values,time_vals=d.getTime()[:],time_bnds=d.getTime().genGenericBounds()) ; # Not valid for time
            cmor.write(varid,values,time_vals=time[:],time_bnds=time.getBounds())
            f.close()
            cmor.close()
            # Cleanup
            del(outFile,var,f,d,lat,lon,time) ; gc.collect()

        # Write areacello variable
        if 'amipbc_sic' in filePath:
            cmor.setup(inpath='CMOR/input4MIPs-cmor-tables/Tables',netcdf_file_action=cmor.CMOR_REPLACE_4)
            cmor.dataset_json("CMOR/drive_input4MIPs_obs.json")
            f       = cdm.open(outFileA)
            d       = f['areacello']
            lat     = d.getLatitude()
            lon     = d.getLongitude()
            time    = d.getTime()
            cmor.set_cur_dataset_attribute('history',f.history) ; # Force local file attribute as history
            table   = 'input4MIPs.json'
            cmor.load_table(table)
            #cmor.set_table(tables[1]) ; # Change back to input4MIPs table
            axes    = [ {'table_entry': 'latitude',
                         'units': 'degrees_north',
                         'coord_vals': lat[:],
                         'cell_bounds': lat.getBounds()},
                        {'table_entry': 'longitude',
                         'units': 'degrees_east',
                         'coord_vals': lon[:],
                         'cell_bounds': lon.getBounds()},
                      ]
            axis_ids = list()
            for axis in axes:
                axis_id = cmor.axis(**axis)
                axis_ids.append(axis_id)
            varid   = cmor.variable('areacello',d.units,axis_ids,missing_value=1e20)
            values  = np.array(d[:],np.float32)
            cmor.set_deflate(varid,1,1,1) ; # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
            cmor.write(varid,values)
            f.close()
            cmor.close()
            # Cleanup
            del(outFileA,f,d,lat,lon,time) ; gc.collect()
        # Cleanup
        del(table,axes,axis_ids,varid,values) ; gc.collect()

        #%% Reset variables
        varComp     = np.ma.zeros([fileCount*12,180,360])
        timeComp    = np.zeros([fileCount*12])
        count       = 0