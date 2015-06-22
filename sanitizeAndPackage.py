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

@author: durack1
"""

import datetime,gc,glob,os,pytz,re,tarfile
import cdms2 as cdm
import cdutil as cdu
from durolib import globalAttWrite,mkDirNoOSErr

#%% Set process options
sanitize = True
package = True

#%% Set cdms format criteria
cdm.setNetcdfClassicFlag(1)
cdm.setNetcdfDeflateFlag(0)

#%% Set directories
homePath    = '/work/durack1/Shared/150219_AMIPForcingData/'
sanPath     = os.path.join(homePath,'360x180_san')

#%% Get files
newList = sorted(glob.glob(os.path.join(homePath,'360x180/*.nc')))

#%% Get variable into memory
if sanitize:
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
        fH  = cdm.open(filePath)
        var = fH(varLoad)
    
        #%% Cleanup coord atts
        # time
        time                    = var.getAxis(0)
        time.standard_name      = 'time'
        time.long_name          = 'time'
        time.axis               = 'T'
        cdu.setTimeBoundsMonthly(time) ; # Resolve issues with bounds being mid-time values rather than month-end/start values
        time.toRelativeTime('days since 1870-1-1') ; # Fix negative values
        if BC == 'bcs':
            time._bounds_ = None ; # Required to purge bounds created by cdu call above
        # latitude
        latitude                = var.getAxis(1)
        latitude.standard_name  = 'latitude'
        latitude.axis           = 'Y'
        # longitude
        longitude               = var.getAxis(2)
        longitude.standard_name = 'longitude'
        longitude.axis          = 'X'
    
        #%% Assign new coord atts
        var.setAxis(0,time)
        var.setAxis(1,latitude)
        var.setAxis(2,longitude)
    
        #%% Cleanup variable atts - dependent on source data
        delattr(var,'comments')
        delattr(var,'grid_name')
        delattr(var,'grid_type')
        delattr(var,'level_description')
        delattr(var,'time_statistic')
        if BC == 'bcs':
            var.cell_methods    = 'time: point'
            longTxt             = 'constructed mid-month'
            dataUsageTips       = 'The mid-month data should be linearly interpolated in time and then clipped for use as boundary conditions \
to drive AMIP simulations as described at: http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php'
            refTxt              = 'Taylor, K.E., D. Williamson and F. Zwiers, 2000: The sea surface temperature and sea ice concentration boundary conditions \
for AMIP II simulations. PCMDI Report 60, Program for Climate Model Diagnosis and Intercomparison, Lawrence Livermore National Laboratory, 25 pp. \
Available online: http://www-pcmdi.llnl.gov/publications/pdf/60.pdf'
        else:
            var.cell_methods    = 'time: mean'
            longTxt             = 'observed monthly mean'
            dataUsageTips       = 'The observed monthly-mean data should *NOT* be used to drive AMIP simulations. For further information see: \
http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php'
            refTxt              = 'Hurrell, J. W., J. J. Hack, D. Shea, J. M. Caron, and J. Rosinski (2008) A New Sea Surface Temperature and Sea Ice Boundary \
Dataset for the Community Atmosphere Model. J. Climate, 22 (19), pp 5145-5153. doi: 10.1175/2008JCLI2292.1'
    
        if varName == 'sst':
            var.name            = ''.join(['tos',BC])
            var.long_name       = ' '.join(['AMIP',longTxt,'sea surface temperature'])
            var.standard_name   = 'sea_surface_temperature'
            var.units           = 'K'
            realmTxt            = 'ocean'
        elif varName == 'sic':
            var.name            = ''.join(['sic',BC])
            var.long_name       = ' '.join(['AMIP',longTxt,'sea ice area fraction'])
            var.standard_name   = 'sea_ice_area_fraction'
            var.units           = '%'
            realmTxt            = 'seaIce'
    
        #%% Create output file
        outFile = os.path.join(sanPath,varPath,filePath.split('/')[-1])
        print 'Processing: ',outFile
        if not os.path.exists(os.path.join(sanPath,varPath)):
            mkDirNoOSErr(os.path.join(sanPath,varPath))
        if os.path.exists(outFile):
            os.remove(outFile)
        fO = cdm.open(outFile,'w')
        # global atts
        globalAttWrite(fO,options='noid'); # 'noid' option prevents data_contact and institution being written
        fO.Conventions      = 'CF-1.6' ; fO.sync()
        fO.comment          = 'Based on Hurrell SST/sea ice consistency criteria applied to merged HadISST (1870-01 1981-10) & NCEP-0I2 (1981-11 to 2015-03)' ; fO.sync()
        fO.contact          = 'pcmdi-cmip@lists.llnl.gov' ; fO.sync() ; # Overwritten globalAttWrite
        local               = pytz.timezone("America/Los_Angeles")
        time_now            = datetime.datetime.now();
        local_time_now      = time_now.replace(tzinfo = local)
        utc_time_now        = local_time_now.astimezone(pytz.utc)
        time_format         = utc_time_now.strftime("%Y-%m-%dT%H:%M:%SZ")
        fO.creation_date    = time_format
        fO.data_structure   = 'grid'
        fO.data_usage_tips  = dataUsageTips ; fO.sync()
        fO.frequency        = 'mon' ; fO.sync()
        fO.institute_id     = 'PCMDI' ; fO.sync()
        fO.institution      = 'Program for Climate Model Diagnosis and Intercomparison (LLNL), Livermore, CA, U.S.A.' ; fO.sync()
        fO.mip_specs        = 'AMIP CMIP5 CMIP6'
        if BC == '':
            fO.product      = 'observations'
        fO.project_id       = 'AMIP'
        fO.realm            = realmTxt
        fO.references       = refTxt ; fO.sync()
        fO.source           = 'Merged SST based on UK MetOffice HadISST and NCEP OI2'
        fO.source_id        = 'PCMDI-AMIPBCS-1.0'
        fO.write(var.astype('float32'));
        del(var,time,outFile,latitude,longitude) ; gc.collect()
        fO.close()
        fH.close()

#%% Now wrap files into tarballs for distribution
if package:
    newList = sorted(glob.glob(os.path.join(homePath,'360x180_san/*/*.nc')))
    sicList,sicbcList,tosList,tosbcList = [[] for _ in range(4)]
    for i,j in enumerate(newList):
        if 'sic/amipbc' in j:
            sicbcList.append(j)
        elif 'sic/amipobs' in j:
            sicList.append(j)
        elif 'tos/amipbc' in j:
            tosbcList.append(j)
        elif 'tos/amipobs' in j:
            tosList.append(j)
    del(newList) ; gc.collect()
    
    fileLists = [sicbcList,sicList,tosbcList,tosList]
    fileNames = ['amipbcs_sic','amipobs_sic','amipbcs_tos','amipobs_tos']
    
    test1800s = re.compile('18[0-9]{2}.nc')
    test1900s = re.compile('19[0-4][0-9].nc')
    test1950s = re.compile('19[5-9][0-9].nc')
    test2000s = re.compile('20[0-9]{2}.nc')
    
    tests       = [test1800s,test1900s,test1950s,test2000s]
    testNames   = ['1870-1899','1900-1949','1950-1999','2000-2015']
    
    # Change to output dir
    os.chdir(os.path.join(homePath,'360x180_san'))
    print os.getcwd()
    for count1,test in enumerate(tests):
        for count2,files in enumerate(fileLists):
            #print ''.join([testNames[count1],': ',test.pattern])
            #print ''.join([fileNames[count2],'_',testNames[count1],'.tar.gz'])
            tarFile = ''.join([fileNames[count2],'_',testNames[count1],'.tar.gz'])
            filesToTar = []
            for filePath in files:
                if re.search(test,filePath):
                    filesToTar.append(os.path.join(filePath.split('/')[-2],filePath.split('/')[-1]))
                    #print os.path.join(filePath.split('/')[-2],filePath.split('/')[-1])
            print ''.join(['Processing: ',tarFile])
            if os.path.exists(tarFile):
                os.remove(tarFile)
            tarH = tarfile.open(tarFile, "w:gz")
            for fileName in filesToTar:
                #print 'Adding:',fileName
                tarH.add(fileName)
            tarH.close()
            print ''.join(['Complete:   ',tarFile])
