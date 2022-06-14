#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:08:03 2015

Paul J. Durack 16th Jun 2015

This script sanitizes nc files written by the lats library

"""
# 2015/16
"""
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
PJD  8 Jun 2016     - Corrected unit issue with areacello - km2 -> m2
PJD  6 Sep 2016     - Update to out-of-band 1.1.0a update (Peter Rensch - CSIRO)
PJD  6 Sep 2016     - Updated to deal with partial years (varLen)
PJD  6 Sep 2016     - Deal with amipbc_sst_360x180_v1.1.0a.out
PJD  7 Sep 2016     - Deal with makeCalendar quirks - off by one
PJD 20 Oct 2016     - Update to 1.1.1
"""
# 2017/18
"""
PJD 13 Apr 2017     - Update to 1.1.2 and generate sftof
PJD 17 Apr 2017     - Corrected sftof mask to be correct indexes
PJD 18 Apr 2017     - Updated CMOR to "2017.04.18.3.2.3" and corrected SImon table
PJD 19 Apr 2017     - Updated tos units from K -> degC
PJD  9 Oct 2017     - Update to 1.1.3
PJD 23 Oct 2017     - Updated input4MIPs tables
PJD 24 Oct 2017     - Further updates to input4MIPs tables
PJD 25 Oct 2017     - Further updates to input4MIPs tables (region)
PJD 25 Oct 2017     - Update to generate json publication files
PJD 30 Oct 2017     - Added 'a' to sanPath for testing
PJD 31 Oct 2017     - Removed 'a' for testing; json publication file tweaks
PJD 16 Apr 2018     - Update to 1.1.4
PJD 16 Apr 2018     - Updated masks file 170306 to 170425
PJD 16 Apr 2018     - Updated dataset_version_number -> source_version
PJD 17 Apr 2018     - Updated registered source_id, and drive_input4MIPs_obs2.json (migrate content into registration)
PJD 19 Apr 2018     - Updated to drive_input4MIPs_obs.json and fixed globalAttWrite bug
PJD 19 Apr 2018     - Updated input*.json files to point directly to css03
PJD 20 Apr 2018     - Updated input*.json files "outpath" points to /p/user_pub/work/input4MIPs
PJD 20 Apr 2018     - Updated pytz call; try CMOR331 env as realm issue in v332
PJD 20 Apr 2018     - Turn off json descriptor creation
PJD 26 Apr 2018     - Updated to use jsonWriteFile,washPerms
PJD 27 Apr 2018     - Updated CMOR input to write to /p/user_pub/work - opened host dir to world readable (ames4)
PJD 27 Apr 2018     - Updated input4MIPsFuncs.py library - utc/pytz update
"""
# 2019
"""
PJD  4 Jan 2019     - Updated to 1.1.5
PJD 17 Jan 2019     - Added last_month (and if statement) variable to deal with half vs full years
PJD 18 Jan 2019     - Added output file list for publication steps
PJD 22 Jan 2019     - Finalized json and output file logging for ESGF publication
PJD 22 Jan 2019     - Updated CMOR/drive_input4MIPs_bcs/obs.json for final production run
PJD 23 Jan 2019     - Updated input4MIPsFuncs.py to deal with ESGF publication log files
PJD 24 Jan 2019     - Updated to use revised input4MIPsFuncs.py (added dataVersion to json name)
PJD  2 Jul 2019     - Updated to 1.2.0
PJD  2 Jul 2019     - Updated durolib path
PJD 10 Jul 2019     - Updated input4MIPsFuncs path
PJD 10 Jul 2019     - Added source_version to CMOR/drive_input4MIPs_bcs/*.json files (removed this and target_mip)
PJD 15 Jul 2019     - Update jsonId to deal with revised input4MIPsFuncs update
PJD 20 Nov 2019     - Update to write out 2018 data to v1.1.6
PJD 20 Nov 2019     - Updated prints for py3
PJD 20 Nov 2019     - /p/user_pub/work needed perm updates to allow DRS writing (drwxrwxr-x - climatew)
"""
# 2021
"""
PJD 28 Jul 2021     - Update to latest data; Convert to use mkhurrell_wrapper.py; Update home path
PJD  1 Sep 2021     - Added history attribute to replicate previous files
PJD  2 Sep 2021     - Further testing and migration to /p/user_pub/climate_work/durack1
PJD  7 Sep 2021     - Updated source URL from ftp://ftp.emc.ncep.noaa.gov/cmb/sst/oimonth_v2/ to
                      ftp://ftp.cpc.ncep.noaa.gov/precip/PORT/sst/oimonth_v2/
PJD  9 Sep 2021     - Updated to latest August 2021 data
PJD  3 Nov 2021     - Update for latest code; Code pads 12-months to beginning and end so no truncation required
PJD  4 Nov 2021     - Update for latest September 2021 data
PJD 17 Nov 2021     - Switchout for latest October 2021 data
PJD 18 Nov 2021     - Evaluating impact of 202109 vs 202110 input data
PJD  2 Dec 2021     - Updates to hurrellfx.py and *sub.f - obs range checks
PJD  2 Dec 2021     - Renamed mkhurrell -> pcmdiAmipBcs; hurrellfx.py -> pcmdiAmipBcsFx.py
"""
"""
PJD  1 Feb 2022     - Updated to reflect v1.1.7 data not v1.2.0 (CMIP6 not CMIP6Plus)
PJD 14 Jun 2022     - Updated to reflect v1.1.8 data not v1.2.0 (CMIP6 not CMIP6Plus)
PJD 14 Jun 2022     - Corrected license to reflect CC BY 4.0 (was garbled before)
                    - TODO:
                    - Always check for group membership to climatew before running this, otherwise problems occur



OIv2 info
ftp://ftp.emc.ncep.noaa.gov/cmb/sst/oimonth_v2/OIv2_monthly.info.asc
Last Day of             Analysis
Month Falls on        Completed on
--------------------------------
Sunday                   8th
Monday                   7th
Tuesday                  6th
Wednesday                5th
Thursday                11th
Friday                  10th
Saturday                 9th

@author: durack1
"""

# import gc
import numpy as np
import MV2 as mv
import cmor
import datetime
import glob
import os
import pytz
import sys
import pdb
import cdms2 as cdm
import cdat_info as cdatInfo
import cdutil as cdu
from socket import gethostname
sys.path.insert(0, "/home/durack1/git/durolib/durolib")
from durolib import makeCalendar  # globalAttWrite, mkDirNoOSErr
sys.path.append("/home/durack1/git/input4MIPs-cmor-tables/src/")
from input4MIPsFuncs import createPubFiles, jsonWriteFile, washPerms
sys.path.insert(0, "pcmdiAmipBcs")
import pcmdiAmipBcsFx

# %% Kludge for json/encoder warning
# import warnings
##warnings.filterwarnings('ignore', category=DeprecationWarning)
# with warnings.catch_warnings():
#    warnings.filterwarnings('ignore', category=DeprecationWarning)
#    import vcs

# %% Set cdms format criteria
cdm.setNetcdfClassicFlag(0)  # 1 = amipbc files 437.2MB; amipobs files 437.3MB
# 1 = amipbc files sic 75.2MB, tos 239.2MB; amipobs files sic 29.3MB, tos 225.9MB
cdm.setNetcdfDeflateFlag(1)

# %% Set version info
activity_id = "input4MIPs"
# WILL REQUIRE UPDATING
contact = "pcmdi-cmip@lists.llnl.gov"
dataVerNum = "1.1.8"  # WILL REQUIRE UPDATING
dataVer = "PCMDI-AMIP-XX".replace("XX", dataVerNum.replace(".", "-"))
dataVerSht = "".join(["v", dataVerNum])
data_structure = "grid"
frequency = "mon"
# WILL REQUIRE UPDATING - point to GMD paper when available
further_info_url = "https://pcmdi.llnl.gov/mips/amip/"
institution_id = "PCMDI"
institution = "Program for Climate Model Diagnosis and Intercomparison, Lawrence Livermore National Laboratory, Livermore, CA 94550, USA"
last_year = "2021"  # WILL REQUIRE UPDATING
last_month = 12  # WILL REQUIRE UPDATING
comment = "Based on Hurrell SST/sea ice consistency criteria applied to merged HadISST (1870-01 to 1981-10) & NCEP-0I2 (1981-11 to 2021-06)"
comment = "".join(["Based on Hurrell SST/sea ice consistency criteria applied to ",
                   "merged HadISST (1870-01 to 1981-10) & NCEP-0I2 (1981-11 to ",
                   last_year, "-", "{:0>2}".format(last_month), ")"])
license_txt = " ".join(
    [
        "AMIP boundary condition data produced by PCMDI is licensed",
        "under a Creative Commons Attribution 4.0 International (CC BY 4.0;",
        "https://creativecommons.org/licenses/by/4.0/) License.",
        "The data producers and data providers make no warranty,",
        "either express or implied, including but not limited to,",
        "warranties of merchantability and fitness for a particular",
        "purpose. All liabilities arising from the supply of the",
        "information (including any liability arising in negligence)",
        "are excluded to the fullest extent permitted by law.",
    ]
)
mip_era = "CMIP6"  # "CMIP6Plus"
mip_specs = "AMIP CMIP5 CMIP6 CMIP6Plus"
project_id = "AMIP"
ref_obs = " ".join(
    [
        "Hurrell, J. W., J. J. Hack, D. Shea, J. M. Caron, and J. Rosinski",
        "(2008) A New Sea Surface Temperature and Sea Ice Boundary Dataset",
        "for the Community Atmosphere Model. J. Climate, 22 (19), pp",
        "5145-5153. doi: 10.1175/2008JCLI2292.1",
    ]
)
ref_bcs = " ".join(
    [
        "Taylor, K.E., D. Williamson and F. Zwiers, 2000: The sea surface",
        "temperature and sea ice concentration boundary conditions for",
        "AMIP II simulations. PCMDI Report 60, Program for Climate Model",
        "Diagnosis and Intercomparison, Lawrence Livermore National",
        "Laboratory, 25 pp. Available online: http://www-pcmdi.llnl.gov/publications/pdf/60.pdf",
    ]
)
source = "PCMDI-AMIP XX: Merged SST based on UK MetOffice HadISST and NCEP OI2".replace(
    "XX", dataVerNum
)
target_mip = "CMIP"
time_period = "".join(["187001-", last_year, "{:0>2}".format(last_month)])
destPath = "/p/user_pub/climate_work/durack1"
# For CMOR this is set in the CMOR/drive_input4MIPs*.json files
# destPath = '/work/durack1/Shared/150219_AMIPForcingData'  # USE FOR TESTING

# %% Get time/history info
utcNow = datetime.datetime.utcnow()
utcNow = utcNow.replace(tzinfo=pytz.utc)
timeFormat = utcNow.strftime("%d-%m-%Y %H:%M:%S %p")
localTz = pytz.timezone("America/Los_Angeles")
localNow = utcNow.astimezone(localTz)
cdatVerInfo = cdatInfo.version()
# Deal with quirky formats
if len(cdatVerInfo) > 2:
    cdatVerInfo = ".".join(["%s" % el for el in cdatInfo.version()])
else:
    cdatVerInfo = cdatInfo.version()[-1].strip("v")
    # Trim off the v
history = "".join(["File processed: ", timeFormat,
                  " UTC; San Francisco, CA, USA"])
host = "".join(
    [
        "Host: ",
        gethostname(),
        "; CDAT version: ",
        cdatVerInfo,
        "; Python version: ",
        sys.version.replace("\n", "; ").replace(") ;", ");"),
    ]
)
history = "".join([history, "; \n", host])
print(history)

# %% Set directories and input data
homePath = os.path.join(destPath, "Shared/150219_AMIPForcingData/")
# sanPath     = os.path.join(homePath,'_'.join(['360x180',dataVerSht,'san']))
sanPath = os.path.join(
    homePath, "".join(["SST_", dataVerNum.replace(".", "-")])
)
dataEnd = "202205"
# sanPath = os.path.join(homePath, "SST_1-2-0_old4")
# dataEnd = "202109"
# sanPath = os.path.join(homePath, "SST_1-2-0_old3")
# dataEnd = "202108"
# sanPath = os.path.join(homePath, "SST_1-2-0_old2")
# dataEnd = "202106"
# sanPath = os.path.join(homePath, "SST_1-2-0-1-1-6")
# dataEnd = "201903"
print("sanPath:", sanPath)

# %% Process each variable
varList = {}
varList["tos"] = {
    "varName": "SST",
    "fileVar": "SST",
    "ftype": "sst",
    "units": "degC",
    "outVar": "tosbcs",
}  # units = 'CELCIUS'
varList["siconc"] = {
    "varName": "SEAICE",
    "fileVar": "ICE",
    "ftype": "ice",
    "units": "%",
    "outVar": "siconcbcs",
}
# REPLACE READ WITH ORIGINAL DATA
for varId in ["siconc", "tos"]:
    # for varId in ['tos']:  # Testing
    varName = varList[varId]["varName"]
    fileVar = varList[varId]["fileVar"]
    ftype = varList[varId]["ftype"]
    units = varList[varId]["units"]
    outVar = varList[varId]["outVar"]
    inFile = "".join(
        ["MODEL.", fileVar, ".HAD187001-198110.OI198111-", dataEnd, ".nc"])
    fH = cdm.open(os.path.join(sanPath, inFile), "r")
    var = fH(varName)
    print("var.shape", var.shape)
    # print('var:', var)
    # print('var.getTime():', var.getTime())
    time = var.getTime().asComponentTime()
    lastYr = time[-1].year
    lastMn = time[-1].month
    print("lastYr:", lastYr, "lastMn:", lastMn)
    del time
    # Create calendar
    if lastMn == 6:
        endYr = str(int(lastYr))  # Half year/Same year
        time = makeCalendar("1870", endYr, monthEnd=(
            lastMn + 1), calendarStep="months")
    elif lastMn == 12:
        endYr = str(int(lastYr))  # Half year/Same year
        time = makeCalendar("1870", endYr + 1, monthEnd=1,
                            calendarStep="months")
    else:  # Case of full year; last_month = 12
        endYr = str(int(lastYr))  # Correct off by one, full year
        if lastMn == 12:
            # Dec (1) 2017 completion; June (6) 2016 completion
            time = makeCalendar(
                "1870", endYr, monthEnd=lastMn + 1, calendarStep="months"
            )
        else:
            print("Some calendar error, passing to pdb")
            # pdb.set_trace()
            # Dec (1) 2017 completion; June (6) 2016 completion
            time = makeCalendar(
                "1870", endYr, monthEnd=(lastMn + 1), calendarStep="months"
            )
    print("first:", time.asComponentTime()[0])
    print("last: ", time.asComponentTime()[-1])
    print("time.units:", time.units)
    print("time len:", len(time))
    print("varName:", varName)
    print("outVar:", outVar)
    print("var.units:", var.units)
    print("ftype:", ftype)

    # Reassign correct calendar/time axis to variable
    var.setAxis(0, time)

    # Get target grid
    targetGrid = var.getGrid()
    print("grid generated..")

    # Add mask as input
    fn_mask = "".join(
        [
            "/p/user_pub/work/input4MIPs/CMIP6/CMIP/PCMDI/",
            "PCMDI-AMIP-1-1-6/ocean/fx/sftof/gn/v20191121/",
            "sftof_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-6",
            "_gn.nc",
        ]
    )
    f = cdm.open(fn_mask)
    sftof = f("sftof")
    f.close()
    print("sftof read..")

    # run pcmdiAmipBcs/compile to refresh binaries (ensure conda env is consistent!)
    # create tos midpoint values
    print("Entering createMonthlyMidpoints function..")
    nyears = 10  # Buffer ~24-month climatology calculated over nyears
    varBcs = pcmdiAmipBcsFx.createMonthlyMidpoints(
        var, ftype, units, nyears, outVar
    )  # , grid=targetGrid, mask=sftof)
    print("Exiting createMonthlyMidpoints function..")

    # check input file and and output times
    print("inputFile:", dataEnd)
    print("".join([varId, ".shape:"]), var.shape)
    print("".join([varId, "bcs.shape:"]), varBcs.shape)
    time = var.getTime()
    print(time.asComponentTime()[-1])

    # Cleanup partial year data - always end on full or half years (12/6)
    endInd = np.mod(time.asComponentTime()[-1].month, 6)
    var = var[
        :-endInd,
    ]
    varBcs = varBcs[
        :-endInd,
    ]
    print("".join([varId, ".shape:"]), var.shape)
    print("".join([varId, "bcs.shape:"]), varBcs.shape)
    time = var.getTime()
    print(time.asComponentTime()[-1])
    if time.asComponentTime()[-1].month not in (6, 12):
        pdb.set_trace()

    # 1. Write siconc/tos and siconcBcs/tosBcs
    for varToCMOR in ["obs", "obsBcs"]:
        if varToCMOR == "obs":
            dataSetJson = "CMOR/drive_input4MIPs_obs.json"
            dataSetTime = "time"
            dHandle = "var"
            cmorVarId = varId
        else:
            dataSetJson = "CMOR/drive_input4MIPs_bcs.json"
            dataSetTime = "time1"
            dHandle = "varBcs"
            cmorVarId = "".join([varId, "bcs"])
        # Start CMORising
        cmor.setup(
            inpath="CMOR/input4MIPs-cmor-tables/Tables",
            netcdf_file_action=cmor.CMOR_REPLACE_4,
        )
        cmor.dataset_json(dataSetJson)
        d = eval(dHandle)
        lat = d.getLatitude()
        lon = d.getLongitude()
        time = d.getTime()
        # Force local file attribute as history
        cmor.set_cur_dataset_attribute("history", history)
        if "sic" in varId:
            table = "input4MIPs_SImon.json"
        else:
            table = "input4MIPs_Omon.json"
        cmor.load_table(table)
        axes = [
            {"table_entry": dataSetTime, "units": "days since 1870-01-01"},
            {
                "table_entry": "latitude",
                "units": "degrees_north",
                "coord_vals": lat[:],
                "cell_bounds": lat.getBounds(),
            },
            {
                "table_entry": "longitude",
                "units": "degrees_east",
                "coord_vals": lon[:],
                "cell_bounds": lon.getBounds(),
            },
        ]
        axis_ids = list()
        for axis in axes:
            axis_id = cmor.axis(**axis)
            axis_ids.append(axis_id)
        print("varName:", cmorVarId, "units:", units, "axis_ids:", axis_ids)
        varid = cmor.variable(cmorVarId, units, axis_ids)
        values = np.array(d[:], np.float32)
        # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
        cmor.set_deflate(varid, 1, 1, 1)
        cmor.write(varid, values,
                   time_vals=time[:], time_bnds=time.getBounds())
        f.close()
        cmor.close()

# 2. Create areacello and sftof and write

# areacello
areacello = cdu.area_weights(
    var[
        0,
    ]
)
# areacello.sum() = 1.0
earthSurfaceArea = 510.1
# million km2
earthSurfaceAreaM2 = earthSurfaceArea * 1e12
# m2
areacelloM2 = areacello * earthSurfaceAreaM2
areacelloM2.standard_name = "cell_area"
areacelloM2.long_name = "Ocean Grid-Cell Area"
areacelloM2.units = "m2"
areacelloM2.id = "areacello"
areacello = areacelloM2
del areacelloM2
# sftof
maskFile = os.path.join(
    destPath, "Shared/obs_data/WOD13/170425_WOD13_masks_1deg.nc")
fMask = cdm.open(maskFile)
landSea1deg = fMask("landsea")
# Fix longitude
a = landSea1deg[:, 0:180]
b = landSea1deg[:, 180:]
c = np.concatenate((b, a), axis=1)
del (a, b)
landSea1degTmp = cdm.createVariable(c, type="int16")
del c
landSea1degTmp.setAxis(0, areacello.getAxis(0))
# Impose identical axes to areacello
landSea1degTmp.setAxis(1, areacello.getAxis(1))
landSea1deg = landSea1degTmp
del landSea1degTmp
landSea1deg = mv.where(landSea1deg > 1.0, 0.0, landSea1deg)
# sea=0, land=1
landSea1deg = mv.where(landSea1deg == 1.0, 2.0, landSea1deg)
# Change land > 2.
landSea1deg = mv.where(landSea1deg == 0.0, 100.0, landSea1deg)
# Change sea > 100.
landSea1deg = mv.where(landSea1deg == 2.0, 0.0, landSea1deg)
# Change land > 0.
# Need to tweak some cells
sftof = cdm.createVariable(landSea1deg)
sftof.standard_name = "sea_area_fraction"
sftof.long_name = "Sea Area Fraction"
sftof.units = "%"
sftof.id = "sftof"

fxFiles = ["areacello", "sftof"]
for fxVar in fxFiles:
    cmor.setup(
        inpath="CMOR/input4MIPs-cmor-tables/Tables",
        netcdf_file_action=cmor.CMOR_REPLACE_4,
    )
    cmor.dataset_json("CMOR/drive_input4MIPs_obs.json")
    d = eval(fxVar)
    lat = d.getLatitude()
    lon = d.getLongitude()
    time = d.getTime()
    # Force local file attribute as history
    cmor.set_cur_dataset_attribute("history", history)
    # cmor.set_cur_dataset_attribute('frequency', 'fx')  # <-- test? no good
    table = "input4MIPs_Ofx.json"  # <-- doesn't overwrite source_id value
    cmor.load_table(table)
    axes = [
        {
            "table_entry": "latitude",
            "units": "degrees_north",
            "coord_vals": lat[:],
            "cell_bounds": lat.getBounds(),
        },
        {
            "table_entry": "longitude",
            "units": "degrees_east",
            "coord_vals": lon[:],
            "cell_bounds": lon.getBounds(),
        },
    ]
    axis_ids = list()
    for axis in axes:
        axis_id = cmor.axis(**axis)
        axis_ids.append(axis_id)
    varid = cmor.variable(fxVar, d.units, axis_ids, missing_value=1e20)
    values = np.array(d[:], np.float32)
    # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
    cmor.set_deflate(varid, 1, 1, 1)
    cmor.write(varid, values)
    cmor.close()

# %% Generate json files for publication step
jsonFilePaths, variableFilePaths = [[] for _ in range(2)]

# Get list of new files
dataVersion = datetime.datetime.now().strftime("v%Y%m%d")
files = glob.glob(
    os.path.join(
        destPath,
        activity_id,
        mip_era,
        target_mip,
        institution_id,
        dataVer,
        "*/*/*/gn",
        dataVersion,
        "*.nc",
    )
)
# OLD: 150219_AMIPForcingData/CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/mon/ocean/PCMDI-AMIP-1-1-3/tos/gn/v20171024/tos_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc
# NEW: /p/user_pub/work/input4MIPs/CMIP6/CMIP/PCMDI/PCMDI-AMIP-1-1-3/ocean/mon/tos/gn/v20171031/tos_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc
for filePath in files:
    print("filePath:", filePath)
    fH = cdm.open(filePath, "r")
    conventions = fH.Conventions
    activityId = fH.activity_id
    contact = fH.contact
    creationDate = fH.creation_date
    datasetCategory = fH.dataset_category
    sourceVersion = fH.source_version
    frequency = fH.frequency
    furtherInfoUrl = fH.further_info_url
    gridLabel = fH.grid_label
    institution = fH.institution
    institutionId = fH.institution_id
    mipEra = fH.mip_era
    nominalResolution = fH.nominal_resolution
    realm = fH.realm
    source = fH.source
    sourceId = fH.source_id
    targetMip = fH.target_mip
    targetMipJson = [fH.target_mip]
    title = fH.title
    trackingIdList = [fH.tracking_id]
    variableId = fH.variable_id
    outPath = filePath.replace(homePath, "")
    fileName = filePath.split("/")[-1]
    # Other vars
    activityId = "input4MIPs"
    deprecated = False
    jsonId = "CMIP-PaulDurack"
    jsonWriteFile(
        conventions,
        activityId,
        contact,
        creationDate,
        datasetCategory,
        sourceVersion,
        frequency,
        furtherInfoUrl,
        gridLabel,
        institution,
        institutionId,
        mipEra,
        nominalResolution,
        realm,
        source,
        sourceId,
        targetMip,
        targetMipJson,
        title,
        variableId,
        filePath,
        trackingIdList,
        deprecated,
        dataVersion,
        destPath,
        jsonId,
    )

    # Save json file list for publication
    destFilePath = os.path.join(
        activityId,
        mipEra,
        targetMip,
        institutionId,
        sourceId,
        realm,
        frequency,
        variableId,
        gridLabel,
        dataVersion,
        fileName,
    )
    variableFilePaths.append(os.path.join(
        destPath, destFilePath.replace(fileName, "")))
    jsonFilePath = os.path.join(
        destPath,
        activityId,
        mipEra,
        targetMip,
        institutionId,
        "".join(
            [
                "_".join(
                    [
                        institutionId,
                        frequency,
                        sourceId,
                        variableId,
                        gridLabel,
                        dataVersion,
                    ]
                ),
                ".json",
            ]
        ),
    )
    jsonFilePaths.append(jsonFilePath)

# Clean up permissions - hardcoded to /p/
# washPerms(
#     destPath,
#     activityId,
#     mipEra,
#     targetMip,
#     institutionId,
#     sourceId,
#     realm,
#     frequency,
#     gridLabel,
#     dataVersion,
# )
# Create output files for publication
createPubFiles(destPath, jsonId, jsonFilePaths, variableFilePaths)
