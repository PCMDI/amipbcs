#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 15:39:03 2025

Paul J. Durack 21st Jul 2025

This script cmorizes nc files

PJD 23 Jul 25 - updates to map xr functions to replace durolib
PJD 24 Jul 25 - remapped all dependencies to xcdat/array (remove
                cdms2 dependence)
PJD 28 Jul 25 - updated CMOR time_units to remove HH:MM:SS.x does
                this fix the ~6 hrs temporal offset
PJD 28 Jul 25 - update to xc.open_dataset($file, decode_times=False)
PJD 29 Jul 25 - further updates further cleaning up redundant code
                and correctly assigns obs/bcs vars for CMOR writes
                indexing variables across the obs and bcs variants
PJD 29 Jul 25 - updated to replace time_bnds with generated calendar
                xarray.date_range
PJD  4 Aug 25 - updated for perlmutter and paths
"""

# %% imports
import cftime as cft
import cmor
import datetime
import numpy as np
import os
import socket
import sys
import xcdat as xc
import xarray as xr

sys.path.insert(0, "pcmdiAmipBcs")
import pcmdiAmipBcsFx


# %% set data version info
activity_id = "input4MIPs"
contact = "pcmdi-cmip@llnl.gov"
dataVerNum = "1.1.10"  # WILL REQUIRE UPDATING
dataVer = "PCMDI-AMIP-XX".replace("XX", dataVerNum.replace(".", "-"))
dataVerSht = "".join(["v", dataVerNum])
data_structure = "grid"
frequency = "mon"
further_info_url = "https://pcmdi.llnl.gov/mips/amip/"  # WILL REQUIRE UPDATING - point to GMD paper when available
institution_id = "PCMDI"
institution = " ".join(
    [
        "Program for Climate Model Diagnosis and Intercomparison,",
        "Lawrence Livermore National Laboratory, Livermore, CA 94550, USA",
    ]
)
last_year = "2022"  # WILL REQUIRE UPDATING
last_month = 12  # WILL REQUIRE UPDATING
comment = "".join(
    [
        "Based on Hurrell SST/sea ice consistency criteria applied to ",
        "merged HadISST (1870-01 to 1981-10) & NCEP-0I2 (1981-11 to ",
        last_year,
        "-",
        "{:0>2}".format(last_month),
        ")",
    ]
)
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
mip_era = "CMIP7"  # "CMIP6Plus"  # "CMIP6"
mip_specs = "AMIP CMIP5 CMIP6 CMIP6Plus CMIP7"
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
        "Laboratory, 25 pp. Available online: https://pcmdi.llnl.gov/report/pdf/60.pdf",
    ]
)
source = "PCMDI-AMIP XX: Merged SST based on UK MetOffice HadISST and NCEP OI2".replace(
    "XX", dataVerNum
)
target_mip = "CMIP"
time_period = "".join(["187001-", last_year, "{:0>2}".format(last_month)])

# %% get time/history/host info
utcNow = datetime.datetime.now(datetime.timezone.utc)
timeFormat = utcNow.strftime("%d-%m-%Y %H:%M:%S %p")
xcVersion = xc.__version__
history = " ".join(["File processed:", timeFormat, "UTC; San Francisco, CA, USA"])
host = "".join(
    [
        "Host: ",
        socket.gethostname(),
        "; xCDAT version: ",
        xcVersion,
        "; Python version: ",
        sys.version.split(" |")[0],
        ";",
    ]
)
history = "".join([history, "; \n", host])
print(history)

# %% Set directories and input data
# destPath = "/p/user_pub/climate_work/durack1"  ## LLNL/detect
destPath = "/pscratch/sd/d/durack1/"

homePath = os.path.join(destPath, "Shared/150219_AMIPForcingData/")
sanPath = os.path.join(homePath, "".join(["SST_", dataVerNum.replace(".", "-")]))
dataEnd = "202301"
print("sanPath:", sanPath)
print("os.getcwd():", os.getcwd())

# %% create replacement calendar/time_bnds
newCal = xr.date_range(
    start="1870", end="2024", freq="MS", calendar="gregorian", use_cftime=True
)
newCal187001to202301 = newCal[:-12]  # trim to end of 2023-01
time_bnds = np.stack((newCal187001to202301[:-1], newCal187001to202301[1:]), axis=1)

# %% preload data, iterate over variables, diddle, and pass to CMOR
varList = {}
varList["tos"] = {
    "varName": "SST",
    "fileVar": "SST",
    "ftype": "sst",
    "units": "degC",
    "outVar": "tosbcs",
}
varList["siconc"] = {
    "varName": "SEAICE",
    "fileVar": "ICE",
    "ftype": "ice",
    "units": "%",
    "outVar": "siconcbcs",
}

# %% REPLACE READ WITH ORIGINAL DATA
for varId in ["siconc", "tos"]:
    # for varId in ['tos']:  # Testing
    varName = varList[varId]["varName"]
    fileVar = varList[varId]["fileVar"]
    ftype = varList[varId]["ftype"]
    units = varList[varId]["units"]
    outVar = varList[varId]["outVar"]
    inFile = "".join(["MODEL.", fileVar, ".HAD187001-198110.OI198111-", dataEnd, ".nc"])
    fH = xc.open_dataset(os.path.join(sanPath, inFile))
    # , decode_times=False)
    fH = fH.bounds.add_time_bounds(method="midpoint")  # add time bounds
    xrVar = ".".join(["fH", varName])
    print("xrVar:", xrVar)
    var = eval(xrVar)

    # run pcmdiAmipBcs/compile - refresh binaries (ensure env consistent!)
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
    print(fH.time)

    # Cleanup partial year data - always end on full or half years (12/6)
    endInd = np.mod(fH.time.dt.month[-1].data, 6)
    var = var[:-endInd,]
    varBcs = varBcs[:-endInd,]
    print("".join([varId, ".shape:"]), var.shape)
    print("".join([varId, "bcs.shape:"]), varBcs.shape)
    print(var.time[-1])
    if var.time.dt.month[-1] not in (6, 12):
        print("Catch case of bad data..")
        pdb.set_trace()

    # %% CMORize
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
            inpath="Tables",
            set_verbosity=cmor.CMOR_NORMAL,
            netcdf_file_action=cmor.CMOR_REPLACE_4,
        )
        cmor.dataset_json(dataSetJson)

        # Force local file attribute as history
        cmor.set_cur_dataset_attribute("history", history)
        if "sic" in varId:
            table = "input4MIPs_SImon.json"
        else:
            table = "input4MIPs_Omon.json"

        # Load relevant table file
        tablePath = "Tables"
        tablePath = os.path.join(tablePath, table)
        print("tablePath:", tablePath)
        cmor.load_table(tablePath)

        axes = [
            {"table_entry": dataSetTime, "units": "days since 1870-01-01"},
            {
                "table_entry": "latitude",
                "units": "degrees_north",
                "coord_vals": fH.lat.data,
                "cell_bounds": fH["lat_bnds"].data,
            },
            {
                "table_entry": "longitude",
                "units": "degrees_east",
                "coord_vals": fH.lon.data,
                "cell_bounds": fH["lon_bnds"].data,
            },
        ]
        axis_ids = list()
        for axis in axes:
            axis_id = cmor.axis(**axis)
            axis_ids.append(axis_id)
        print("varName:", cmorVarId, "units:", units, "axis_ids:", axis_ids)
        varid = cmor.variable(cmorVarId, units, axis_ids)
        values = np.array(eval(dHandle), np.float32)  # output either obs/bcs
        # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
        cmor.set_deflate(varid, 1, 1, 1)

        cmor.write(
            varid,
            values,
            time_vals=cft.date2num(var.cf["time"], "days since 1870-1-1"),
            time_bnds=cft.date2num(time_bnds, "days since 1870-1-1"),
        )
        del values  # explicitly purge so a new copy is generated

        cmor.close()

# 2. Create areacello and sftof and write

# areacello
areacello = fH.spatial.get_weights(axis="Y")
areacello, _ = xr.broadcast(areacello, var[0])
areacello = areacello / 720.0  # need to remap back to cdutil values
# areacello.sum() = 1.0  # validated
earthSurfaceArea = 510.1
# million km2
earthSurfaceAreaM2 = earthSurfaceArea * 1e12
# m2
areacelloM2 = areacello * earthSurfaceAreaM2
areacello = areacelloM2
del areacelloM2

# sftof
maskPath = os.path.join(destPath, "Shared/obs_data/WOD13")
maskFile = os.path.join(maskPath, "170425_WOD13_masks_1deg.nc")
fM = xc.open_dataset(maskFile)
landSea1deg = fM["landsea"]
# Fix longitude
aMat = landSea1deg[:, 0:180]
bMat = landSea1deg[:, 180:]
cMat = np.concatenate((bMat, aMat), axis=1)
del (aMat, bMat)

landSea1degTmp = xr.DataArray(
    cMat.astype(np.int16),
    dims=("lat", "lon"),
    coords={k: var.coords[k] for k in ["lat", "lon"]},
)
del cMat

landSea1deg = landSea1degTmp
del landSea1degTmp
landSea1deg = np.where(landSea1deg > 1.0, 0.0, landSea1deg)
# sea=0, land=1
landSea1deg = np.where(landSea1deg == 1.0, 2.0, landSea1deg)
# Change land > 2.
landSea1deg = np.where(landSea1deg == 0.0, 100.0, landSea1deg)
# Change sea > 100.
landSea1deg = np.where(landSea1deg == 2.0, 0.0, landSea1deg)
# Change land > 0.
sftof = landSea1deg
del landSea1deg

fxFiles = ["areacello", "sftof"]
for fxVar in fxFiles:
    if fxVar == "areacello":
        units = "m2"
    elif fxVar == "sftof":
        units = "%"

    # load relevant variable
    d = eval(fxVar)

    # setup cmor
    cmor.setup(
        inpath="Tables",
        set_verbosity=cmor.CMOR_NORMAL,
        netcdf_file_action=cmor.CMOR_REPLACE_4,
    )
    cmor.dataset_json("CMOR/drive_input4MIPs_obs.json")

    # Force local file attribute as history
    cmor.set_cur_dataset_attribute("history", history)

    # Load relevant table file
    tablePath = "Tables"
    table = "input4MIPs_Ofx.json"  # <-- doesn't overwrite source_id value
    cmor.load_table(os.path.join(tablePath, table))

    axes = [
        {
            "table_entry": "latitude",
            "units": "degrees_north",
            "coord_vals": fH.lat.data,
            "cell_bounds": fH["lat_bnds"].data,
        },
        {
            "table_entry": "longitude",
            "units": "degrees_east",
            "coord_vals": fH.lon.data,
            "cell_bounds": fH["lon_bnds"].data,
        },
    ]
    axis_ids = list()
    for axis in axes:
        axis_id = cmor.axis(**axis)
        axis_ids.append(axis_id)
    varid = cmor.variable(fxVar, units, axis_ids, missing_value=1e20)
    values = np.array(d[:], np.float32)
    # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
    cmor.set_deflate(varid, 1, 1, 1)
    cmor.write(varid, values)
    cmor.close()
