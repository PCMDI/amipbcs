#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 09:46:03 2025

Paul J. Durack 15th Aug 2025

This script cmorizes CMIP7Plus nc files

PJD 15 Aug 25 - copied from cmorize.py and updated input
PJD 27 Aug 25 - updated to run on NERSC with m4581/zelinka1 paths
"""

# %% imports
import datetime
import json
import os
import pdb
import subprocess
import sys
import cftime as cft
import cmor
import numpy as np
import xcdat as xc
import xarray as xr

sys.path.insert(0, "pcmdiAmipBcs")
import pcmdiAmipBcsFx

# %% get time/history/host info
utcNow = datetime.datetime.now(datetime.timezone.utc)
timeFormat = utcNow.strftime("%d-%m-%Y %H:%M:%S %p")
xcVersion = xc.__version__
history = " ".join(["File processed:", timeFormat, "UTC; San Francisco, CA, USA"])
host = "".join(
    [
        "Host: ",
        subprocess.check_output(["hostname", "-f"], text=True).strip(),
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
srcPath = "/global/cfs/projectdirs/m4581/zelinka1"
dataPaths = {
    "PCMDI-AMIP-ERSST5-1-0": {
        "filePath": "NOAA_ERSST_V5/MODEL.SST.HAD187001-198110.OI198111-202301.NOAA_ERSST_V5.nc",
        "sourceId": "PCMDI-AMIP-ERSST5-1-0",
    },
    "PCMDI-AMIP-Had1p1-1-0": {
        "filePath": "HadISST-1.1/MODEL.SST.HAD187001-198110.OI198111-202301.HadISST-1.1.nc",
        "sourceId": "PCMDI-AMIP-Had1p1-1-0",
    },
    "PCMDI-AMIP-OI2p1-1-0": {
        "filePath": "NOAA-OISST-v2.1/MODEL.SST.HAD187001-198110.OI198111-202301.NOAA-OISST-v2.1.nc",
        "sourceId": "PCMDI-AMIP-OI2p1-1-0",
    },
}

# destPath = "." # test
destPath = "/global/cfs/projectdirs/m4581/durack1"  # NERSC

for count, dataset in enumerate(dataPaths.keys()):
    # set file paths
    # homePath = os.path.join(destPath, "Shared/150219_AMIPForcingData/")
    # sanPath = os.path.join(homePath, "".join(["SST_", dataVerNum.replace(".", "-")]))
    print(dataPaths[dataset]["filePath"])
    sanPath = os.path.join(srcPath, dataPaths[dataset]["filePath"])
    dataEnd = "202301"
    print("sanPath:", sanPath)
    print("os.getcwd():", os.getcwd())

    # %% create replacement calendar/time_bnds
    newCal = xr.date_range(
        start="1870", end="2024", freq="MS", calendar="gregorian", use_cftime=True
    )
    newCal187001to202301 = newCal[:-12]  # trim to end of 2023-01
    time_bnds = np.stack((newCal187001to202301[:-1], newCal187001to202301[1:]), axis=1)

    # %% REPLACE READ WITH ORIGINAL DATA

    varId = "tos"
    varName = "SST"
    fileVar = "SST"
    ftype = "sst"
    units = "degC"
    outVar = "tosbcs"
    fH = xc.open_dataset(sanPath)
    # , decode_times=False)
    # add CF-required attributes back in - required by pcmdiAmipBcs
    fH["lat"].attrs["units"] = "degrees_north"
    fH["lat"].attrs["long_name"] = "latitude"
    fH["lon"].attrs["units"] = "degrees_east"
    fH["lon"].attrs["long_name"] = "longitude"
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

    # Write tos and tosBcs
    for varToCMOR in ["obs", "obsBcs"]:
        match varToCMOR:
            case "obs":
                dataSetTime = "time"
                dHandle = "var"
                product = "observations"
                cmorVarId = varId
            case "obsBcs":
                dataSetTime = "time1"
                dHandle = "varBcs"
                product = "derived"
                cmorVarId = "".join([varId, "bcs"])

        # create user_input.json
        tmp = {}
        tmp["activity_id"] = "input4MIPs"
        tmp["source_id"] = dataset
        tmp["product"] = product
        tmp["outpath"] = destPath
        tmp["_history_template"] = (
            "%s; CMOR rewrote data to be consistent with <activity_id>, CMIP6, CMIP6Plus and <Conventions> standards"
        )
        tmp["output_path_template"] = (
            "<activity_id><mip_era><target_mip><institution_id><source_id><realm><frequency><variable_id><grid_label><version>"
        )
        tmp["output_file_template"] = (
            "<variable_id><activity_id><dataset_category><target_mip><source_id><grid_label>"
        )
        tmp["tracking_prefix"] = "hdl:21.14100"
        tmp["_controlled_vocabulary_file"] = "input4MIPs_CV.json"
        tmp["_AXIS_ENTRY_FILE"] = "input4MIPs_coordinate.json"
        tmp["_FORMULA_VAR_FILE"] = "input4MIPs_formula_terms.json"
        # cleanup
        if os.path.exists("tmp.json"):
            os.remove("tmp.json")
        with open("tmp.json", "w") as f:
            json.dump(
                tmp,
                f,
                ensure_ascii=True,
                sort_keys=True,
                indent=4,
                separators=(",", ":"),
            )

        # Start CMORising
        cmor.setup(
            inpath="Tables",
            set_verbosity=cmor.CMOR_NORMAL,
            netcdf_file_action=cmor.CMOR_REPLACE_4,
        )
        cmor.dataset_json("tmp.json")
        os.remove("tmp.json")

        # Force local file attribute as history
        cmor.set_cur_dataset_attribute("history", history)

        # Toggle data and appropriate table
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
