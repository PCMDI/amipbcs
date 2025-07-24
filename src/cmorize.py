#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 15:39:03 2025

Paul J. Durack 21st Jul 2025

This script cmorizes nc files

PJD 23 Jul 25 - updates to map xr functions to replace durolib
TODO: remap all cdms2 dependencies to xcdat
TODO: remap makeCalendar to https://docs.xarray.dev/en/latest/generated/xarray.date_range.html
"""

# %% imports
##%%time
import cftime
import cmor
import datetime

# import json
import numpy as np
import os
import pdb
import socket
import sys
import xcdat as xc

# import xarray as xr
sys.path.insert(0, "pcmdiAmipBcs")
import pcmdiAmipBcsFx

# %% set data version info
##%%time
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
destPath = "/p/user_pub/climate_work/durack1"
# For CMOR this is set in the CMOR/drive_input4MIPs*.json files
# destPath = '/p/user_pub/climate_work/durack1/Shared/150219_AMIPForcingData'  # USE FOR TESTING

# %% get time/history/host info
##%%time
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
##%%time
homePath = os.path.join(destPath, "Shared/150219_AMIPForcingData/")
homePath = "./"
sanPath = os.path.join(homePath, "".join(["SST_", dataVerNum.replace(".", "-")]))
dataEnd = "202301"
print("sanPath:", sanPath)
print("os.getcwd():", os.getcwd())

# %% preload data iterating over each variable, fix calendar, diddle and pass to CMOR
##%%time
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
    # fH = cdm.open(os.path.join(sanPath, inFile), "r")
    fH = xc.open_dataset(os.path.join(sanPath, inFile))
    fH = fH.bounds.add_time_bounds(method="midpoint")  # add time bounds
    xrVar = ".".join(["fH", varName])
    print("xrVar:", xrVar)
    var = eval(xrVar)

    """
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
        time = makeCalendar("1870", endYr, monthEnd=(lastMn + 1), calendarStep="months")
    elif lastMn == 12:
        endYr = str(int(lastYr))  # Half year/Same year
        time = makeCalendar("1870", endYr + 1, monthEnd=1, calendarStep="months")
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
    #targetGrid = var.getGrid()
    #print("grid generated..")

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
"""

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
    # time = var.getTime()
    # print(time.asComponentTime()[-1])
    print(fH.time)

    # Cleanup partial year data - always end on full or half years (12/6)
    ###endInd = np.mod(time.asComponentTime()[-1].month, 6)
    endInd = np.mod(fH.time.dt.month[-1].data, 6)
    var = var[:-endInd,]
    varBcs = varBcs[:-endInd,]
    print("".join([varId, ".shape:"]), var.shape)
    print("".join([varId, "bcs.shape:"]), varBcs.shape)
    ###time = var.getTime()
    ###print(time.asComponentTime()[-1])
    print(var.time[-1])
    ###if time.asComponentTime()[-1].month not in (6, 12):
    if var.time.dt.month[-1] not in (6, 12):
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
            netcdf_file_action=cmor.CMOR_REPLACE_4,
        )
        # inpath="CMOR/input4MIPs-cmor-tables/Tables",
        cmor.dataset_json(dataSetJson)
        ###d = eval(dHandle)
        ###lat = d.getLatitude()
        lat = var.cf["latitude"]
        ###lon = d.getLongitude()
        lon = var.cf["longitude"]
        ###time = d.getTime()
        time = var.cf["time"]
        # Force local file attribute as history
        cmor.set_cur_dataset_attribute("history", history)
        if "sic" in varId:
            table = "input4MIPs_SImon.json"
        else:
            table = "input4MIPs_Omon.json"

        # Fudge table files to force project=CMIP6Plus
        tablePath = "Tables"
        # with open(os.path.join(tablePath, table)) as fh:
        #    tmp = json.load(fh)
        # tmp["Header"]["mip_era"] = "CMIP6Plus"
        # oH = open(os.path.join(tablePath, table)
        # json.dump(
        #    tmp, oH, ensure_ascii=True, sort_keys=True, indent=4, separators=(",", ":")
        # )
        # oH.close()
        # cmor.load_table("tmp.json")
        tablePath = os.path.join(tablePath, table)
        print("tablePath:", tablePath)
        cmor.load_table(tablePath)
        # os.remove("tmp.json")

        ###pdb.set_trace()

        axes = [
            {"table_entry": dataSetTime, "units": "days since 1870-01-01 0:0:0.0"},
            {
                "table_entry": "latitude",
                "units": "degrees_north",
                "coord_vals": lat.data,
                "cell_bounds": fH["lat_bnds"].data,
            },
            {
                "table_entry": "longitude",
                "units": "degrees_east",
                "coord_vals": lon.data,
                "cell_bounds": fH["lon_bnds"].data,
            },
        ]
        axis_ids = list()
        for axis in axes:
            axis_id = cmor.axis(**axis)
            axis_ids.append(axis_id)
        print("varName:", cmorVarId, "units:", units, "axis_ids:", axis_ids)
        varid = cmor.variable(cmorVarId, units, axis_ids)
        values = np.array(var[:], np.float32)
        # shuffle=1,deflate=1,deflate_level=1 ; CMOR 3.0.6+
        cmor.set_deflate(varid, 1, 1, 1)
        cmor.write(
            varid,
            values,
            time_vals=cftime.date2num(time, "days since 1870-1-1 0:0:0.0"),
            time_bnds=cftime.date2num(fH["time_bnds"], "days since 1870-1-1 0:0:0.0"),
        )
        cmor.close()

# %%
