#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 09:46:03 2025

Paul J. Durack 15st Aug 2025

This script cmorizes nc files

PJD 15 Aug 25 - copied from cmorize.py and updated input
"""

# %% imports
import cftime as cft
import cmor
import datetime
import json
import numpy as np
import os
import pdb
import subprocess
import sys
import xcdat as xc
import xarray as xr

sys.path.insert(0, "pcmdiAmipBcs")
import pcmdiAmipBcsFx


# %% set data version info
"""
activity_id = "input4MIPs"
contact = "pcmdi-cmip@llnl.gov"
dataUpdateNotes = "".join(
    [
        "v1.1.9 and v1.1.10 differences: this ",
        "update changes a single month (Dec-22) ",
        "erroneous sea ice concentration ",
        "(siconc). Due to the tapering affect ",
        "of the 'diddling' method, some very ",
        "small changes (<1 percent) can be seen ",
        "starting in August 2022 in diddled ",
        "fields (siconcbcs). For v1.1.10 a ",
        "climatology-anomaly infill was ",
        "undertaken, replacing the Dec-22 ",
        "problem values. For more details, see ",
        "https://nbviewer.org/github/durack1/",
        "notebooks/blob/main/jlnbs/PCMDI-AMIP-",
        "queryOISST2-0Data.ipynb; There are no ",
        "changes to either the SST (tos) or ",
        "diddled SST (tosbcs) fields; NOAA OISST ",
        "v2.0 data was deprecated in February ",
        "2023, and no further PCMDI-AMIP-1-x-y ",
        "updates will be produced. Ongoing ",
        "discussions focused on a v2.0 product ",
        "continue, see https://github.com/PCMDI/",
        "amipbcs/issues/6.",
    ]
)
dataRepo = "https://github.com/PCMDI/amipbcs"
# WILL REQUIRE UPDATING
dataVerNum = "1.1.10"
dataVer = "PCMDI-AMIP-XX".replace("XX", dataVerNum.replace(".", "-"))
dataVerSht = "".join(["v", dataVerNum])
data_structure = "grid"
doi = "https://doi.org/10.25981/ESGF.input4MIPs.CMIP7/2575015"
frequency = "mon"
# WILL REQUIRE UPDATING - point to GMD paper when available
further_info_url = "https://pcmdi.llnl.gov/mips/amip/"
institution_id = "PCMDI"
institution = " ".join(
    [
        "Program for Climate Model Diagnosis and Intercomparison,",
        "Lawrence Livermore National Laboratory, Livermore, CA",
        "94550, USA (ROR: https://ror.org/02k3nmd98)",
    ]
)
last_year = "2022"  # WILL REQUIRE UPDATING
last_month = 12  # WILL REQUIRE UPDATING
comment = "".join(
    [
        "Based on Hurrell SST/sea ice consistency criteria applied to ",
        "merged HadISST v1.0 (1870-01 to 1981-10) & NCEP-0ISST v2.0 ",
        "(1981-11 to ",
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
        "Laboratory, 25 pp. Available online:",
        "https://pcmdi.llnl.gov/report/pdf/60.pdf",
    ]
)
source = " ".join(
    [
        "PCMDI-AMIP XX: Merged SST based on UK MetOffice HadISST v1.0",
        "and NCEP OISST v2.0",
    ]
)
source.replace("XX", dataVerNum)
target_mip = "CMIP"
time_period = "".join(["187001-", last_year, "{:0>2}".format(last_month)])
"""

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
dataPaths = {
    "PCMDI-AMIP-ERSST5-1-0": {
        "filePath": "newProto/NOAA-ERSST-V5/MODEL.SST.HAD187001-198110.OI198111-202301.NOAA_ERSST_V5.nc",
        "sourceId": "PCMDI-AMIP-ERSST5-1-0",
    },
    "PCMDI-AMIP-Had1p1-1-0": {
        "filePath": "newProto/HadISST-1p1/MODEL.SST.HAD187001-198110.OI198111-202301.HadISST-1.1.nc",
        "sourceId": "PCMDI-AMIP-Had1p1-1-0",
    },
    "PCMDI-AMIP-OI2p1-1-0": {
        "filePath": "newProto/NOAA-OISST-v2.1/MODEL.SST.HAD187001-198110.OI198111-202301.NOAA-OISST-v2.1.nc",
        "sourceId": "PCMDI-AMIP-OI2p1-1-0",
    },
}
for count, dataset in enumerate(dataPaths.keys()):
    # destPath = "/p/user_pub/climate_work/durack1"  ## LLNL/detect
    # destPath = "/pscratch/sd/d/durack1/"
    destPath = "."

    # homePath = os.path.join(destPath, "Shared/150219_AMIPForcingData/")
    # sanPath = os.path.join(homePath, "".join(["SST_", dataVerNum.replace(".", "-")]))
    print(dataPaths[dataset]["filePath"])
    sanPath = os.path.join(destPath, dataPaths[dataset]["filePath"])
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
    # inFile = "".join(["MODEL.", fileVar, ".HAD187001-198110.OI198111-", dataEnd, ".nc"])
    fH = xc.open_dataset(sanPath)
    # , decode_times=False)
    # add CF-required attributes back in
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
