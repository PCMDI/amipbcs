#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 13:21:49 2021

PJD 20 Sep 2021 - Started
PJD 30 Sep 2021 - Updated for v20210930 data
PJD 27 Oct 2021 - Updated to generate diff as % maps
PJD  2 Nov 2021 - Updated following tweaks in https://github.com/PCMDI/amipbcs/issues/23#issuecomment-958164331
PJD  3 May 2023 - Updates for the v1.1.9 data
PJD  3 May 2023 - Updated for cdms2 -> xcdat
PJD  3 May 2023 - Added plotter function
PJD  4 May 2023 - Hitting issue with 2002-11 timestep and xarray DataArray plotting
PJD  9 May 2023 - Add transform_first=True to contourf call
PJD  9 May 2023 - Added +1 for last year, off by one PCMDI-AMIP-1-1-8 finishes in 2021-12
PJD  9 May 2023 - Added ffmpeg call - installed ffmpeg-python
PJD 10 May 2023 - Added statsStr to diff plot; updated diff scale <2%; corrected denom da1 vs s1 ref
PJD 10 May 2023 - Add statsStr; update contour levels to target data
PJD 12 May 2023 - Updated for latest v1.1.9 data run
PJD 18 May 2023 - Compare v1.1.9 versions (released 20230512, and new mamba env 20230518)
PJD 24 Jul 2025 - Updating for mac-local v1.1.10 vs v1.1.9
PJD 29 Jul 2025 - Updating for mac-local v1.1.10 v20250724 -> 20250729
PJD  4 Aug 2025 - Updating for perlmutter v1.1.10
PJD  4 Aug 2025 - Updated output *mp4 to take verId as filename arg

@author: durack1
"""

# %% imports
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import ffmpeg
import glob
import numpy as np
import os
import shutil
from xcdat import open_dataset

# %% function defs


def statsStr(da):
    """
    Create stats string to add to plot text box

    """
    fmtStr = "{:7.4f}"
    statsStr = " ".join(
        [
            "min:",
            fmtStr.format(da.min()),
            "\n1pc:",
            fmtStr.format(np.percentile(da, 1)),
            "\nmean:",
            fmtStr.format(da.mean()),
            "\n99pc:",
            fmtStr.format(np.percentile(da, 99)),
            "\nmax:",
            fmtStr.format(da.max()),
        ]
    )

    return statsStr


def plotter(
    da1,
    da2,
    da1Str,
    da2Str,
    lev1,
    lev2,
    cmap,
    timeStr,
    titleString,
    varColStr,
    path,
    var,
    fileName,
):
    """
    Generate generic plotting function
    """
    # Open canvas
    fig = plt.figure(figsize=(10, 15))
    plt.axis("off")
    plt.ioff()  # turn off interactive plots - background mode
    plt.title(titleString)

    # prepare lon, lat
    lon = np.tile(da1.lon.data, (180, 1))
    lat = np.tile(da1.lat.data, (360, 1)).transpose()

    # Start subplots
    ax1 = fig.add_subplot(
        3,
        1,
        1,
        projection=ccrs.Robinson(
            central_longitude=centralLon,
            globe=None,
            false_easting=None,
            false_northing=None,
        ),
    )
    # print("type(da1.lon):", type(da1.lon))
    # print("type(da1.lat):", type(da1.lat))
    # print("type(da1):", type(da1))
    # print("type(da2):", type(da2))
    # print("type(lev1):", type(lev1))
    # print("type(cmap):", type(cmap))
    # pdb.set_trace()
    # x1 = np.squeeze(np.array(da1.data))
    # la1 = np.array(da1.lat.data)
    # lo1 = np.array(da1.lon.data)
    # Failing line - only when an xarray DataArray is sent
    # cs1 = ax1.contourf(da1.lon, da1.lat, da1[0,],
    # cs1 = ax1.contourf(da1.lon.data, da1.lat.data, da1.squeeze().data,
    cs1 = ax1.contourf(
        lon,
        lat,
        da1[0,],
        lev1,  # 20
        transform=ccrs.PlateCarree(),
        transform_first=True,
        cmap=cmap,
    )
    tx1 = plt.text(
        labX,
        labY,
        da1Str,
        fontsize=fntsz,
        horizontalalignment="center",
        transform=ccrs.Geodetic(),
    )
    # create da1 dob variables
    tx2 = plt.text(
        labX,
        labY - 20,
        statsStr(da1),
        fontsize=fntsz,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ccrs.Geodetic(),
    )
    ax2 = fig.add_subplot(
        3,
        1,
        2,
        projection=ccrs.Robinson(
            central_longitude=centralLon,
            globe=None,
            false_easting=None,
            false_northing=None,
        ),
    )
    cs2 = ax2.contourf(
        lon,
        lat,
        da2.squeeze().data,
        # cs2 = ax2.contourf(lo1, la1, x2,
        lev1,
        transform=ccrs.PlateCarree(),
        transform_first=True,
        cmap=cmap,
    )
    tx3 = plt.text(
        labX,
        labY,
        da2Str,
        fontsize=fntsz,
        horizontalalignment="center",
        transform=ccrs.Geodetic(),
    )
    # create da2 dob variables
    tx4 = plt.text(
        labX,
        labY - 20,
        statsStr(da2),
        fontsize=fntsz,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ccrs.Geodetic(),
    )
    ax3 = fig.add_subplot(
        3,
        1,
        3,
        projection=ccrs.Robinson(
            central_longitude=centralLon,
            globe=None,
            false_easting=None,
            false_northing=None,
        ),
    )
    # Generate % change
    diff = da1[0,] - da2[0,]
    inds = np.nonzero(diff.data)
    diffnew = np.ma.zeros(diff.shape)
    denom = (np.abs(da1[0,]) + np.abs(da2[0,])) / 2
    np.squeeze(denom).shape
    diffnew[inds] = diff.data[inds] / denom.data[inds]

    cs3 = ax3.contourf(
        lon,
        lat,
        diffnew,
        lev2,
        transform=ccrs.PlateCarree(),
        transform_first=True,
        cmap=cmap,
    )
    tx5 = plt.text(
        labX,
        labY,
        " ".join([da1Str, "-", da2Str]),
        fontsize=fntsz,
        horizontalalignment="center",
        transform=ccrs.Geodetic(),
    )
    # create diff dob variables
    tx6 = plt.text(
        labX,
        labY - 20,
        statsStr(diffnew),
        fontsize=fntsz,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ccrs.Geodetic(),
    )

    # make the map global rather than have it zoom in to the extents of
    # any plotted data ax.set_global()
    # ax1.stock_img()
    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()
    # ax.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
    # ax.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.Geodetic())

    # https://matplotlib.org/stable/gallery/axes_grid1/demo_colorbar_with_inset_locator.html
    axin1 = inset_axes(
        ax1,
        width="5%",  # width = 5% of parent_bbox width
        height="50%",  # height : 50%
        loc="lower left",
        bbox_to_anchor=(1.05, -1.17, 1, 4.3),
        bbox_transform=ax1.transAxes,
        borderpad=0,
    )

    axin3 = inset_axes(
        ax3,
        width="5%",  # width = 5% of parent_bbox width
        height="50%",  # height : 50%
        loc="lower left",
        bbox_to_anchor=(1.05, 0.0, 1, 2),
        bbox_transform=ax3.transAxes,
        borderpad=0,
    )

    # cax1 = plt.axes([0.1, 0.63, 0.75, 0.02])
    # fig.colorbar(ax1, cax=cax2, orientation='horizontal', cmap='RdBu')
    rot = 270
    lblpd = 15
    cax1 = fig.colorbar(cs1, cax=axin1)
    cax1.ax.set_ylabel(varColStr, rotation=rot, labelpad=lblpd)
    cax2 = fig.colorbar(cs3, cax=axin3)
    cax2.ax.set_ylabel("% difference", rotation=rot, labelpad=lblpd)

    # Resize plots
    plt.subplots_adjust(
        bottom=0.005, left=0.01, right=0.84, top=0.985, hspace=0.01, wspace=0.01
    )

    # pdb.set_trace()
    # plt.show()
    testPath = os.path.join(path)
    if not os.path.exists(testPath):
        os.mkdir(testPath)
    if not os.path.exists(os.path.join(path, var)):
        os.mkdir(os.path.join(path, var))
    fig.savefig(os.path.join(path, var, ".".join([fileName, "png"])), dpi=100)
    plt.close()


# %% Variables

outPathVer = "pngs_v1.1.10"
# outPath = "/p/user_pub/climate_work/durack1/Shared/150219_AMIPForcingData/"  # LLNL/detect
outPath = "/global/homes/d/durack1/git/amipbcs"
# New data
verId = "v1.1.10"
# verPath = "/p/user_pub/climate_work/durack1/"  # LLNL/detect
# verPath = "../"
verPath = "/global/homes/d/durack1/git/amipbcs"
ver = "v20250804"  # Update for each run
verPath = os.path.join(verPath, "input4MIPs/CMIP7/CMIP/PCMDI/PCMDI-AMIP-1-1-10/")
print("verPath:", verPath)
# tos_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-2-0_gn_187001-202108.nc
# Old data
verOldId = "v1.1.9"
verOld = "v20230512"  # Update for each run
# verOldPath = "/p/user_pub/work/input4MIPs/CMIP6Plus/CMIP/PCMDI/PCMDI-AMIP-1-1-9/"  # LLNL/detect
verOldPath = "/global/cfs/projectdirs/m4931/gsharing/user_pub_work"
verOldPath = os.path.join(
    verOldPath, "input4MIPs/CMIP6Plus/CMIP/PCMDI/PCMDI-AMIP-1-1-9/"
)
print("verOldPath:", verOldPath)

f1 = glob.glob("".join([verOldPath, "*/mon/siconc/gn/", verOld, "/*.nc"]))[0]
f2 = glob.glob("".join([verPath, "*/mon/siconc/gn/", ver, "/*.nc"]))[0]
f3 = glob.glob("".join([verOldPath, "*/mon/siconcbcs/gn/", verOld, "/*.nc"]))[0]
f4 = glob.glob("".join([verPath, "*/mon/siconcbcs/gn/", ver, "/*.nc"]))[0]
f5 = glob.glob("".join([verOldPath, "*/mon/tos/gn/", verOld, "/*.nc"]))[0]
f6 = glob.glob("".join([verPath, "*/mon/tos/gn/", ver, "/*.nc"]))[0]
f7 = glob.glob("".join([verOldPath, "*/mon/tosbcs/gn/", verOld, "/*.nc"]))[0]
f8 = glob.glob("".join([verPath, "*/mon/tosbcs/gn/", ver, "/*.nc"]))[0]
ds1 = open_dataset(f1)
ds2 = open_dataset(f2)
ds3 = open_dataset(f3)
ds4 = open_dataset(f4)
ds5 = open_dataset(f5)
ds6 = open_dataset(f6)
ds7 = open_dataset(f7)
ds8 = open_dataset(f8)
# https://scitools.org.uk/cartopy/docs/v0.18/crs/projections.html
x = ds1.lon.data
y = ds1.lat.data

# %% Standard plot - actual and diff maps

# Contour levels
levs1 = list(np.arange(-10, 111, 10))  # siconc
levs2 = list(np.arange(-150, 151, 10))  # siconcbcs
levs3 = list(np.arange(-5, 36, 2.5))  # tos diff
# levs3 = list(np.arange(-0.15, 0.1501, 0.05))  # tos diff
levs4 = list(np.arange(0, 2.1, 0.1))  # % change

# Lab x, y
labX = -140.0
centralLon = 202
labY = 0.0
fntsz = "large"
cmap = "cool"  # 'RdBu'
# https://matplotlib.org/stable/tutorials/colors/colormaps.html

# %% start looping
# Do cleanup
tmpPath = os.path.join(outPath, "pngs", outPathVer, ver)
if os.path.exists(tmpPath):
    shutil.rmtree(tmpPath)
else:
    os.makedirs(tmpPath, exist_ok=False)

for var in ["siconc", "siconcbcs", "tos", "tosbcs"]:
    for yr in np.arange(1870, ds1.time.data[-1].year + 1):  # 1870
        for mn in np.arange(1, 13):
            startTime = "-".join([str(yr), "{:02d}".format(mn), "01"])
            endTime = "-".join([str(yr), "{:02d}".format(mn), "28"])
            print("start:", startTime, "end:", endTime)
            # load into arrays
            match var:
                case "siconc":
                    s1 = eval("ds1.siconc.sel(time=slice(startTime, endTime))")
                    s2 = eval("ds2.siconc.sel(time=slice(startTime, endTime))")
                    lev1 = levs1
                    cmap = "RdBu_r"
                    varColStr = "% coverage"
                case "siconcbcs":
                    s1 = eval("ds3.siconcbcs.sel(time=slice(startTime, endTime))")
                    s2 = eval("ds4.siconcbcs.sel(time=slice(startTime, endTime))")
                    lev1 = levs2
                    cmap = "cool"
                    varColStr = "% coverage"
                case "tos":
                    s1 = eval("ds5.tos.sel(time=slice(startTime, endTime))")
                    s2 = eval("ds6.tos.sel(time=slice(startTime, endTime))")
                    lev1 = levs3
                    cmap = "RdBu_r"
                    varColStr = "degree_C"
                case "tosbcs":
                    s1 = eval("ds7.tosbcs.sel(time=slice(startTime, endTime))")
                    s2 = eval("ds8.tosbcs.sel(time=slice(startTime, endTime))")
                    lev1 = levs3
                    cmap = "cool"
                    varColStr = "degree_C"

            # get time from index
            timeString = "{}{:02d}".format(s1.time.data[0].year, s1.time.data[0].month)
            titleString = "{}{:02d}{}{}".format(
                s1.time.data[0].year, s1.time.data[0].month, " ", var
            )

            plotter(
                s1,
                s2,
                verOldId,
                verId,
                lev1,
                levs4,
                cmap,
                timeString,
                titleString,
                varColStr,
                os.path.join(outPath, "pngs", outPathVer, ver),
                var,
                timeString,
            )
            # pdb.set_trace()
    # end of var - plot video
    out, err = (
        ffmpeg.input(
            os.path.join(outPath, "pngs", outPathVer, ver, var, "*.png"),
            pattern_type="glob",
            framerate=25,
        )
        .output(
            os.path.join(
                outPath,
                "pngs",
                "_".join(["AMIPBCS_newVsOld", var, verId, "".join([ver, ".mp4"])]),
            ),
            crf=20,
            preset="slower",
            movflags="faststart",
            pix_fmt="yuv420p",
        )
        .run()
    )
    # .view(filename='filter_graph')
    # .filter('deflicker', mode='pm', size=10)
    # .filter('scale', size='hd1080', force_original_aspect_ratio='increase')
    # ffmpeg -framerate 48 -i %04d_ESGF-PubStatsPB-MSSans.png 230117_output_48.mp4
