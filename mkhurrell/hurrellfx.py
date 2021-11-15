#!/bin/env python
# -*- coding: utf-8 -*-

"""
Stephen Po-Chedley 15 November 2018

Routines used to calculate the monthly midpoint values for AMIP boundary conditions.

PJD 15 Jul 2019     - Updated createMonthlyMidpoints with varOut argument
PJD 16 Jul 2019     - Updated units assignment outside of grid if block (createMonthlyMidpoints)
PJD  8 Aug 2019     - Further updates to deal with merge conflicts
PJD 23 Sep 2021     - Add NaNf output debug
PJD 28 Sep 2021     - Working with @taylor13 on debugging
PJD 30 Sep 2021     - Update to run through complete grid
PJD 15 Nov 2021     - Update addClimo to ensure last trailing month is December/12
                    - See discussion in https://github.com/PCMDI/amipbcs/issues/23#issuecomment-966610924
PJD 15 Nov 2021     - Update createMonthlyMidpoints with edaysl argument from addClimo

@author: pochedls and durack1
"""

import cdms2
import mkhurrell  # mkhurrell.cpython-37m-x86_64-linux-gnu ; # This occurred the first time it was compiled
import numpy as np

# Control debug output format
np.set_printoptions(formatter={"float": lambda x: "{:8.3f}".format(x)})
from calendar import monthrange

# from matplotlib import pyplot as plt


def getNumDays(time):
    """
    ndays = getNumDays(time)

    Function calculated a time series of the number of days in
    each month based off a cdms2 time vector.
    """
    ndays = np.zeros(len(time), dtype=int)
    timeComponent = time.asComponentTime()
    for i in range(len(time)):
        y = timeComponent[i].year
        m = timeComponent[i].month
        fdays = monthrange(y, m)[1]
        ndays[i] = fdays
    return ndays


def getJacobian(ndays):
    """
    aa, cc = getJacobian(ndays)

    Function calculated the off-diagonals for the jacobian
    used by 'solvmid' using a time series of the number
    of days in each month (ndays).
    """
    ndays = np.insert(ndays, 0, 31)
    ndays = np.append(ndays, 31)
    aa = 2 * ndays[1:-1] / (ndays[1:-1] + ndays[0:-2])
    cc = 2 * ndays[1:-1] / (ndays[1:-1] + ndays[2:])
    return aa, cc


def addClimo(tosi, nyears, ndays, ftype):
    """
    tosi, ndaysp = addClimo(tosi, nyears, ndays, ftype)

    Function adds a one-year climatology (calculated from the first nyears of
    the start and end of the time series) to the start and end of the sst/sic
    time series (tosi). The function will update the ndays vector. The
    climatology decays into the 'real' data using a correlation vector set by
    the data type (ftype = 'ice' or 'sst'). The values of these decorrelation
    vectors was calculated from observed data by Karl Taylor - calculation and
    validation of quantities from the data provided to the function would be a
    useful test

    PJD 15 Nov 2021     - Added edaysl; updated pad to 24 months from start and
                          up to 23 months (if January) to end
    """

    # Create decorrel vectors that are 24-months long
    if ftype == "sst":
        decorrel = [0.68, 0.46, 0.33, 0.26, 0.20, 0.17] + [0.0] * 18
    elif ftype == "ice":
        decorrel = [0.53, 0.23, 0.13, 0.08, 0.05, 0.02] + [0.0] * 18
    else:
        decorrel = [0.0] * 24

    time = tosi.getTime()
    lat = tosi.getLatitude()
    lon = tosi.getLongitude()

    # code to test end indexes
    # lastMonth = time.asComponentTime()[-1].month
    # 9 = 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9; + 10, 11 ,12
    # monVec = np.concatenate([np.arange(1, 13), np.arange(1, 13)], axis=0)
    # monInd = np.arange(23, -1, -1)  # 23:0
    # lastMonInd = np.where(monVec == lastMonth)[0][0]  # First index
    # calOrg = np.concatenate(
    #    [np.arange(lastMonth + 1, 13), np.arange(1, lastMonth + 1)], axis=0
    # )

    # get nyear average climo [12 x lat x lon] - start/end month not relevant
    sclimo = np.mean(
        np.reshape(tosi[0 : nyears * 12, :, :], (nyears, 12, len(lat), len(lon))),
        axis=0,
    )
    eclimo = np.mean(
        np.reshape(tosi[-nyears * 12 :, :, :], (nyears, 12, len(lat), len(lon))), axis=0
    )

    # duplicate climatology, extending to 24-months long
    sclimo = np.tile(sclimo, (2, 1, 1))
    eclimo = np.tile(eclimo, (2, 1, 1))

    # get anomaly map - first/last month minus nyear climatology value
    smap = tosi[0, :, :] - sclimo[0, :, :]
    emap = tosi[-1, :, :] - eclimo[-1, :, :]

    # scale climatology by decorrel vector
    # First repeat start/end anomaly maps for 24 months (np.tile)
    # Expand decorrel into a [24, 1, 1] matrix. If working on the starting
    # anomaly series (sanom), reverse the order of decorrel
    # (i.e., decorrel[::-1] to [0., 0., 0., ..., 0.33, 0.046, 0.68]).
    # Then take the product of these two matrices (the singleton dimensions
    # will broadcast automatically). This yields a 24 months start and end
    # anomaly time series. For sanom the first months will be all zeros (and
    # the end months will be all zeros for eanom) with decorrel values will be
    # applied to the 6 adjacent months to valid data
    sanom = np.tile(smap, (24, 1, 1)) * np.array(
        np.expand_dims(np.expand_dims(decorrel[::-1], 1), 2)
    )
    eanom = np.tile(emap, (24, 1, 1)) * np.array(
        np.expand_dims(np.expand_dims(decorrel, 1), 2)
    )

    # add these anomaly time series to the start and ending climatologies
    sclimo = sclimo + sanom
    eclimo = eclimo + eanom

    # concatenate climatology with start/sclimo and end/eclimo
    tosi = np.concatenate((sclimo, tosi, eclimo), axis=0)

    # get the first and last month (e.g., January = 1 and June = 6)
    smonth = time.asComponentTime()[0].month
    emonth = time.asComponentTime()[-1].month

    # create ndays vector with climatological length values (* 2 = 24 months)
    ndaysclimo = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * 2)

    # tosi starts in January - append a climatology of January through December
    # (times two). We want the number of days in January -> December appended
    # (twice) to the ndays time series. In this case, indices 0 through 11 in
    # ndaysclimo = ([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
    #                31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    sdays = list(ndaysclimo[smonth - 1 : smonth + 11]) * 2

    # tosi ends in June - append a climatology of July through June (times
    # two). We want the number of days in July -> June appended (twice) to the
    # ndays time series. In this case, indices 6 through 17 are required
    # ([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
    #   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    edays = list(ndaysclimo[emonth : emonth + 12]) * 2

    # concatenate all day lengths together
    ndaysp = np.concatenate((sdays, ndays, edays))

    # end timeseries in December. If last month is June, remove last 6 months
    # of the 24-month climatology (July, August, ..., December). The below
    # assumes the start in January 1870 never changes
    tosi = tosi[0:-emonth]
    ndaysp = ndaysp[0:-emonth]
    edaysl = len(range(emonth, 24))

    # tosi and ndaysp are the final padded time series. They both start with
    # two years of climatology. The end with up to 23 months of climatology. If
    # your time series ends in June, you end up with 18 months of climatology.
    # If the time series ends in May, you have 19 months of climatology, etc.
    # variable edaysl provides length of end padding

    return tosi, ndaysp, edaysl


def createMonthlyMidpoints(tosi, ftype, units, nyears, varOut, **kargs):
    """tosimp = createMonthlyMidpoints(tosi, ftype, units, nyears, varOut, **kargs)

    Function creates a time series of monthly midpoint values for gridded
    monthly mean time series. The routine will automatically pad the time
    series with one year of climatology. The climatology is computed from
    the first nyears and last nyears of data. Makes use of Taylor et al.
    (2000) algorithm with parameters set by the type of data (ftype) and
    the units (units). The function will also accept a mask, which will
    determine land points (that are skipped) and a target grid to regrid
    sst or sea ice observations.

    Input arguments:
        tosi[time, lat, lon] - cdms transient variable of sst or sic time series
        ftype - flag to declare type of data ('sst' or 'ice')
        units - units (either 'C'/'Celcius' or 'K'/'Kelvin' for SST
                and '1'/'fraction' for sea ice)
        nyears - integer number of years at the start/end of the timeseries
                 with which to take a climatology (for boundary conditions)
        varOut - string with output filename
    Optional arguments:
        mask[lat, lon] - mask such that values less than one will be skipped
                         (must have shape of the target grid)
        grid - cdms2 grid object used to regrid observations if desired

    Output:
        tosimp - cdms2 transient variable with monthly midpoint values

    Key Reference:

    Taylor, K.E., D. Williamson, and F. Zwiers (2000): The sea surface
          temperature and sea-ice concentration boundary conditions for
          AMIP II simulations. PCMDI Report No. 60 and UCRL-MI-125597,
          Lawrence Livermore National Laboratory, Livermore, CA, 25 pp.

    PJD 16 Jul 2019     - Update units assignment outside of if 'grid' block
    PJD  8 Aug 2019     - Update for merged conflicts
    PJD  2 Nov 2021     - Update bbmin .1 -> 0.01 (consistent with conv)
    PJD  2 Nov 2021     - Correct nitertot counter
    PJD  3 Nov 2021     - Comment solvmid debug statements
    PJD 15 Nov 2021     - Added edaysl from addClimo

    """
    # regrid data if needed
    if "grid" in kargs:
        targetGrid = kargs["grid"]
        diag = {}
        tosi = tosi.regrid(
            targetGrid,
            regridTool="esmf",
            regridMethod="linear",
            missing=np.nan,
            coordSys="deg",
            diag=diag,
            periodicity=1,
        )

    # get axis information
    lat = tosi.getLatitude()
    lon = tosi.getLongitude()
    time = tosi.getTime()
    units = tosi.units

    # deal with optional arguments
    if "mask" in kargs:
        mask = kargs["mask"]
    else:
        mask = np.ones((len(lat), len(lon)))

    # set default values for tmid
    maxiter = 200
    bbmin = 0.01  # 0.1  # 0.001
    jcnt = np.array(-1)
    if ftype == "sst":
        conv = 0.001
        tmax = 400.0
        dt = 0.001
        if (
            (units.lower() == "celcius")
            | (units.lower() == "c")
            | (units.lower() == "degc")
        ):
            tmin = -1.8
        elif (units.lower() == "kelvin") | (units.lower() == "k"):
            tmin = 271.35
    elif ftype == "ice":
        tmin = 0.0
        if (units.lower() == "fraction") | (units.lower() == "1"):
            conv = 0.0001
            tmax = 1.0
        elif (units.lower() == "percent") | (units.lower() == "%"):
            conv = 0.01
            tmax = 100.0
        dt = (tmax - tmin) / 100.0

    maxFlag = np.max(np.max(np.max(tosi))) > tmax
    minFlag = np.min(np.min(np.min(tosi))) < tmin
    if minFlag | maxFlag:
        raise ValueError("Underlying data exceeds limits for " + ftype + " data")

    # construct jacobian for solver
    ndays = getNumDays(time)
    tosip, ndaysp, edaysl = addClimo(tosi, nyears, ndays, ftype)
    aa, cc = getJacobian(ndaysp)

    # pre-allocate output matrix - same size as actual data (no addClimo pad)
    tosimp = np.zeros(tosi.shape)

    # loop over each grid cell and create midpoint values
    sumnotconverg = 0.0
    allresidmax = 0.0
    icnttot = 0
    nitertot = 0
    minall = 0
    maxall = 0
    jjall = 0
    for i in range(len(lat)):
        for j in range(len(lon)):
            obsmean = np.array(tosip[:, i, j])  # extract gridpoint time series
            # lat / lon for diagnostics - not needed?
            alat = lat[i]
            alon = lon[j]
            if mask[i, j] < 0:
                ss = obsmean
            else:

                # Debug for example issues
                # if alat == -77.5 and alon == 182.5:
                # # call solver
                #     (ss, icnt, niter, notconverg, jj, resid, residmax, jumps)\
                #     = mkhurrell.solvmid(alon, alat, conv, dt, tmin, tmax, bbmin,
                #                         maxiter, aa, cc, obsmean, jcnt)
                # elif alat == -76.5 and alon == 186.5:
                # # call solver
                #     (ss, icnt, niter, notconverg, jj, resid, residmax, jumps)\
                #     = mkhurrell.solvmid(alon, alat, conv, dt, tmin, tmax, bbmin,
                #                         maxiter, aa, cc, obsmean, jcnt)
                # else:
                #     continue

                # call solver
                (
                    ss,
                    icnt,
                    niter,
                    notconverg,
                    jj,
                    resid,
                    residmax,
                    jumps,
                ) = mkhurrell.solvmid(
                    alon,
                    alat,
                    conv,
                    dt,
                    tmin,
                    tmax,
                    bbmin,
                    maxiter,
                    aa,
                    cc,
                    obsmean,
                    jcnt,
                )
                # Debug solver output
                # inds = np.where(np.isnan(ss))[0]
                # if len(inds) > 0:
                #    plt.plot(ss[inds[0] - 12 : inds[0] + 12]-10, label="output-10")
                #    plt.plot(obsmean[inds[0] - 12 : inds[0] + 12], label="input")
                #    print(' '.join(["lat:", str(alat), "lon:", str(alon)]))
                #    print('input:')
                #    print(obsmean[inds[0] - 12 : inds[0] + 12])
                #    print('output:')
                #    print(ss[inds[0] - 12 : inds[0] + 12])
                #    plt.title(' '.join(["lat:", str(alat), "lon:", str(alon)]))
                #    plt.legend()
                #    plt.show()
                #    print("stepping..")

            # subset time series (remove padded months) and add to array
            # assumes start is padded with 24 months, end padded by len(edaysl)
            tosimp[:, i, j] = ss[24:-edaysl]

            # test for convergence
            if notconverg > 0:
                print(alat, alon, icnt, niter, notconverg, jj, resid, residmax)
            sumnotconverg = sumnotconverg + notconverg
            if residmax > allresidmax:
                allresidmax = residmax
            icnttot = icnttot + icnt
            nitertot = nitertot + niter
            if jj == -2:
                minall = minall + 1
            elif jj == -1:
                maxall = maxall + 1
            else:
                jjall = jjall + jj

    # Print some diagnostics
    ncells = len(lat) * len(lon)
    nontriv = ncells - minall - maxall
    print()
    print("number of grid cells: ", ncells)
    print("number of cells with non-trivial solutions: ", ncells - minall - maxall)
    print("number of cells with values all = tmax: ", maxall)
    print("number of cells with values all = tmin: ", minall)
    print("number of jumps from tmin to tmax or from tmax to tmin: ", icnttot)
    print("number of cells where observations were smoothed: ", jcnt)
    print(
        "mean number of iterations required (non-trivial cells): ",
        float(nitertot) / float(nontriv),
    )
    print("total number of independent samples: ", jjall)
    print("number of cells where calculation failed to converge: ", sumnotconverg)
    print("maximum residual across all cells and months: ", allresidmax)
    print()

    # create cdms transient variable
    tosimp = cdms2.createVariable(tosimp)
    tosimp.id = varOut
    if varOut == "sst":
        tosimp.standard_name = "sea_surface_temperature"
        tosimp.long_name = "Constructed mid-month Sea Surface Temperature"
    elif varOut == "ice":
        tosimp.standard_name = "sea_ice_concentration"
        tosimp.long_name = "Constructed mid-month Sea-ice concentration"
    tosimp.units = units
    tosimp.setAxis(0, time)
    tosimp.setAxis(1, lat)
    tosimp.setAxis(2, lon)

    return tosimp


def fillVoid(data):
    """fill(data)
    assumes data is of the form [time, lat, lon] and infills zonally and then
    meridionally (setting grid cells with no data equal to adjacent grid
    cells with data).

    This is not meant to be completely physical, but at least some GCMs
    crash if there is no valid data in a given grid.
    """
    # from:
    # https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array

    lat = data.getLatitude()
    lon = data.getLongitude()
    time = data.getTime()
    varId = data.id
    missing = data.missing
    data = np.array(data)
    data[data == missing] = np.nan

    shape = data.shape
    dim = len(shape)
    flagAll = np.zeros(shape, dtype=bool)
    flagAll[:, :, :] = True
    flagAll[np.isnan(data)] = False

    slcs = [slice(None)] * 2

    for i in range(len(time)):
        flag = flagAll[i, :, :]
        dataSlice = data[i, :, :]
        # iterCount = 0
        nlflags = 0
        iterCount = 0
        while np.any(~flag):  # as long as there are any False's in flag
            nflags = np.sum(flag)

            # if working zonally doesn't reduce the number of Falses, infill meridionally
            if nlflags == nflags:
                dim = 0
            else:
                dim = 1

            # make slices to shift view one element along the axis
            slcs1 = slcs[:]
            slcs2 = slcs[:]
            slcs1[dim] = slice(0, -1)
            slcs2[dim] = slice(1, None)

            # replace from the right
            repmask = np.logical_and(~flag[slcs1], flag[slcs2])
            dataSlice[slcs1][repmask] = dataSlice[slcs2][repmask]
            flag[slcs1][repmask] = True

            # replace from the left
            repmask = np.logical_and(~flag[slcs2], flag[slcs1])
            dataSlice[slcs2][repmask] = dataSlice[slcs1][repmask]
            flag[slcs2][repmask] = True
            iterCount += 1
            nlflags = nflags

        data[i] = dataSlice

    data = cdms2.createVariable(data)
    data.id = varId
    data.setAxis(0, time)
    data.setAxis(1, lat)
    data.setAxis(2, lon)

    return data
