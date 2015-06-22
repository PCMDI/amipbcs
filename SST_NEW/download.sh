#!/bin/tcsh

# Author: Paul J. Durack : pauldurack@llnl.gov
# Created on Wed Apr 22 15:42:19 2015
# @author: durack1

# File written to download all OISSTv2 files
# PJD  9 Apr 2015	- Adapted from ERSST_V3b data
# PJD 22 Apr 2015	- Updated to complete data generation
# PJD 10 Jun 2015 	- Updated to latest data
# USER WILL NEED TO SET:
# ENVIRONMENT VARIABLES: NCARG_ROOT (and PATH to include ncl path)
# PATHS: oi2path

### USER TO SET ###
setenv NCARG_ROOT /work/durack1/Shared/150219_AMIPForcingData/NCL/ncl_ncarg-6.3.0.Linux_RHEL6.4_x86_64_nodap_gcc447
setenv PATH /work/durack1/Shared/150219_AMIPForcingData/NCL/ncl_ncarg-6.3.0.Linux_RHEL6.4_x86_64_nodap_gcc447/bin:${PATH}
set oi2path=/work/durack1/Shared/150219_AMIPForcingData/SST_NEW3/

set date=`date +%y%m%d`
cd ${oi2path}
\rm -r -f ${date} ; # Purge if exists
\mkdir ${date}
cd ${date}

### Step 1 - get most up-to-date files ###
set currentYear=`date +%y`
set acceptList=oiv2mon.20${currentYear}
echo 'downloading '${acceptList}\*.gz
\wget -o ../${date}_log.txt -nv -nc -nH --cut-dirs=3 -rl1 -A ${acceptList}\*.gz --no-check-certificate ftp://ftp.emc.ncep.noaa.gov/cmb/sst/oimonth_v2/
set previousYear=`expr ${currentYear} - 1`
set previousYear=oiv2mon.20${previousYear}
echo 'downloading '${previousYear}\*.gz
\wget -a ../${date}_log.txt -nv -nc -nH --cut-dirs=3 -rl1 -A ${previousYear}\*.gz --no-check-certificate ftp://ftp.emc.ncep.noaa.gov/cmb/sst/oimonth_v2/

### Step 2 - unzip files ###
\gunzip *.gz


### Step 3 - invoke SSTICE.Update.unf.csh ###
# This step requires NCL installed
cd ${oi2path}
./SSTICE.Update.unf.csh ; # REQUIRES EDITING TO POINT TO NEW DATA PATH ${date}

### Step 4 - Interrogate files to make sure things are ok ###
# SST_COMPARE.ncl will generate some slices to peruse

### Step 5 - Extract months of data to append to previous files ###
# Extract only 'new' months, Here 12 months and 'time' index values: 3,16
ncks -O -h -d time,3,16 MODEL.OI2.ice.mnly.201401-201505.unf.nc ICE.update.nc
ncks -O -h -d time,3,16 MODEL.OI2.sst.mnly.201401-201505.unf.nc SST.update.nc
# Make sure updates went correctly ... only the 12 new months, Check 'date' variable
ncdump -v date ICE.update.nc
ncdump -v date SST.update.nc

### Step 6 - Append new months onto existing data ###
# Purge existing files
rm -f MODEL.ICE.HAD187001-198110.OI198111-201505.nc
rm -f MODEL.SST.HAD187001-198110.OI198111-201505.nc
# Merge old and new files
ncrcat ../MODEL.ICE.HAD187001-198110.OI198111-201403.nc ICE.update.nc MODEL.ICE.HAD187001-198110.OI198111-201505.nc
ncrcat ../MODEL.SST.HAD187001-198110.OI198111-201403.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201505.nc
# Make sure updates went correctly
ncdump -v date MODEL.ICE.HAD187001-198110.OI198111-201505.nc
ncdump -v date MODEL.SST.HAD187001-198110.OI198111-201505.nc
# Purge partial new files
rm -f ICE.update.nc ; # These may not purge do to a file handle being unreleased by ncrcat
rm -f MODEL.OI2.ice.mnly.201401-201505.unf.nc
rm -f SST.update.nc
rm -f MODEL.OI2.sst.mnly.201401-201505.unf.nc

### Zip up source data and archive codes ###
echo "Archive source data and purge temp ${date} directory.."
rm -f ${date}.tar.bz2
tar -cjf ${date}.tar.bz2 ${date} ${date}_log.txt download.sh coast_land_NearNbor.f consistent.f lstags.onedeg.dat SST_COMPARE.ncl SSTICE.Update.unf.csh ; # Archive using bzip2 compression
# Extract using >tar -xjf ${date}.tar.bz2
# Purge directory
rm -rf ${date}
echo "${date}_AMIP.nc download complete.."
# Conditionally clean up existing versions of files and archive
if ( $1 != "" ) then
	\tar -cjf ${1}_archive.tar.bz2 ${1}*.* ; # Archive using bzip2 compression
	\rm -f ${1}.tar.bz2
	\rm -f ${1}_log.txt
endif

: <<--
comments
--
