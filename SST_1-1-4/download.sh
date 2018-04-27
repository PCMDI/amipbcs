#!/bin/tcsh

# Author: Paul J. Durack : pauldurack@llnl.gov
# Created on Wed Apr 22 15:42:19 2015
# @author: durack1

# File written to download all OISSTv2 files
# PJD  9 Apr 2015	- Adapted from ERSST_V3b data
# PJD 22 Apr 2015 	- Updated to complete data generation
# PJD 23 Apr 2015 	- Updated URL and DOI
# PJD 14 Apr 2016 	- Updated for V1.0.1 - April 2015 to March 2016 extension
# PJD 25 May 2016 	- Updated for V1.0.1 using correct V1.0.0 input data - May 2015 to April 2016 extension
# PJD 26 May 2016 	- Added wrapit77 to PATH
# PJD 20 Oct 2016 	- Updated for V1.1.1 using V1.1.0 input data - May 2016 to September 2016 extension
# PJD 10 Apr 2017 	- Updated for V1.1.2 using V1.1.1 input data - July 2016 to December 2016 extension
# PJD 10 Apr 2017 	- Oceanonly rebuild required a reinstall of nco.x86_64
# PJD 11 Apr 2017 	- Updated to use new conda environment cdatcmornclnco
# PJD  9 Oct 2017 - Updated for V1.1.3 and using conda environment cdat212cmor327nclnco
# PJD 16 Apr 2018 - Updated for V1.1.4 and using conda env cdat80cmor332nclnco

# USER WILL NEED TO SET:
# ENVIRONMENT VARIABLES: NCARG_ROOT (and PATH to include ncl path)
# PATHS: oi2path (and diri, diro and dirm in SSTICE.Update.unf.csh)

# More info: https://climatedataguide.ucar.edu/climate-data/merged-hadley-noaaoi-sea-surface-temperature-sea-ice-concentration-hurrell-et-al-2008
# Doc: http://doi.org/10.1175/2008JCLI2292.1 (Hurrell et al., 2008)

### USER TO SET ###
setenv NCARG_ROOT /export/durack1/anaconda2/envs/cdat80cmor332nclnco
setenv PATH /export/durack1/anaconda2/envs/cdat80cmor332nclnco/bin:${PATH} ; # Add wrapit77, ncl, nco to PATH

###### UPDATE : PATH REQUIRES UPDATING ##
set oi2path=/work/durack1/Shared/150219_AMIPForcingData/SST_1-1-4/ ; ## UPDATE : PATH REQUIRES UPDATING ##
set oi2oldpath=/work/durack1/Shared/150219_AMIPForcingData/SST_1-1-3/ ; ## UPDATE : PATH REQUIRES UPDATING ##
######

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
###### UPDATE : REQUIRES UPDATING ##
./SSTICE.Update.unf.csh ; ## UPDATE : REQUIRES EDITING TO POINT TO NEW DATA PATH ${date} ##
######


### Step 4 - Interrogate files to make sure things are ok ###
# SST_COMPARE.ncl will generate some slices to peruse


### Step 5 - Extract months of data to append to previous files ###
# Extract only 'new' months, Here 12 months and 'time' index values:
# e.g 3,14 (April 2015 through March 2016)
#     5,15 (June 2015 through April 2016)
#     16,20 (May 2016 through September 2016)
#     9,14 (October 2016 through March 2017) ; Note 0 indexing
#     15,20 (April 2017 through September 2017)
#     9,14 (Sept 2017 through March 2018)
###### UPDATE : YEARS REQUIRE UPDATING ##
ncks -O -h -d time,9,14 MODEL.OI2.ice.mnly.201701-201803.unf.nc ICE.update.nc ; # Extract months ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
ncks -O -h -d time,9,14 MODEL.OI2.sst.mnly.201701-201803.unf.nc SST.update.nc ; ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
# Make sure updates went correctly ... only the 12 new months, Check 'date' variable
######
# Check times of new updates
#ncdump -v date ICE.update.nc
ncdump -v date SST.update.nc
# Check times of older files
ncdump -v date ${oi2oldpath}MODEL.SST.HAD187001-198110.OI198111-201709.nc


### Step 6 - Append new months onto existing data ###
# Purge existing files
#rm -f MODEL.ICE.HAD187001-198110.OI198111-201503.nc ; # Files copied from initial SST_NEW directory - require updating
#rm -f MODEL.SST.HAD187001-198110.OI198111-201503.nc
# Merge old and new files
# V1.0.0 - extend to 201503
#ncrcat ../MODEL.ICE.HAD187001-198110.OI198111-201403.nc ICE.update.nc MODEL.ICE.HAD187001-198110.OI198111-201503.nc
#ncrcat ../MODEL.SST.HAD187001-198110.OI198111-201403.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201503.nc

###### UPDATE : FILENAME/YEARS REQUIRE UPDATING ##
# V1.1.2 - extend to 201703
ncrcat ${oi2oldpath}MODEL.ICE.HAD187001-198110.OI198111-201709.nc ICE.update.nc MODEL.ICE.HAD187001-198110.OI198111-201803.nc ; ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
ncrcat ${oi2oldpath}MODEL.SST.HAD187001-198110.OI198111-201709.nc SST.update.nc MODEL.SST.HAD187001-198110.OI198111-201803.nc ; ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
# Make sure updates went correctly
#ncdump -v date MODEL.ICE.HAD187001-198110.OI198111-201609.nc ; ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
ncdump -v date MODEL.SST.HAD187001-198110.OI198111-201803.nc ; ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
# Purge partial new files
rm -f ICE.update.nc ; # These may not purge do to a file handle being unreleased by ncrcat
rm -f MODEL.OI2.ice.mnly.201701-201803.unf.nc ; ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
rm -f SST.update.nc
rm -f MODEL.OI2.sst.mnly.201701-201803.unf.nc ; ## UPDATE : YEARS IN FILENAMES REQUIRE UPDATING ##
######


### Step 7 - Zip up source data and archive codes ###
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