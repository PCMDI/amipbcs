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
# PJD  9 Oct 2017   - Updated for V1.1.3 and using conda environment cdat212cmor327nclnco
# PJD 16 Apr 2018   - Updated for V1.1.4 and using conda env cdat80cmor332nclnco
# PJD 18 Oct 2018   - Updated for V1.1.5 and using conda env cdat80cmor333nclnco
# PJD 18 Oct 2018   - Updated from 201803 to 201809
# PJD 18 Jan 2019   - Updated to reflect latest conda envs, latest downloads and relevant nco indexing
# PJD  2 Jul 2019   - Updated to reflect latest conda envs, latest downloads and relevant nco indexing
# PJD 27 Jul 2021   - Updated download URLs, latest conda envs etc
# PJD  9 Sep 2021   - Updated again for the latest August 2021 data availability
# PJD  4 Nov 2021   - Updated again for the latest September 2021 data availability
# PJD 14 Jun 2022   - Updated belatedly for April 2020 data release; updated to add WORKPATH
# PJD 12 Apr 2023   - Updated for April 2023 data release; CMIP6Plus 1-2-0
# PJD 17 Apr 2023   - Updated for Jan 2023 update; CMIP6Plus 1x1 1.2.0 release
# PJD  3 May 2023   - Updated version number to 1-1-9 after end of line OISSTv2 data identified - will publish as CMIP6Plus

# 1.0 degree data no longer updated
# https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html 1.0 deg
# https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html 0.25 deg

# USER WILL NEED TO SET:
# ENVIRONMENT VARIABLES: NCARG_ROOT (and PATH to include ncl path)
# PATHS: oi2path (and diri, diro and dirm in SSTICE.Update.unf.csh)

# More info: https://climatedataguide.ucar.edu/climate-data/merged-hadley-noaaoi-sea-surface-temperature-sea-ice-concentration-hurrell-et-al-2008
# Doc: http://doi.org/10.1175/2008JCLI2292.1 (Hurrell et al., 2008)

### USER TO SET ###
set condaEnv=230412
setenv NCARG_ROOT /home/durack1/mambaforge/envs/amipbcs${condaEnv}
setenv PATH /home/durack1/mambaforge/envs/amipbcs${condaEnv}/bin:${PATH} ; # Add wrapit77, ncl, nco to PATH
setenv WORKPATH /p/user_pub/climate_work/durack1/Shared/150219_AMIPForcingData

###### UPDATE : PATH REQUIRES UPDATING ##
## NEW DATA ##
set lastYrMn=202301 ; ## UPDATE : NEW DATA END YEAR-MONTH REQUIRES UPDATING ##
set prevYrMn=202201 ; ## UPDATE : NEW DATA START YEAR-MONTH REQUIRES UPDATING ##
set oi2path=${WORKPATH}/SST_1-1-9/ ; ## UPDATE : PATH REQUIRES UPDATING ##
set oi2icefile=MODEL.ICE.HAD187001-198110.OI198111-${lastYrMn}
set oi2sstfile=MODEL.SST.HAD187001-198110.OI198111-${lastYrMn}
set oi2unficefile=MODEL.OI2.ice.mnly.${prevYrMn}-${lastYrMn}.unf
set oi2unfsstfile=MODEL.OI2.sst.mnly.${prevYrMn}-${lastYrMn}.unf
## OLD DATA - are more than two years of data being downloaded? ##
set prevLastYrMn=202205 ; ## UPDATE : PREVIOUS END YEAR-MONTH REQUIRES UPDATING ##
set oi2oldpath=${WORKPATH}/SST_1-1-8/ ; ## UPDATE : PATH REQUIRES UPDATING ##
set oi2oldicefile=MODEL.ICE.HAD187001-198110.OI198111-${prevLastYrMn}
set oi2oldsstfile=MODEL.SST.HAD187001-198110.OI198111-${prevLastYrMn}
######

set date=`date +%y%m%d`
cd ${oi2path}

\rm -r -f ${date} ; # Purge if exists
\mkdir ${date}
cd ${date}

### Step 1 - get most up-to-date files, most bi-yearly updates will only require a single call ###
# 2023
set currentYear=`date +%y`
set acceptList=oiv2mon.20${currentYear}
#set url1=ftp://ftp.emc.ncep.noaa.gov/cmb/sst/oimonth_v2/
set url2=ftp://ftp.cpc.ncep.noaa.gov/precip/PORT/sst/oimonth_v2/
echo 'downloading '${acceptList}\*.gz
\wget -o ../${date}_log.txt -nv -nc -nH --cut-dirs=4 -rl1 -A ${acceptList}\*.gz --no-check-certificate ${url2}
# 2022
set previousYear=`expr ${currentYear} - 1`
set previousYear=oiv2mon.20${previousYear}
echo 'downloading '${previousYear}\*.gz
\wget -a ../${date}_log.txt -nv -nc -nH --cut-dirs=4 -rl1 -A ${previousYear}\*.gz --no-check-certificate ${url2}
# 2021 - not needed, recent updated only includes two current and previous year
#set previousYear=`expr ${currentYear} - 2`
#set previousYear=oiv2mon.20${previousYear}
#echo 'downloading '${previousYear}\*.gz
#\wget -a ../${date}_log.txt -nv -nc -nH --cut-dirs=3 -rl1 -A ${previousYear}\*.gz --no-check-certificate ${url1}

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
#     15,20 (April 2018 through September 2018)
#     3,8 (April 2018 through September 2018) note 0 indexing ; Run 190118
#     9,14 (October 2018 through March 2019) note 0 indexing ; Run 190702
#     3,31 (April 2019 through August 2021) note 0 indexing ; Run 210909
#     3,32 (April 2019 through September 2021) note 0 indexing ; Run 211104
#     3,33 (April 2019 through October 2021) note 0 indexing ; Run 211115
#     10, 16 (November 2021 through May 2022); Run 220614
#     5, 12 (June 2022 through Jan 2023); Run 230412


###### UPDATE : YEARS REQUIRE UPDATING ##
ncks -O -h -d time,5,12 ${oi2unficefile}.nc ICE.update.nc ; # Extract months ## UPDATE : INDEXED MONTHS REQUIRE UPDATING ##
ncks -O -h -d time,5,12 ${oi2unfsstfile}.nc SST.update.nc ; ## UPDATE : INDEXED MONTHS REQUIRE UPDATING - LAST EDIT REQUIRED ##
# Make sure updates went correctly ... only the 12 new months, Check 'date' variable
######
# Check times of new updates
echo '**********'
echo 'Update file time extent'
#ncdump -v date ICE.update.nc
ncdump -v date SST.update.nc
# Check times of older files
echo '**********'
echo 'Previous file time extent'
#ncdump -v date ${oi2oldpath}MODEL.SST.HAD187001-198110.OI198111-201903.nc
ncdump -v date ${oi2oldpath}${oi2oldsstfile}.nc

### Step 6 - Append new months onto existing data ###
# Purge existing files and merge old and new files
rm -f MODEL.ICE.HAD187001-198110.OI198111-*.nc
echo "if rm: No match. - no cleanup required"
ncrcat ${oi2oldpath}${oi2oldicefile}.nc ICE.update.nc ${oi2icefile}.nc
rm -f MODEL.SST.HAD187001-198110.OI198111-*.nc
echo "if rm: No match. - no cleanup required"
ncrcat ${oi2oldpath}${oi2oldsstfile}.nc SST.update.nc ${oi2sstfile}.nc
# Make sure updates went correctly
echo '**********'
echo 'Updated file time extent'
ncdump -v date ${oi2sstfile}.nc
# Purge partial new files
rm -f ICE.update.nc ; # These may not purge do to a file handle being unreleased by ncrcat
rm -f ${oi2unficefile}.nc
rm -f SST.update.nc
rm -f ${oi2unfsstfile}.nc
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