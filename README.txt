Wed 01 Nov 2017 10:59:39 AM PDT

The following documents the steps to update the boundary conditions:

1. Install Anaconda
> Visit https://www.continuum.io/downloads and follow the prompts for your OS

2. Install CDAT, cd77 (ezget/lats/drs support), CMOR3, NCL and NCO
[durack1@oceanonly ~]$ conda create -n cdat212cmor328Nightlynclnco -c conda-forge -c uvcdat uvcdat
# uvcdat; 2.12-h5103ec1_0; cdms2: 2.12-np113py27_0; cdutil: 2.12-py27_0; vcs: 2.12-py27had75802_0, vtk-cdat: 7.1.0.2.12-py27hd6333ce_0
[durack1@oceanonly ~]$ bash
durack1@oceanonly:[~]:[6835]> source activate cdat212cmor328Nightlynclnco
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[~]:[6819]> ls /export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/lib/libmfhdf.*
/export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/lib/libmfhdf.a
/export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/lib/libmfhdf.la
/export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/lib/libmfhdf.so
/export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/lib/libmfhdf.so.0
/export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/lib/libmfhdf.so.0.0.0
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[~]:[6824]> conda install -c uvcdat cd77
# cd77: 1.0.0-py27_2 uvcdat; ezget: 1.0.1-h2f67867_1 uvcdat; lats: 1.0.1-h9630090_1 uvcdat
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[~]:[6821]> conda install -c pcmdi/label/nightly cmor
# cmor: 2017.10.19.3.2.8-np113py27hd55a59f_0 pcmdi/label/nightly
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[~]:[6822]> conda install -c conda-forge ncl
# ncl: 6.4.0-blas_openblas_4 conda-forge [blas_openblas] and other dependencies
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[~]:[6820]> conda install -c conda-forge nco
# nco: 4.6.9-1 conda-forge

3. Copy previous directory contents to new directory
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6674]> cp -R SST_1-1-2 SST_1-1-3

4. Purge generated files
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[SST_1-1-3]:[6624]> rm -f 170411*.* MODEL.*.nc

5. Edit script files (there is a need to update PATH info to include the binaries installed into the conda env)
> download.sh & SSTICE.Update.unf.csh
# Make sure to set environment variables NCARG_ROOT and PATH in SSTICE.Update.unf.csh to /export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/bin:${PATH}

6. Run
> This step requires NCL 6.4.0 or equivalent installed
(cdat212cmor327nclnco) durack1@oceanonly:[SST_1-1-3]:[5110]> download.sh

7. Edit fortran code
> mkhurrell1.f
> increment version number v1.1.x -> v1.1.x+1
> Update start/end dates and array sizes (Should be ok for a while, this update extended to 160 [2030-1870 = 160])

8. Compile and set environment variables
> This step requires an Anaconda installation with UVCDAT (and cd77 support) and gfortran installed
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5111]> which cd77
/export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/bin/cd77
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5111]> gfortran --version
GNU Fortran (GCC) 4.8.5
# Make sure that drsdef.h, ketgrib.parms and lats.inc are located in the ./src subdir
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5087]> chmod 755 mkhurrell1
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5087]> ./mkhurrell1
mkhurrell1: error while loading shared libraries: libmfhdf.so.0: cannot open shared object file: No such file or directory
# KLUDGE Add libmfhdf.so.0 to path
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5087]> export LD_LIBRARY_PATH=/export/durack1/anaconda2/envs/cdat212cmor328Nightlynclnco/lib

9. Run
# Create/prepare output directory
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5087]> mkdir 360x180_v1.1.3 OR rm -rf 360x180_v1.1.3/* ; # Purge contents if directory exists
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5087]> ./mkhurrell1

10. Edit sanitize script
> sanitize.py ; # Be aware of complete vs partial years with calendar creation

11. Clone (or update) the input4MIPs-cmor-table repo
#[durack1@oceanonly 150219_AMIPForcingData/CMOR]$ git clone https://github.com/PCMDI/input4MIPs-cmor-tables
#durack1@oceanonly:[input4MIPs-cmor-tables]:[master]:[5294]> git pull
#[durack1@oceanonly /export/durack1/git/input4MIPs-cmor-tables/src]$ python writeJson.py

12. Make any required tweaks to CMOR json input
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5552]> ls CMOR
drive_input4MIPs_bcs.json  drive_input4MIPs_obs.json  input4MIPs-cmor-tables

13. Load cmor environment and run
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[4937]> python sanitize.py

14. Validate data against previously published version - edit validation script
> make_newVsOldDiffs.py ; # Be wary of any unit changes

15. Source cdat-nox environment and run movie generation script
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[4937]> source activate uvcdat2120
(uvcdat2120) durack1@oceanonly:[150219_AMIPForcingData]:[5968]> python make_newVsOldDiffs.py > 171031_1632_make_newVsOldDiffs.txt

16. Validate output movies
# Pay particular attention to end points (2016-)

17. Push data to crunchy for ESGF publication
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5553]> chmod 755 -R CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5552]> chmod 644 CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/*/*/PCMDI-AMIP-1-1-3/*/gn/v20171031/*.nc
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5553]> ls -al CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/*/*/PCMDI-AMIP-1-1-3/*/gn/v20171031/*.nc
-rw-r--r-- 1 durack1 climate  46K Oct 31 09:44 CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/fx/ocean/PCMDI-AMIP-1-1-3/areacello/gn/v20171031/areacello_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn.nc
-rw-r--r-- 1 durack1 climate  53K Oct 31 09:44 CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/fx/ocean/PCMDI-AMIP-1-1-3/sftof/gn/v20171031/sftof_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn.nc
-rw-r--r-- 1 durack1 climate 232M Oct 31 09:49 CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/mon/ocean/PCMDI-AMIP-1-1-3/tos/gn/v20171031/tos_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc
-rw-r--r-- 1 durack1 climate 240M Oct 31 09:46 CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/mon/ocean/PCMDI-AMIP-1-1-3/tosbcs/gn/v20171031/tosbcs_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc
-rw-r--r-- 1 durack1 climate  31M Oct 31 09:47 CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/mon/seaIce/PCMDI-AMIP-1-1-3/siconc/gn/v20171031/siconc_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc
-rw-r--r-- 1 durack1 climate  74M Oct 31 09:44 CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/mon/seaIce/PCMDI-AMIP-1-1-3/siconcbcs/gn/v20171031/siconcbcs_input4MIPs_SSTsAndSeaIce_CMIP_PCMDI-AMIP-1-1-3_gn_187001-201706.nc
(cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5542]> rsync -vruta --exclude="PCMDI-AMIP-1-1-0a" CMIP6 crunchy:/work/durack1/Shared/160427_CMIP6_Forcing/
# Make sure to clean up previous problems
# (cdat212cmor328Nightlynclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5542]> rsync -vrutan --delete --exclude="PCMDI-AMIP-1-1-0a" CMIP6 crunchy:/work/durack1/Shared/160427_CMIP6_Forcing/

18. Upload/publish data to ESGF
> Will need to involve Sasha/Tony for this step

19. Update github repo
> cd /export/durack1/git/amipbcs
# Update all source files
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.f .
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/src . ; # Copy all files required by mkhurrell1.f
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.py .
# Purge unused files
> rm -f make_monthsMP.py make_newVsOldDiffsMP.py
# Update CMOR info
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/CMOR/ CMOR/
# Update data contributed (*.nc files are not hosted in git - .gitignore)
> #rm -rf SST_1-1-2 ; # Purge local dir (if offical release replaces it)
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/SST_1-1-3 . ; # Copy remote dir
> chmod 755 SST_1-1-3 ; # Update file permissions
> chmod 755 SST_1-1-3/*.csh SST_1-1-3/*.dat SST_1-1-3/*.f SST_1-1-3/*.ncl SST_1-1-3/*.sh
> chmod 644 SST_1-1-3/*.bz2 SST_1-1-3/*.nc SST_1-1-3/*.txt
# Update this file once everything is updated
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/README.txt .
# Update README.md
# Commit version
> git commit -am 'Updated for new release - 1.1.3; CMOR3.2.8'
> git push