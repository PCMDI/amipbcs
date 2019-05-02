Thu 24 Jan 2019 02:30:52 PM PST

The following documents the steps to update the boundary conditions:

1. Install or update Anaconda
> Visit https://www.continuum.io/downloads and follow the prompts for your OS - 2018.12 is latest

[durack1@oceanonly ~]$ conda update --prefix /export/durack1/anaconda2 anaconda
Solving environment: done

# All requested packages already installed.


2a. Install cdms2, cdutil, cd77 (ezget/lats/drs support), CMOR3, NCL, NCO and flex (for wrapit77)
[durack1@oceanonly ~]$ conda create -n cdmsNcd77Ncmor340nclncoflex -c cdat/label/nightly -c conda-forge cdms2 cdutil cd77 cmor ncl nco flex 'python<=3'
# cdms2: 3.1.0.2018.11.29.19.18.g5a83dd4-py27h481b005_0 cdat/label/nightly; cdat_info: 8.0.2018.11.28.19.22.gb6530fa-py_0 cdat/label/nightly;
# cdutil: 8.0.2018.11.12.19.20.g3bd6383-py_0 cdat/label/nightly; cd77: 1.1-py27_0 cdat/label/nightly; cmor: 3.4.0-py27ha570855_0 conda-forge;
# ncl: 6.5.0-blas_openblash04324b8_3 conda-forge [blas_openblas]; nco: 4.7.8-he007dc3_0 conda-forge; flex: 2.6.4-hf484d3e_1004 conda-forge

Installing flex solves issues
***
This solves issues with wrapit77

Hit issues with wrapit77 - 

WRAPIT Version: 120209                                                  
OPERATING SYSTEM: Linux                                                 
nbits = 64                                                              
FORTRAN COMPILER (f90c): gfortran                                       
FORTRAN COMPILER OPTIONS (fopts):  -m64 -fPIC -v -c -fno-second-underscore
/export/durack1/anaconda2/envs/cdat80cmor333nclnco/bin/wrapit77: error while loading shared libraries: libfl.so.2: cannot open shared object file: No such file or directory
FATAL ERROR: wrapit77 failed
***

[durack1@oceanonly ~]$ bash
ANACONDA:/usr/local/anaconda2

You are logged into a Linux Machine...(Version 2.6.32-754.3.5.el6.x86_64)
 Host => oceanonly.llnl.gov

 Hardware: x86_64....
 Using Emacs Bindings...
 11:23:03 up 111 days, 23:55, 10 users,  load average: 0.28, 0.61, 0.60

(psst... may want to install "fortune")

/export/durack1/anaconda2/etc/profile.d/conda.sh
bash: conda: command not found
/export/durack1/anaconda2/etc/profile.d/conda.sh
durack1@oceanonly:[~]:[13752]> source activate cdmsNcd77Ncmor340nclncoflex
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[~]:[13753]>

2b. Validate hdf libraries are correctly installed
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[~]:[13753]> ls -al /export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/lib/libmfhdf.*
-rw-rw-r-- 3 durack1 climate 295K Oct 18 04:20 /export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/lib/libmfhdf.a
lrwxrwxrwx 1 durack1 climate   17 Jan 17 11:23 /export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/lib/libmfhdf.so -> libmfhdf.so.0.0.0
lrwxrwxrwx 1 durack1 climate   17 Jan 17 11:23 /export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/lib/libmfhdf.so.0 -> libmfhdf.so.0.0.0
-rwxrwxr-x 3 durack1 climate 187K Oct 18 04:20 /export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/lib/libmfhdf.so.0.0.0


3. Copy previous directory contents to new directory (NOTE THIS WAS COMPLETED IN OCTOBER, so just replicating what was done back then along with edited scripts)
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9632]> cp -R SST_1-1-4 SST_1-1-5b
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9633]> cd SST_1-1-5b
/work/durack1/Shared/150219_AMIPForcingData/SST_1-1-5b


4. Purge previous files and edit contents to reflect latest downloads
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[SST_1-1-5b]:[9634]> rm -f 180416* ; # Cleanup old files
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[SST_1-1-5b]:[9638]> rm -f MODEL.*.nc ; # Cleanup old files
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[SST_1-1-5b]:[9634]> rsync ../SST_1-1-5/download.sh . ; # Copy files from previous v1.1.5 update
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[SST_1-1-5b]:[9634]> rsync ../SST_1-1-5/SSTICE.Update.unf.csh . ; # Copy files from previous v1.1.5 update


5. Edit script files (there is a need to update PATH info to include the binaries installed into the conda env)
> download.sh & SSTICE.Update.unf.csh
# Make sure to set environment variables NCARG_ROOT and PATH in SSTICE.Update.unf.csh to /export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/bin:${PATH}
# Correct all indexing to reflect download updates


6. Run
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[SST_1-1-5b]:[9642]> download.sh


7. Jump back up to master directory; Edit fortran code
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[SST_1-1-5b]:[9642]> cd ..
> mkhurrell1.f
> increment version number v1.1.4 -> v1.1.5 (for this test use latest data SST_1-1-5b)
> Update start/end dates and array sizes (Should be ok for a while, this update extended to 160 [2030-1870 = 160])


8. Compile
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9643]> cd77 -ezget -lats -cdms -fcray-pointer mkhurrell1.f -o mkhurrell1
f951: Warning: Nonexistent include directory '/export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/include/cdms' [-Wmissing-include-dirs]
f951: Warning: Nonexistent include directory '/export/durack1/anaconda2/envs/cdmsNcd77Ncmor340nclncoflex/lib/libffi-3.1/include' [-Wmissing-include-dirs]
f951: Warning: Nonexistent include directory '/usr/X11R6/include' [-Wmissing-include-dirs]
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9643]> ls -al src
total 44K
drwxr-xr-x  2 durack1 climate 4.0K Apr 19  2017 .
drwxr-xr-x 30 durack1 climate 4.0K Jan 18  2019 ..
-rwxr-xr-x  1 durack1 climate 8.1K Apr 19  2017 drsdef.h
-rwxr-xr-x  1 durack1 climate  14K Jun 16  2015 ketgrib.parms
-rwxr-xr-x  1 durack1 climate 5.5K Apr 19  2017 lats.inc


9. Run
# Create/prepare output directory
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9643]> mkdir 360x180_v1.1.5  
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9639]> mkhurrell1
*****
# If errors are hit, add LD_LIBRARY_PATH to env
mkhurrell1: error while loading shared libraries: libmfhdf.so.0: cannot open shared object file: No such file or directory
# KLUDGE Add libmfhdf.so.0 to path
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5087]> export LD_LIBRARY_PATH=/export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib
*****


10a. Edit sanitize script
> sanitize.py ; # Be aware of complete vs partial years with calendar creation


10b. Register new PCMDI-AMIP-1-1-5 in input4mips-cmor-tables
> Update all details to reflect latest version

10c. Update to latest input4mips-cmor-tables
[durack1@oceanonly input4MIPs-cmor-tables]$ cd
[durack1@oceanonly input4MIPs-cmor-tables]$ cd git/input4MIPs-cmor-tables/
[durack1@oceanonly input4MIPs-cmor-tables]$ git fetch --all -p
[durack1@oceanonly input4MIPs-cmor-tables]$ git checkout master
[durack1@oceanonly input4MIPs-cmor-tables]$ git pull
# Check read perms on input4mips-cmor-tables/Tables subdir


10d. Edit CMOR input jsons
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9639]> ls -al CMOR
total 16K
drwxr-xr-x  2 durack1 climate 4.0K Apr 27  2018 .
drwxr-xr-x 30 durack1 climate 4.0K Jan 18 12:25 ..
-rwxr-xr-x  1 durack1 climate 2.7K Jan  4 13:54 drive_input4MIPs_bcs.json
-rwxr-xr-x  1 durack1 climate 2.5K Jan  4 13:55 drive_input4MIPs_obs.json
lrwxrwxrwx  1 durack1 climate   42 Jun  2  2016 input4MIPs-cmor-tables -> /export/durack1/git/input4MIPs-cmor-tables
# Migrate as much content from CMOR input jsons to source_id registration


11. Cleanup previous v1.1.5 versions
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9632]> rm -rf SST_1-1-5 SST_1-1-5a
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9634]> rm -rf 360x180_v1.1.5_san 360x180_v1.1.5_old 360x180_v1.1.5_old2
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9635]> rm -rf input4MIPs/CMIP6/CMIP/PCMDI/PCMDI-AMIP-1-1-5


12. Edit and run sanitize.py (Write to local directory for testing)
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9635]> python sanitize.py
# Edits required for new version, dates etc
# Note test before writing direct to publication destination


12a. Validate data against previously published version - edit validation script
> make_newVsOldDiffs.py ; # Be wary of any unit changes
# Source cdat80 environment and run movie generation script
(cdat80py2) durack1@oceanonly:[150219_AMIPForcingData]:[11026]> python make_newVsOldDiffs.py > 190118_1322_make_newVsOldDiffs.txt


12b. Validate output movies
# Pay particular attention to end points (2016-)


12c. Edit path and rerun sanitize.py (Now writing directly to destination for publication)
(cdmsNcd77Ncmor340nclncoflex) durack1@oceanonly:[150219_AMIPForcingData]:[9635]> python sanitize.py


13. Publish data to ESGF
> Will need to involve Sasha/Jiwoo for this step


14. Update github repo
> cd /export/durack1/git/amipbcs
# Update all source files
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.f .
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/src . ; # Copy all files required by mkhurrell1.f
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.py .
> chmod 755 *.py *.f ; # Update file permissions
> chmod 644 *.txt *.md ; # Update file permissions
# Update CMOR info
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/CMOR/ CMOR/
> chmod 644 CMOR/*.json ; # Update file permissions
# Update data contributed (*.nc files are not hosted in git - .gitignore)
> #rm -rf SST_1-1-4 ; # Purge local dir (if offical release replaces it)
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/SST_1-1-5b . ; # Copy remote dir
> chmod 755 SST_1-1-5b ; # Update file permissions
> chmod 755 SST_1-1-5b/*.csh SST_1-1-5b/*.dat SST_1-1-5b/*.f SST_1-1-5b/*.ncl SST_1-1-5b/*.sh
> chmod 644 SST_1-1-5b/*.bz2 SST_1-1-5b/*.nc SST_1-1-5b/*.txt
# Update this file once everything is updated
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/README.txt .
# Update README.md
# Add new dir to repo
> git add SST_1-1-5b
# Commit version
> git commit -am 'Updated for new release - 1.1.5; CMOR3.4.0'
> git push
