Fri 27 Apr 2018 01:06:22 PM PDT

The following documents the steps to update the boundary conditions:

1. Install or update Anaconda
> Visit https://www.continuum.io/downloads and follow the prompts for your OS
> durack1@oceanonly:[150219_AMIPForcingData]:[6913]> conda update --prefix /export/durack1/anaconda2 anaconda

2a. Install CDAT, cd77 (ezget/lats/drs support), CMOR3, NCL and NCO
durack1@oceanonly:[150219_AMIPForcingData]:[6918]> conda create -n cdat80cmor332nclnco -c conda-forge -c cdat cdat 'python<=3'
# cdat: 8.0-0 cdat; cdat_info: 8.0-py27_0 conda-forge; cdms2: 3.0-py27_5 conda-forge; cdutil: 8.0-py27_1 conda-forge
durack1@oceanonly:[150219_AMIPForcingData]:[6919]> conda activate cdat80cmor332nclnco
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6920]>

2b. Validate hdf libraries are correctly installed
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6920]> ls /export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib/libmfhdf.*
/export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib/libmfhdf.a
/export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib/libmfhdf.la
/export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib/libmfhdf.so
/export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib/libmfhdf.so.0
/export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib/libmfhdf.so.0.0.0

2c. Install cd77, cmor, ncl, nco
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6921]> conda install -c uvcdat cd77
# cd77: 1.0.0-py27_2 uvcdat
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6918]> conda install -c pcmdi cmor
# cmor: 3.3.2.npy1.13-py27hb8bc26d_0 pcmdi
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6917]> conda install -c conda-forge ncl
# ncl: 6.4.0-blas_openblas_9 conda-forge [blas_openblas] AND OTHER DEPENDENCIES
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6918]> conda install -c conda-forge nco
# nco: 4.7.4-0 conda-forge

3. Copy previous directory contents to new directory
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6909]> cp -R SST_1-1-3 SST_1-1-4

4. Purge generated files and set perms to execute
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6912]> cd SST_1-1-4/
/work/durack1/Shared/150219_AMIPForcingData/SST_1-1-4
(cdat80cmor332nclnco) durack1@oceanonly:[SST_1-1-4]:[6913]> rm -f 171023*.* MODEL.*.nc
(cdat80cmor332nclnco) durack1@oceanonly:[SST_1-1-4]:[6878]> chmod 755 *.*

5. Edit script files (there is a need to update PATH info to include the binaries installed into the conda env)
> download.sh & SSTICE.Update.unf.csh
# Make sure to set environment variables NCARG_ROOT and PATH in SSTICE.Update.unf.csh to /export/durack1/anaconda2/envs/cdat80cmor332nclnco/bin:${PATH}

6. Run
(cdat80cmor332nclnco) durack1@oceanonly:[SST_1-1-4]:[6880]> download.sh

7. Edit fortran code
> mkhurrell1.f
> increment version number v1.1.3 -> v1.1.4
> Update start/end dates and array sizes (Should be ok for a while, this update extended to 160 [2030-1870 = 160])

8. Compile and set environment variables
(cdat80cmor332nclnco) durack1@oceanonly:[SST_1-1-4]:[6880]> which cd77
/export/durack1/anaconda2/envs/cdat80cmor332nclnco/bin/cd77
(cdat80cmor332nclnco) durack1@oceanonly:[SST_1-1-4]:[6864]> gfortran --version
GNU Fortran (GCC) 4.8.5
...
(cdat80cmor332nclnco) durack1@oceanonly:[SST_1-1-4]:[6860]> ls ../src/
drsdef.h  ketgrib.parms  lats.inc
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6832]> cd ..
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6832]> ls -al mkhurrell1.f
-rwxr-xr-x 1 durack1 climate 214K Apr 16 17:06 mkhurrell1.f
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6831]> cd77 -ezget -lats -cdms -fcray-pointer mkhurrell1.f -o mkhurrell1
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6834]> chmod 755 mkhurrell1

9. Run
# Create/prepare output directory
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6825]> mkdir 360x180_v1.1.4
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[6826]> ./mkhurrell1
# If errors are hit, add LD_LIBRARY_PATH to env
mkhurrell1: error while loading shared libraries: libmfhdf.so.0: cannot open shared object file: No such file or directory
# KLUDGE Add libmfhdf.so.0 to path
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5087]> export LD_LIBRARY_PATH=/export/durack1/anaconda2/envs/cdat80cmor332nclnco/lib

10. Edit sanitize script
> sanitize.py ; # Be aware of complete vs partial years with calendar creation

11. Update (or clone) the input4MIPs-cmor-table repo
[durack1@oceanonly input4MIPs-cmor-tables]$ cd
[durack1@oceanonly input4MIPs-cmor-tables]$ cd git/input4MIPs-cmor-tables/
[durack1@oceanonly input4MIPs-cmor-tables]$ git fetch --all -p
[durack1@oceanonly input4MIPs-cmor-tables]$ git checkout master
[durack1@oceanonly input4MIPs-cmor-tables]$ git pull

12. Make any required tweaks to CMOR json input
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[5552]> ls CMOR
drive_input4MIPs_bcs.json  drive_input4MIPs_obs.json  input4MIPs-cmor-tables

13. Load cmor environment and run (this now writes direct to destination for publication)
(cdat80cmor332nclnco) durack1@oceanonly:[150219_AMIPForcingData]:[4937]> python sanitize.py

14. Validate data against previously published version - edit validation script
> make_newVsOldDiffs.py ; # Be wary of any unit changes

15. Source cdat80 environment and run movie generation script
(cdat80) durack1@oceanonly:[150219_AMIPForcingData]:[8079]> python make_newVsOldDiffs.py > 180427_1217_make_newVsOldDiffs.txt

16. Validate output movies
# Pay particular attention to end points (2016-)

18. Publish data to ESGF
> Will need to involve Sasha/Tony for this step

19. Update github repo
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
> #rm -rf SST_1-1-3 ; # Purge local dir (if offical release replaces it)
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/SST_1-1-4 . ; # Copy remote dir
> chmod 755 SST_1-1-4 ; # Update file permissions
> chmod 755 SST_1-1-4/*.csh SST_1-1-4/*.dat SST_1-1-4/*.f SST_1-1-4/*.ncl SST_1-1-4/*.sh
> chmod 644 SST_1-1-4/*.bz2 SST_1-1-4/*.nc SST_1-1-4/*.txt
# Update this file once everything is updated
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/README.txt .
# Update README.md
# Add new dir to repo
> git add SST_1-1-4
# Commit version
> git commit -am 'Updated for new release - 1.1.4; CMOR3.3.2'
> git push
