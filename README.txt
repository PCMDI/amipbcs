Thu 08 Sep 2016 09:29:11 AM PDT

The following documents the steps to update the boundary conditions:

1. Install Anaconda
> Visit https://www.continuum.io/downloads and follow the prompts for your OS

2. Install UVCDAT with cd77 (ezget/lats/drs support), CMOR3 and the input4MIPS-cmor-tables
[durack1@oceanonly ~]$ conda create -n uvcdat -c uvcdat/label/nightly -c uvcdat --show-channel-urls uvcdat cd77
[durack1@oceanonly ~]$ conda create -n cmor -c pcmdi cmor -c uvcdat
[durack1@oceanonly 150219_AMIPForcingData/CMOR]$ git clone https://github.com/PCMDI/input4MIPs-cmor-tables

3. Copy the previous directory
[durack1@oceanonly 150219_AMIPForcingData]$ cp -R SST_NEW5 SST_1-1-0a

4. Purge generated files
[durack1@oceanonly SST_1-1-0a]$ rm -f 160526*.* MODEL.*.nc download.sh~ SSTICE.Update.unf.csh~

5. Edit script files
> download.sh
> SSTICE.Update.unf.csh

6. Run
> This step requires NCL 6.3.0 or equivalent is installed
[durack1@oceanonly SST_1-1-0a]$ download.sh

7. Edit fortran code
> mkhurrell1.f
> increment version number v1.1.x -> v1.1.x+1
> Update start/end dates and array sizes

8. Compile
> This step requires an Anaconda installation with UVCDAT (and cd77 support) and gfortran installed
[durack1@oceanonly SST_1-1-0a]$ bash
c durack1@oceanonly:[SST_1-1-0a]:[215]> source activate uvcdat
c (uvcdat)durack1@oceanonly:[SST_1-1-0a]:[215]> which cd77
c /export/durack1/anaconda2/envs/uvcdat/bin/cd77
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[203]> export LD_LIBRARY_PATH=/export/durack1/anaconda2/envs/uvcdat/lib
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[203]> echo $LD_LIBRARY_PATH
c /export/durack1/anaconda2/envs/uvcdat/lib
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[205]> cd77 -ezget -lats -cdms -fcray-pointer mkhurrell1.f -o mkhurrell1
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[204]> rm -rf 360x180_v1.1.0a/*

9. Edit sanitize script
> sanitize.py
> Be aware of complete vs partial years with calendar creation

10. Load cmor environment and run
[durack1@oceanonly SST_1-1-0a]$ bash
c durack1@oceanonly:[SST_1-1-0a]:[215]> source activate cmor
c (cmor)durack1@oceanonly:[150219_AMIPForcingData]:[203]> sanitize.py

11. Validate data against previously published version - edit validation script
> make_monthsMP.py and make_newVsOldDiffsMP.py

12. Run
[durack1@oceanonly 150219_AMIPForcingData]$ make_monthsMP.py

13. Validate output movies

14. Upload/publish data to ESGF
> Will need to involve Sasha for this step

15. Update github repo
> cd /export/durack1/git/amipbcs
# Update all source files
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.f .
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/ketgrib.parms .
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.py .
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.txt .
# Update CMOR info
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/CMOR/ CMOR/
# Update data contributed (*.nc files are not hosted in git - .gitignore)
> #rm -rf SST_NEW5 ; # Purge local dir (if offical release replaces it)
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/SST_1-1-0a . ; # Copy remote dir
# Commit version
> git commit -am 'Updated for new release - 1.1.x'
> git push
