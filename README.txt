Wed 08 Jun 2016 02:38:45 PM PDT

The following documents the steps to update the boundary conditions:

1. Install Anaconda
> Visit https://www.continuum.io/downloads and follow the prompts for your OS

2. Install UVCDAT with cd77 (ezget/lats/drs support) and CMOR
[durack1@oceanonly ~]$ conda create -n uvcdat -c uvcdat/label/nightly -c uvcdat --show-channel-urls uvcdat cd77 ; # Need to check this
[durack1@oceanonly ~]$ conda create -n cmor -c pcmdi cmor -c uvcdat ; # Need to check this

3. Copy the previous directory
[durack1@oceanonly 150219_AMIPForcingData]$ cp -R SST_NEW4/ SST_NEW5/

4. Purge generated files
[durack1@oceanonly SST_NEW5]$ rm -f 160414.tar.bz2 160414_log.txt MODEL*.nc

5. Edit script files
> download.sh
> SSTICE.Update.unf.csh

6. Run
> This step requires NCL 6.3.0 or equivalent is installed
[durack1@oceanonly SST_NEW5]$ download.sh

7. Edit fortran code
> mkhurrell1.f
> increment version number v1.1.x -> v1.1.x+1

8. Compile
> This step requires an Anaconda installation with UVCDAT (and cd77 support) and gfortran installed
[durack1@oceanonly SST_NEW5]$ bash
c durack1@oceanonly:[SST_NEW5]:[215]> source activate uvcdat
c (uvcdat)durack1@oceanonly:[SST_NEW5]:[215]> which cd77
c /export/durack1/anaconda2/envs/uvcdat/bin/cd77
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[203]> export LD_LIBRARY_PATH=/export/durack1/anaconda2/envs/uvcdat/lib
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[203]> echo $LD_LIBRARY_PATH
c /export/durack1/anaconda2/envs/uvcdat/lib
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[205]> cd77 -ezget -lats -cdms -fcray-pointer mkhurrell1.f -o mkhurrell1
c (uvcdat)durack1@oceanonly:[150219_AMIPForcingData]:[204]> rm -rf 360x180_v1.1.x/*

9. Edit sanitize script
> sanitize.py

10. Load cmor environment and run
[durack1@oceanonly SST_NEW5]$ bash
c durack1@oceanonly:[SST_NEW5]:[215]> source activate cmor
c (cmor)durack1@oceanonly:[150219_AMIPForcingData]:[203]> sanitize.py

11. Validate data against previously published version - edit validation script
> make_monthsMP.py and make_newVsOldDiffsMP.py

12. Run
[durack1@oceanonly 150219_AMIPForcingData]$ make_monthsMP.py

13. Validate output movies

14. Upload/publish data to ESGF
> Will need to involve Sasha for this step
