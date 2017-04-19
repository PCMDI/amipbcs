Tue 18 Apr 2017 05:23:53 PM PDT

The following documents the steps to update the boundary conditions:

1. Install Anaconda
> Visit https://www.continuum.io/downloads and follow the prompts for your OS

2. Install CDAT, cd77 (ezget/lats/drs support), CMOR3, NCL and NCO
[durack1@oceanonly ~]$ conda create -n cdatcmornclnco -c conda-forge -c uvcdat uvcdat cd77
[durack1@oceanonly ~]$ bash
durack1@oceanonly:[~]:[3438]> source activate cdatcmornclnco
(cdatcmornclnco) durack1@oceanonly:[~]:[3439]> conda install -c conda-forge -c pcmdi cmor ; # 3.2.2
(cdatcmornclnco) durack1@oceanonly:[~]:[3440]> conda install -c conda-forge ncl ; # 6.4.0
(cdatcmornclnco) durack1@oceanonly:[~]:[3441]> conda install -c conda-forge nco ; # 4.6.5
(cdatcmornclnco) durack1@oceanonly:[~]:[3442]> conda install spyder ; # 3.1.3

3. Clone the input4MIPs-cmor-table repo
[durack1@oceanonly 150219_AMIPForcingData/CMOR]$ git clone https://github.com/PCMDI/input4MIPs-cmor-tables

4. Copy the previous directory
[durack1@oceanonly 150219_AMIPForcingData]$ cp -R SST_1-1-1 SST_1-1-2

5. Purge generated files
[durack1@oceanonly SST_1-1-2]$ rm -f 161020*.* MODEL.*.nc download.sh~ SSTICE.Update.unf.csh~

6. Edit script files (there is a need to update PATH info to include the binaries installed into the conda env)
> download.sh & SSTICE.Update.unf.csh

7. Run
> This step requires NCL 6.4.0 or equivalent installed
[durack1@oceanonly SST_1-1-2]$ download.sh

8. Edit fortran code
> mkhurrell1.f
> increment version number v1.1.x -> v1.1.x+1
> Update start/end dates and array sizes

9. Compile
> This step requires an Anaconda installation with UVCDAT (and cd77 support) and gfortran installed
[durack1@oceanonly SST_1-1-2]$ bash
durack1@oceanonly:[SST_1-1-2]:[215]> source activate cdatcmornclnco
(cdatcmornclnco)durack1@oceanonly:[SST_1-1-1]:[215]> which cd77
/export/durack1/anaconda2/envs/cdatcmornclnco/bin/cd77
(cdatcmornclnco)durack1@oceanonly:[SST_1-1-1]:[215]> gfortran --version
GNU Fortran (GCC) 4.4.7 20120313 (Red Hat 4.4.7-18)
(cdatcmornclnco)durack1@oceanonly:[150219_AMIPForcingData]:[203]> export LD_LIBRARY_PATH=/export/durack1/anaconda2/envs/uvcdat/lib
(cdatcmornclnco)durack1@oceanonly:[150219_AMIPForcingData]:[203]> echo $LD_LIBRARY_PATH
/export/durack1/anaconda2/envs/cdatcmornclnco/lib
(cdatcmornclnco)durack1@oceanonly:[150219_AMIPForcingData]:[204]> cd77 -ezget -lats -cdms -fcray-pointer mkhurrell1.f -o mkhurrell1
(cdatcmornclnco)durack1@oceanonly:[150219_AMIPForcingData]:[205]> mkdir 360x180_v1.1.2
(cdatcmornclnco)durack1@oceanonly:[150219_AMIPForcingData]:[205]> rm -rf 360x180_v1.1.2/* ; # Purge contents if directory exists

10. Run
(cdatcmornclnco)durack1@oceanonly:[150219_AMIPForcingData]:[205]> mkhurrell1

11. Edit sanitize script
> sanitize.py ; # Be aware of complete vs partial years with calendar creation

12. Load cmor environment and run
[durack1@oceanonly SST_1-1-2]$ bash
durack1@oceanonly:[SST_1-1-2]:[215]> source activate cdatcmornclnco
(cdatcmornclnco)durack1@oceanonly:[150219_AMIPForcingData]:[203]> sanitize.py

13. Validate data against previously published version - edit validation script
> make_monthsMP.py and make_newVsOldDiffsMP.py

14. Source cdat-nox environment and run movie generation script
(cdatcmornclnco) durack1@oceanonly:[150219_AMIPForcingData]:[3239]> source activate uvcdatNightlynox
(uvcdatNightlynox) durack1@oceanonly:[150219_AMIPForcingData]:[3238]>
(uvcdatNightlynox) durack1@oceanonly:[150219_AMIPForcingData]:[3238]> make_monthsMP.py

15. Validate output movies

16. Push data to crunchy for ESGF publication
[durack1@oceanonly 150219_AMIPForcingData]$ rsync -vrutan --exclude="PCMDI-AMIP-1-1-0a" CMIP6 crunchy:/work/durack1/Shared/160427_CMIP6_Forcing/

17. Upload/publish data to ESGF
> Will need to involve Sasha/Tony for this step

18. Update github repo
> cd /export/durack1/git/amipbcs
# Update all source files
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.f .
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/ketgrib.parms .
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/*.py .
# Update CMOR info
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/CMOR/ CMOR/
# Update data contributed (*.nc files are not hosted in git - .gitignore)
> #rm -rf SST_1-1-1 ; # Purge local dir (if offical release replaces it)
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/SST_1-1-2 . ; # Copy remote dir
> chmod 755 SST_1-1-2 ; # Update file permissions
> chmod 755 SST_1-1-2/*.csh SST_1-1-2/*.dat SST_1-1-2/*.f SST_1-1-2/*.ncl SST_1-1-2/*.sh
> chmod 644 SST_1-1-2/*.bz2 SST_1-1-2/*.nc SST_1-1-2/*.txt
# Update this file once everything is updated
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/README.txt .
# Commit version
> git commit -am 'Updated for new release - 1.1.x'
> git push
