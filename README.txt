Tue 16 May 2023 04:36:55 PM PDT 

The following documents the steps to update the boundary conditions:

1a. Install or update Anaconda and mamba
> Visit https://www.continuum.io/downloads and follow the prompts for your OS - 2023.03-1 is latest

bash-4.2$ pwd
/home/durack1/git/amipbcs
bash-4.2$ date
Tue May 16 16:44:23 PDT 2023
bash-4.2$ mamba update --prefix /home/durack1/mambaforge conda

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.4.2) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['conda']

conda-forge/noarch                                  12.3MB @   4.9MB/s  2.8s
conda-forge/linux-64                                31.4MB @   4.4MB/s  7.8s

Pinned packages:
  - python 3.9.*


Transaction

  Prefix: /home/durack1/mambaforge

  All requested packages already installed

bash-4.2$ mamba update --prefix /home/durack1/mambaforge mamba

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.4.2) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['mamba']

conda-forge/linux-64                                        Using cache
conda-forge/noarch                                          Using cache

Pinned packages:
  - python 3.9.*


Transaction

  Prefix: /home/durack1/mambaforge

  All requested packages already installed


2a. Install cdms2, cdutil, CMOR3, flex (for wrapit77), gfortran, joblib, matplotlib, pytz, NCL, NCO and ffmpeg-python

bash-4.2$ date
Tue May 16 16:50:29 PDT 2023
bash-4.2$ mamba create -n amipbcs230516 -c conda-forge cmor cdms2 cdutil flex ffmpeg-python gfortran joblib matplotlib pytz ncl nco ++ cartopy & xcdat

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.4.2) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['cmor', 'cdms2', 'cdutil', 'flex', 'ffmpeg-python', 'gfortran', 'joblib', 'matplotlib', 'pytz', 'ncl', 'nco']

conda-forge/linux-64                                        Using cache
conda-forge/noarch                                          Using cache
Transaction

  Prefix: /home/durack1/mambaforge/envs/amipbcs230516

  Updating specs:

   - cmor
   - cdms2
   - cdutil
   - flex
   - ffmpeg-python
   - gfortran
   - joblib
   - matplotlib
   - pytz
   - ncl
   - nco


  Package                           Version  Build                    Channel                    Size
───────────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
───────────────────────────────────────────────────────────────────────────────────────────────────────

  + _libgcc_mutex                       0.1  conda_forge              conda-forge/linux-64     Cached
...
  + cdat_info                         8.2.1  pyhd8ed1ab_2             conda-forge/noarch       Cached
  + cdms2                             3.1.5  py310hb9168da_16         conda-forge/linux-64     Cached
  + cdtime                            3.1.4  py310h87e304a_8          conda-forge/linux-64     Cached
  + cdutil                            8.2.1  pyhd8ed1ab_4             conda-forge/noarch       Cached
...
  + cftime                            1.6.2  py310hde88566_1          conda-forge/linux-64     Cached
...
  + cmor                              3.7.1  py310ha3dcad8_1          conda-forge/linux-64     Cached
...
  + esmf                              8.4.1  mpi_mpich_h7b33e6e_100   conda-forge/linux-64       25MB
  + esmpy                             8.4.1  pyhc1e730c_0             conda-forge/noarch       Cached
...
  + ffmpeg                            5.1.2  gpl_h8dda1f0_106         conda-forge/linux-64     Cached
  + ffmpeg-python                     0.2.0  py_0                     conda-forge/noarch       Cached
 ...
  + gcc                              12.2.0  h26027b1_13              conda-forge/linux-64       27kB
  + gcc_impl_linux-64                12.2.0  hcc96c02_19              conda-forge/linux-64     Cached
  + genutil                           8.2.1  py310h96516ba_2          conda-forge/linux-64     Cached
...
  + gfortran                         12.2.0  h8acd90e_13              conda-forge/linux-64       26kB
  + gfortran_impl_linux-64           12.2.0  h55be85b_19              conda-forge/linux-64     Cached
...
  + libcdms                           3.1.2  h9366c0b_120             conda-forge/linux-64     Cached
...
  + libgcc-devel_linux-64            12.2.0  h3b97bd3_19              conda-forge/linux-64     Cached
  + libgcc-ng                        12.2.0  h65d4601_19              conda-forge/linux-64     Cached
...
  + libgfortran-ng                   12.2.0  h69a702a_19              conda-forge/linux-64     Cached
  + libgfortran5                     12.2.0  h337968e_19              conda-forge/linux-64     Cached
...
  + matplotlib                        3.7.1  py310hff52083_0          conda-forge/linux-64     Cached
...
  + ncl                               6.6.2  ha851dbf_44              conda-forge/linux-64     Cached
  + nco                               5.1.5  h2649ec8_0               conda-forge/linux-64     Cached
...
  + netcdf-fortran                    4.6.0  mpi_mpich_ha3603da_3     conda-forge/linux-64      417kB
  + netcdf4                           1.6.3  nompi_py310h0feb132_100  conda-forge/linux-64     Cached
...
  + python                          3.10.11  he550d4f_0_cpython       conda-forge/linux-64       26MB
...
  + pytz                             2023.3  pyhd8ed1ab_0             conda-forge/noarch       Cached
...
  + udunits2                         2.2.28  hc3e0081_0               conda-forge/linux-64     Cached
...
  + x264                         1!164.3095  h166bdaf_2               conda-forge/linux-64     Cached
  + x265                                3.5  h924138e_3               conda-forge/linux-64     Cached
...
  + zstd                              1.5.2  h3eb15da_6               conda-forge/linux-64     Cached

  Summary:

  Install: 280 packages

  Total download: 66MB



2b. Validate hdf libraries are correctly installed

bash-4.2$ ls -al /home/durack1/mambaforge/envs/amipbcs230516/lib/libmfhdf.*
-rw-rw-r-- 5 durack1 climate 346478 Oct 31  2022 /home/durack1/mambaforge/envs/amipbcs230516/lib/libmfhdf.a
lrwxrwxrwx 1 durack1 climate     17 May 16 16:53 /home/durack1/mambaforge/envs/amipbcs230516/lib/libmfhdf.so -> libmfhdf.so.0.0.0
lrwxrwxrwx 1 durack1 climate     17 May 16 16:53 /home/durack1/mambaforge/envs/amipbcs230516/lib/libmfhdf.so.0 -> libmfhdf.so.0.0.0
-rwxrwxr-x 5 durack1 climate 189872 Oct 31  2022 /home/durack1/mambaforge/envs/amipbcs230516/lib/libmfhdf.so.0.0.0



3. Copy previous directory contents to new directory and purge redundant files (and archive previous v1.1.8 version)

(amipbcs230516) bash-4.2$ date
Thu May 18 09:36:35 PDT 2023
(amipbcs230516) bash-4.2$ pwd
/home/durack1/p-work/Shared/150219_AMIPForcingData
(amipbcs230516) bash-4.2$ mv SST_1-1-9 SST_1-1-9-release
(amipbcs230516) bash-4.2$ mkdir SST_1-1-9
(amipbcs230516) bash-4.2$ cp -R SST_1-1-9-release/* SST_1-1-9/
(amipbcs230516) bash-4.2$ cd SST_1-1-9
(amipbcs230516) bash-4.2$ ls -al
total 931107
drwxr-x---  2 durack1 climate      4096 May 18 09:37 .
drwxrwxr-x. 9 durack1 climate    262144 May 18 09:37 ..
-rw-r-----  1 durack1 climate      2167 May 18 09:37 230503_log.txt
-rw-r-----  1 durack1 climate   2614276 May 18 09:37 230503.tar.bz2
-rwxr-x---  1 durack1 climate      1957 May 18 09:37 coast_land_NearNbor.f
-rwxr-x---  1 durack1 climate      3726 May 18 09:37 consistent.f
-rwxr-x---  1 durack1 climate      8392 May 18 09:37 download.sh
-rwxr-x---  1 durack1 climate    259200 May 18 09:37 lstags.onedeg.dat
-rw-r-----  1 durack1 climate 476202584 May 18 09:37 MODEL.ICE.HAD187001-198110.OI198111-202301.nc
-rw-r-----  1 durack1 climate 476202604 May 18 09:37 MODEL.SST.HAD187001-198110.OI198111-202301.nc
-rwxr-x---  1 durack1 climate      3334 May 18 09:37 SST_COMPARE.ncl
-rwxr-x---  1 durack1 climate     38436 May 18 09:37 SSTICE.Update.unf.csh
(amipbcs230516) bash-4.2$ rm -rf 230503*
(amipbcs230516) bash-4.2$ rm -rf MODEL*.nc
(amipbcs230516) bash-4.2$ ls -al
total 867
drwxr-x---  2 durack1 climate   4096 May 18 09:38 .
drwxrwxr-x. 9 durack1 climate 262144 May 18 09:37 ..
-rwxr-x---  1 durack1 climate   1957 May 18 09:37 coast_land_NearNbor.f
-rwxr-x---  1 durack1 climate   3726 May 18 09:37 consistent.f
-rwxr-x---  1 durack1 climate   8392 May 18 09:37 download.sh
-rwxr-x---  1 durack1 climate 259200 May 18 09:37 lstags.onedeg.dat
-rwxr-x---  1 durack1 climate   3334 May 18 09:37 SST_COMPARE.ncl
-rwxr-x---  1 durack1 climate  38436 May 18 09:37 SSTICE.Update.unf.csh



4. Edit script files (there is a need to update PATH info to include the binaries installed into the conda env)

-> update paths, update file start/end pairs, update env/download date
-> update conda env info in the *.csh file
> download.sh
-> Environment variables set in download.sh, removed from NCL script
-> Correct all indexing for new downloaded input data, and previous version files
> SSTICE.Update.unf.csh



5. Run

(amipbcs230516) bash-4.2$ date
Thu May 18 09:44:37 PDT 2023
(amipbcs230516) bash-4.2$ ./download.sh



6. Edit sanitize and CMORize files

-> Update all files
-> update all version info, input filenames
> ~git/amipbcs/sanitize.py
-> update all version info and path to write data to; version/source_id
-> will need to be registered in the CVs
> ~git/amipbcs/CMOR/drive_input4MIPs_bcs.json
> ~git/amipbcs/CMOR/drive_input4MIPs_obs.json



7. Run sanitize.py

(amipbcs230516) bash-4.2$ date
Thu May 18 09:45:27 PDT 2023
(amipbcs230516) bash-4.2$ pwd
/home/durack1/git/amipbcs
(amipbcs230516) bash-4.2$ cd pcmdiAmipBcs/
(amipbcs230516) bash-4.2$ ./compile 
gfortran compiler on path
--fcompiler=gnu95 GNU Fortran 95 compiler (12.2.0)
Cleanup and compile will begin in 2 seconds..
..
-> run script that will process and CMORize data
(amipbcsTesting220307) bash-4.2$ cd ..
(amipbcsTesting220307) bash-4.2$ python sanitize.py



8. Validate data against previously published version - edit validation script AND install cartopy, xcdat dependencies (Add to step 2a)

(amipbcs230516) bash-4.2$ date
Thu May 18 10:10:01 PDT 2023
(amipbcs230516) bash-4.2$ mamba install -c conda-forge cartopy

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.4.2) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['cartopy']

conda-forge/noarch                                  12.3MB @   4.0MB/s  3.3s
conda-forge/linux-64                                31.5MB @   4.6MB/s  7.5s

Pinned packages:
  - python 3.10.*


Transaction

  Prefix: /home/durack1/mambaforge/envs/amipbcs230516

  Updating specs:

   - cartopy
   - ca-certificates
   - certifi
   - openssl


  Package    Version  Build            Channel                    Size
────────────────────────────────────────────────────────────────────────
  Install:
────────────────────────────────────────────────────────────────────────

  + cartopy   0.21.1  py310hcb7e713_0  conda-forge/linux-64        2MB
  + pooch      1.7.0  pyha770c72_3     conda-forge/noarch       Cached
  + pyproj     3.5.0  py310h15e2413_0  conda-forge/linux-64      447kB
  + pyshp      2.3.1  pyhd8ed1ab_0     conda-forge/noarch       Cached
  + scipy     1.10.1  py310ha4c1d20_3  conda-forge/linux-64       15MB
  + shapely    2.0.1  py310h8b84c32_0  conda-forge/linux-64      422kB

  Summary:

  Install: 6 packages

  Total download: 18MB

(amipbcs230516) bash-4.2$ mamba install -c conda-forge xcdat

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (1.4.2) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['xcdat']

conda-forge/linux-64                                        Using cache
conda-forge/noarch                                          Using cache

Pinned packages:
  - python 3.10.*


Transaction

  Prefix: /home/durack1/mambaforge/envs/amipbcs230516

  Updating specs:

   - xcdat
   - ca-certificates
   - certifi
   - openssl


  Package                  Version  Build                 Channel                    Size
───────────────────────────────────────────────────────────────────────────────────────────
  Install:
───────────────────────────────────────────────────────────────────────────────────────────

  + arrow-cpp               11.0.0  ha770c72_16_cpu       conda-forge/linux-64     Cached
  + aws-c-auth              0.6.26  h2c7c9e7_6            conda-forge/linux-64     Cached
  + aws-c-cal               0.5.26  h71eb795_0            conda-forge/linux-64     Cached
  + aws-c-common            0.8.17  hd590300_0            conda-forge/linux-64     Cached
  + aws-c-compression       0.2.16  h4f47f36_6            conda-forge/linux-64     Cached
  + aws-c-event-stream      0.2.20  h69ce273_6            conda-forge/linux-64     Cached
  + aws-c-http               0.7.7  h7b8353a_3            conda-forge/linux-64     Cached
  + aws-c-io               0.13.21  h2c99d58_4            conda-forge/linux-64      140kB
  + aws-c-mqtt               0.8.6  h3a1964a_15           conda-forge/linux-64     Cached
  + aws-c-s3                 0.2.8  h0933b68_4            conda-forge/linux-64     Cached
  + aws-c-sdkutils           0.1.9  h4f47f36_1            conda-forge/linux-64     Cached
  + aws-checksums           0.1.14  h4f47f36_6            conda-forge/linux-64     Cached
  + aws-crt-cpp             0.19.9  h85076f6_5            conda-forge/linux-64     Cached
  + aws-sdk-cpp            1.10.57  hf40e4db_10           conda-forge/linux-64     Cached
  + bokeh                    3.1.1  pyhd8ed1ab_0          conda-forge/noarch          6MB
  + cf_xarray                0.8.1  pyhd8ed1ab_0          conda-forge/noarch         52kB
  + click                    8.1.3  unix_pyhd8ed1ab_2     conda-forge/noarch       Cached
  + cloudpickle              2.2.1  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + cytoolz                 0.12.0  py310h5764c6d_1       conda-forge/linux-64     Cached
  + dask                  2023.5.0  pyhd8ed1ab_0          conda-forge/noarch          7kB
  + dask-core             2023.5.0  pyhd8ed1ab_0          conda-forge/noarch        845kB
  + distributed           2023.5.0  pyhd8ed1ab_0          conda-forge/noarch        768kB
  + fsspec                2023.5.0  pyh1a96a4e_0          conda-forge/noarch       Cached
  + gflags                   2.2.2  he1b5a44_1004         conda-forge/linux-64     Cached
  + glog                     0.6.0  h6f12383_0            conda-forge/linux-64     Cached
  + importlib_metadata       6.6.0  hd8ed1ab_0            conda-forge/noarch       Cached
  + jinja2                   3.1.2  pyhd8ed1ab_1          conda-forge/noarch       Cached
  + libabseil           20230125.0  cxx17_hcb278e6_1      conda-forge/linux-64     Cached
  + libarrow                11.0.0  h8dc56a0_16_cpu       conda-forge/linux-64     Cached
  + libcrc32c                1.1.2  h9c3ff4c_0            conda-forge/linux-64     Cached
  + libgoogle-cloud          2.8.0  hac9eb74_2            conda-forge/linux-64     Cached
  + libgrpc                 1.54.2  hcf146ea_0            conda-forge/linux-64        6MB
  + libnuma                 2.0.16  h0b41bf4_1            conda-forge/linux-64     Cached
  + libprotobuf            3.21.12  h3eb15da_0            conda-forge/linux-64     Cached
  + libthrift               0.18.1  h5e4af38_0            conda-forge/linux-64     Cached
  + libutf8proc              2.8.0  h166bdaf_0            conda-forge/linux-64     Cached
  + libxslt                 1.1.37  h873f0b0_0            conda-forge/linux-64     Cached
  + locket                   1.0.0  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + lxml                     4.9.2  py310hbdc0903_0       conda-forge/linux-64     Cached
  + lz4                      4.3.2  py310h0cfdcf0_0       conda-forge/linux-64     Cached
  + markupsafe               2.1.2  py310h1fa729e_0       conda-forge/linux-64     Cached
  + msgpack-python           1.0.5  py310hdf3cbec_0       conda-forge/linux-64     Cached
  + orc                      1.8.3  hfdbbad2_0            conda-forge/linux-64     Cached
  + pandas                   2.0.1  py310h7cbd5c2_1       conda-forge/linux-64     Cached
  + parquet-cpp              1.5.1  2                     conda-forge/noarch       Cached
  + partd                    1.4.0  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + psutil                   5.9.5  py310h1fa729e_0       conda-forge/linux-64     Cached
  + pyarrow                 11.0.0  py310he6bfd7f_16_cpu  conda-forge/linux-64     Cached
  + python-tzdata           2023.3  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + pyyaml                     6.0  py310h5764c6d_5       conda-forge/linux-64     Cached
  + re2                 2023.02.02  hcb278e6_0            conda-forge/linux-64     Cached
  + s2n                     1.3.44  h06160fa_0            conda-forge/linux-64      368kB
  + sortedcontainers         2.4.0  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + tblib                    1.7.0  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + toolz                   0.12.0  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + ucx                     1.14.0  h3484d09_2            conda-forge/linux-64     Cached
  + xarray                2023.4.2  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + xcdat                    0.5.0  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + xyzservices           2023.2.0  pyhd8ed1ab_0          conda-forge/noarch       Cached
  + yaml                     0.2.5  h7f98852_2            conda-forge/linux-64     Cached
  + zict                     3.0.0  pyhd8ed1ab_0          conda-forge/noarch       Cached

  Summary:

  Install: 61 packages

  Total download: 14MB

(amipbcs230516) bash-4.2$ pwd
/home/durack1/git/amipbcs
-> edit validation script
> matPlot.py ; # Be wary of any unit changes

(amipbcs230516) bash-4.2$ python matPlot.py



9. Validate output movies
# Pay particular attention to end points (2022-)



10. Publish data to ESGF
> Alert Sasha for this step
Previous v20230512 data has been published



11. Update github repo
> cd /home/durack1/git/amipbcs
# Delete old and update data contributed (*.nc files are not hosted in git - .gitignore)
bash-4.2$ date
Thu May 18 10:22:34 PDT 2023
bash-4.2$ pwd
/home/durack1/git/amipbcs
bash-4.2$ rm -rf SST_1-2-0_1-1-7/ ; # Purge local dir (if offical release replaces it)
bash-4.2$ rsync -vrutaz /p/user_pub/climate_work/durack1/Shared/150219_AMIPForcingData/SST_1-1-9-release . ; # Copy remote dir
sending incremental file list
SST_1-1-9-release/
SST_1-1-9-release/230503.tar.bz2
SST_1-1-9-release/230503_log.txt
SST_1-1-9-release/MODEL.ICE.HAD187001-198110.OI198111-202301.nc
SST_1-1-9-release/MODEL.SST.HAD187001-198110.OI198111-202301.nc
SST_1-1-9-release/SSTICE.Update.unf.csh
SST_1-1-9-release/SST_COMPARE.ncl
SST_1-1-9-release/coast_land_NearNbor.f
SST_1-1-9-release/consistent.f
SST_1-1-9-release/download.sh
SST_1-1-9-release/lstags.onedeg.dat

sent 320,131,521 bytes  received 210 bytes  14,228,076.93 bytes/sec
total size is 955,336,676  speedup is 2.98
bash-4.2$ chmod 755 SST_1-1-9-release ; # Update file permissions
bash-4.2$ chmod 755 SST_*/*.csh SST_*/*.dat SST_*/*.f SST_*/*.ncl SST_*/*.sh
bash-4.2$ chmod 644 SST_*/*.bz2 SST_*/*.nc SST_*/*.txt

# Update this file once everything is updated
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/README.txt .
# Update README.md
# Add new dir to repo
> git add SST_1-1-9-release
# Commit version
> git commit -am 'Updated for new release - 1.1.9; CMOR3.7.1'
> git push
