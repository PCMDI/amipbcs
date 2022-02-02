Thu 04 Nov 2021 12:31:27 PM PDT

The following documents the steps to update the boundary conditions:

1a. Install or update Anaconda
> Visit https://www.continuum.io/downloads and follow the prompts for your OS - 2021.11 is latest

(base) bash-4.2$ date
Thu Dec  2 14:45:22 PST 2021
(base) bash-4.2$ conda update --prefix /home/durack1/anaconda3 anaconda
Collecting package metadata (current_repodata.json): done
Solving environment: | 
Warning: 2 possible package resolutions (only showing differing packages):
  - defaults/linux-64::libllvm11-11.1.0-h3826bc1_0, defaults/linux-64::llvmlite-0.37.0-py38h295c915_1, defaults/linux-64::numba-0.54.1-py38h51133e4_0, defaults/linux-64::numpy-1.20.3-py38hf144106_0, defaults/linux-64::numpy-base-1.20.3-py38h74d4b33_0, defaults/linux-64::statsmodels-0.12.2-py38h27cfd23_0, defaults/linux-64::tbb-2021.4.0-hd09550d_0
  - defaults/linux-64::libllvm10-10.0.1-hbcb73fb_5, defaults/linux-64::llvmlite-0.36.0-py38h612dafd_4, defaults/linux-64::numba-0.53.1-py38ha9443f7_0, defaults/linux-64::numpy-1.21.2-py38h20f2e39_0, defaults/linux-64::numpy-base-1.21.2-py38h79a1101_0, defaults/linux-64::statsmodels-0.13.0-py38h7f8727e_0, defaults/linux-64::tbb-2020.3-hfd86e86done

## Package Plan ##

  environment location: /home/durack1/anaconda3

  added / updated specs:
    - anaconda


The following packages will be downloaded:

    package                    |            build
    ---------------------------|-----------------
    anaconda-project-0.10.2    |     pyhd3eb1b0_0         218 KB
    arrow-0.13.1               |           py38_0          81 KB
    cycler-0.11.0              |     pyhd3eb1b0_0          12 KB
    debugpy-1.5.1              |   py38h295c915_0         1.7 MB
    fsspec-2021.10.1           |     pyhd3eb1b0_0          96 KB
    idna-3.3                   |     pyhd3eb1b0_0          49 KB
    imagesize-1.3.0            |     pyhd3eb1b0_0           9 KB
    jupyter_client-7.0.6       |     pyhd3eb1b0_0          90 KB
    jupyter_core-4.9.1         |   py38h06a4308_0          75 KB
    libllvm11-11.1.0           |       h3826bc1_0        23.6 MB
    libnghttp2-1.46.0          |       hce63b2e_0         680 KB
    llvmlite-0.37.0            |   py38h295c915_1         431 KB
    matplotlib-3.5.0           |   py38h06a4308_0          28 KB
    matplotlib-base-3.5.0      |   py38h3ed280b_0         5.7 MB
    more-itertools-8.12.0      |     pyhd3eb1b0_0          49 KB
    ncurses-6.3                |       h7f8727e_2         782 KB
    notebook-6.4.6             |   py38h06a4308_0         4.2 MB
    numba-0.54.1               |   py38h51133e4_0         3.3 MB
    openssl-1.1.1l             |       h7f8727e_0         2.5 MB
    packaging-21.3             |     pyhd3eb1b0_0          36 KB
    prometheus_client-0.12.0   |     pyhd3eb1b0_0          47 KB
    pycparser-2.21             |     pyhd3eb1b0_0          94 KB
    pyrsistent-0.18.0          |   py38heee7806_0          94 KB
    pyzmq-22.3.0               |   py38h295c915_2         476 KB
    qtawesome-1.0.3            |     pyhd3eb1b0_0         761 KB
    rope-0.21.1                |     pyhd3eb1b0_0         129 KB
    snowballstemmer-2.2.0      |     pyhd3eb1b0_0          61 KB
    soupsieve-2.3.1            |     pyhd3eb1b0_0          34 KB
    sphinx-4.2.0               |     pyhd3eb1b0_1         1.2 MB
    sqlalchemy-1.4.27          |   py38h7f8727e_0         1.9 MB
    tbb-2021.4.0               |       hd09550d_0         169 KB
    toolz-0.11.2               |     pyhd3eb1b0_0          49 KB
    traitlets-5.1.1            |     pyhd3eb1b0_0          84 KB
    xlsxwriter-3.0.2           |     pyhd3eb1b0_0         111 KB
    zstd-1.5.0                 |       ha4553b6_1         524 KB
    ------------------------------------------------------------
                                           Total:        49.2 MB

The following NEW packages will be INSTALLED:

  libllvm11          pkgs/main/linux-64::libllvm11-11.1.0-h3826bc1_0

The following packages will be REMOVED:

  brunsli-0.1-h2531618_0
  c-blosc2-2.0.4-h5f21a17_1
  cfitsio-3.470-hf0d0db6_6
  charls-2.2.0-h2531618_0
  imagecodecs-2021.8.26-py38hb5ce8f7_1
  jxrlib-1.1-h7b6447c_2
  libaec-1.0.6-h9c3ff4c_0
  libbrotlicommon-1.0.9-h7f98852_6
  libbrotlidec-1.0.9-h7f98852_6
  libbrotlienc-1.0.9-h7f98852_6
  libllvm10-10.0.1-hbcb73fb_5
  libzopfli-1.0.3-he6710b0_0
  openjpeg-2.4.0-h3ad879b_0
  zfp-0.5.5-h2531618_6

The following packages will be UPDATED:

  anaconda-project                      0.10.1-pyhd3eb1b0_0 --> 0.10.2-pyhd3eb1b0_0
  ca-certificates    conda-forge::ca-certificates-2021.10.~ --> pkgs/main::ca-certificates-2021.10.26-h06a4308_2
  cycler             pkgs/main/linux-64::cycler-0.10.0-py3~ --> pkgs/main/noarch::cycler-0.11.0-pyhd3eb1b0_0
  debugpy                              1.4.1-py38h295c915_0 --> 1.5.1-py38h295c915_0
  fsspec                              2021.8.1-pyhd3eb1b0_0 --> 2021.10.1-pyhd3eb1b0_0
  idna                                     3.2-pyhd3eb1b0_0 --> 3.3-pyhd3eb1b0_0
  imagesize                              1.2.0-pyhd3eb1b0_0 --> 1.3.0-pyhd3eb1b0_0
  jupyter_client                         7.0.1-pyhd3eb1b0_0 --> 7.0.6-pyhd3eb1b0_0
  jupyter_core                         4.8.1-py38h06a4308_0 --> 4.9.1-py38h06a4308_0
  libnghttp2         conda-forge::libnghttp2-1.43.0-h812cc~ --> pkgs/main::libnghttp2-1.46.0-hce63b2e_0
  llvmlite                            0.36.0-py38h612dafd_4 --> 0.37.0-py38h295c915_1
  matplotlib                           3.4.3-py38h06a4308_0 --> 3.5.0-py38h06a4308_0
  matplotlib-base                      3.4.3-py38hbbc1b5f_0 --> 3.5.0-py38h3ed280b_0
  more-itertools                        8.10.0-pyhd3eb1b0_0 --> 8.12.0-pyhd3eb1b0_0
  ncurses                                    6.3-h7f8727e_0 --> 6.3-h7f8727e_2
  notebook                             6.4.5-py38h06a4308_0 --> 6.4.6-py38h06a4308_0
  numba                               0.53.1-py38ha9443f7_0 --> 0.54.1-py38h51133e4_0
  packaging                               21.0-pyhd3eb1b0_0 --> 21.3-pyhd3eb1b0_0
  prometheus_client                     0.11.0-pyhd3eb1b0_0 --> 0.12.0-pyhd3eb1b0_0
  pycparser                                       2.20-py_2 --> 2.21-pyhd3eb1b0_0
  pyrsistent                          0.17.3-py38h7b6447c_0 --> 0.18.0-py38heee7806_0
  pyzmq                               22.2.1-py38h295c915_1 --> 22.3.0-py38h295c915_2
  qtawesome                              1.0.2-pyhd3eb1b0_0 --> 1.0.3-pyhd3eb1b0_0
  rope                                  0.19.0-pyhd3eb1b0_0 --> 0.21.1-pyhd3eb1b0_0
  snowballstemmer                        2.1.0-pyhd3eb1b0_0 --> 2.2.0-pyhd3eb1b0_0
  soupsieve                              2.2.1-pyhd3eb1b0_0 --> 2.3.1-pyhd3eb1b0_0
  sphinx                                 4.2.0-pyhd3eb1b0_0 --> 4.2.0-pyhd3eb1b0_1
  sqlalchemy                          1.4.22-py38h7f8727e_0 --> 1.4.27-py38h7f8727e_0
  tbb                                     2020.3-hfd86e86_0 --> 2021.4.0-hd09550d_0
  toolz                                 0.11.1-pyhd3eb1b0_0 --> 0.11.2-pyhd3eb1b0_0
  traitlets                              5.1.0-pyhd3eb1b0_0 --> 5.1.1-pyhd3eb1b0_0
  xlsxwriter                             3.0.1-pyhd3eb1b0_0 --> 3.0.2-pyhd3eb1b0_0
  zstd                   conda-forge::zstd-1.5.0-ha95c52a_0 --> pkgs/main::zstd-1.5.0-ha4553b6_1

The following packages will be SUPERSEDED by a higher-priority channel:

  arrow                                                main --> pkgs/main
  certifi            conda-forge::certifi-2021.10.8-py38h5~ --> pkgs/main::certifi-2021.10.8-py38h06a4308_0
  conda              conda-forge::conda-4.10.3-py38h578d9b~ --> pkgs/main::conda-4.10.3-py38h06a4308_0
  openssl            conda-forge::openssl-1.1.1l-h7f98852_0 --> pkgs/main::openssl-1.1.1l-h7f8727e_0
  tifffile           pkgs/main/noarch::tifffile-2021.7.2-p~ --> pkgs/main/linux-64::tifffile-2020.10.1-py38hdd07704_2

The following packages will be DOWNGRADED:

  numpy                               1.21.2-py38h20f2e39_0 --> 1.20.3-py38hf144106_0
  numpy-base                          1.21.2-py38h79a1101_0 --> 1.20.3-py38h74d4b33_0
  statsmodels                         0.13.0-py38h7f8727e_0 --> 0.12.2-py38h27cfd23_0

(base) bash-4.2$ conda update --prefix /home/durack1/anaconda3 anaconda
Collecting package metadata (current_repodata.json): done
Solving environment: done

# All requested packages already installed.



1b. Install or update mamba
> Visit https://anaconda.org/conda-forge/mamba - 0.19.0 is latest

(base) bash-4.2$ date
Thu Dec  2 14:49:21 PST 2021
(base) bash-4.2$ conda install -n base -c conda-forge mamba="0.19.0"
Collecting package metadata (current_repodata.json): done
Solving environment: done

## Package Plan ##

  environment location: /home/durack1/anaconda3

  added / updated specs:
    - mamba=0.19.0


The following packages will be downloaded:

    package                    |            build
    ---------------------------|-----------------
    conda-4.11.0               |   py38h578d9bd_0        16.9 MB  conda-forge
    curl-7.80.0                |       h2574ce0_0         153 KB  conda-forge
    libcurl-7.80.0             |       h2574ce0_0         337 KB  conda-forge
    libmamba-0.19.0            |       h3985d26_0         1.5 MB  conda-forge
    libmambapy-0.19.0          |   py38h908000c_0         313 KB  conda-forge
    mamba-0.19.0               |   py38h1abaa86_0          41 KB  conda-forge
    pybind11-abi-4             |       hd8ed1ab_3          10 KB  conda-forge
    yaml-cpp-0.6.3             |       he1b5a44_4         208 KB  conda-forge
    ------------------------------------------------------------
                                           Total:        19.4 MB

The following NEW packages will be INSTALLED:

  libmamba           conda-forge/linux-64::libmamba-0.19.0-h3985d26_0
  libmambapy         conda-forge/linux-64::libmambapy-0.19.0-py38h908000c_0
  pybind11-abi       conda-forge/noarch::pybind11-abi-4-hd8ed1ab_3
  yaml-cpp           conda-forge/linux-64::yaml-cpp-0.6.3-he1b5a44_4

The following packages will be UPDATED:

  certifi            pkgs/main::certifi-2021.10.8-py38h06a~ --> conda-forge::certifi-2021.10.8-py38h578d9bd_1
  conda              pkgs/main::conda-4.10.3-py38h06a4308_0 --> conda-forge::conda-4.11.0-py38h578d9bd_0
  curl                                    7.79.1-h2574ce0_1 --> 7.80.0-h2574ce0_0
  libcurl                                 7.79.1-h2574ce0_1 --> 7.80.0-h2574ce0_0
  mamba                               0.17.0-py38h2aa5da1_0 --> 0.19.0-py38h1abaa86_0

The following packages will be SUPERSEDED by a higher-priority channel:

  ca-certificates    pkgs/main::ca-certificates-2021.10.26~ --> conda-forge::ca-certificates-2021.10.8-ha878542_0
  openssl              pkgs/main::openssl-1.1.1l-h7f8727e_0 --> conda-forge::openssl-1.1.1l-h7f98852_0



2a. Install cdms2, cdutil, CMOR3, flex (for wrapit77), gfortran, joblib, matplotlib, pytz, NCL, and NCO
-> Also think about adding vcs

(base) bash-4.2$ date
Fri Dec  3 09:27:32 PST 2021
(base) bash-4.2$ mamba create  -n amipbcs211203 -c conda-forge cdms2 cdutil cmor flex gfortran joblib matplotlib pytz ncl nco

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

        mamba (0.19.0) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████


Looking for: ['cdms2', 'cdutil', 'cmor', 'flex', 'gfortran', 'joblib', 'matplotlib', 'pytz', 'ncl', 'nco']

pkgs/r/linux-64          [====================] (00m:00s) No change
pkgs/r/noarch            [====================] (00m:00s) No change
pkgs/main/noarch         [====================] (00m:00s) Done
pkgs/main/linux-64       [====================] (00m:00s) Done
conda-forge/noarch       [====================] (00m:02s) Done
conda-forge/linux-64     [====================] (00m:05s) Done
Transaction

  Prefix: /home/durack1/anaconda3/envs/amipbcs211203

  Updating specs:

   - cdms2
   - cdutil
   - cmor
   - flex
   - gfortran
   - joblib
   - matplotlib
   - pytz
   - ncl
   - nco


  Package                           Version  Build                        Channel                    Size
───────────────────────────────────────────────────────────────────────────────────────────────────────────
  Install:
───────────────────────────────────────────────────────────────────────────────────────────────────────────

  + _libgcc_mutex                       0.1  conda_forge                  conda-forge/linux-64     Cached
  + _openmp_mutex                       4.5  1_gnu                        conda-forge/linux-64     Cached
  + alsa-lib                          1.2.3  h516909a_0                   conda-forge/linux-64     Cached
  + attrs                            21.2.0  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + binutils_impl_linux-64           2.36.1  h193b22a_2                   conda-forge/linux-64      10 MB
  + binutils_linux-64                  2.36  hf3e587d_1                   conda-forge/linux-64      23 KB
  + boost-cpp                        1.74.0  h312852a_4                   conda-forge/linux-64     Cached
  + brotli                            1.0.9  h7f98852_6                   conda-forge/linux-64     Cached
  + brotli-bin                        1.0.9  h7f98852_6                   conda-forge/linux-64     Cached
  + brotlipy                          0.7.0  py310h6acc77f_1003           conda-forge/linux-64     Cached
  + bzip2                             1.0.8  h7f98852_4                   conda-forge/linux-64     Cached
  + c-ares                           1.18.1  h7f98852_0                   conda-forge/linux-64     Cached
  + ca-certificates               2021.10.8  ha878542_0                   conda-forge/linux-64     Cached
  + cairo                            1.16.0  h6cf1ce9_1008                conda-forge/linux-64     Cached
  + cdat_info                         8.2.1  pyhd8ed1ab_2                 conda-forge/noarch       Cached
  + cdms2                             3.1.5  py310h65abc81_14             conda-forge/linux-64     Cached
  + cdtime                            3.1.4  py310h1a17f1e_7              conda-forge/linux-64     Cached
  + cdutil                            8.2.1  pyhd8ed1ab_1                 conda-forge/noarch       Cached
  + certifi                       2021.10.8  py310hff52083_1              conda-forge/linux-64     Cached
  + cffi                             1.15.0  py310h0fdd8cc_0              conda-forge/linux-64     Cached
  + cfitsio                           4.0.0  h9a35b8e_0                   conda-forge/linux-64     Cached
  + cftime                          1.5.1.1  py310h96516ba_1              conda-forge/linux-64     Cached
  + charset-normalizer                2.0.8  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + cmor                              3.6.1  py310hc9dadf1_5              conda-forge/linux-64     Cached
  + cryptography                     36.0.0  py310h685ca39_0              conda-forge/linux-64     Cached
  + curl                             7.80.0  h2574ce0_0                   conda-forge/linux-64     Cached
  + cycler                           0.11.0  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + dbus                             1.13.6  h48d8840_2                   conda-forge/linux-64     Cached
  + distarray                        2.12.2  pyhd8ed1ab_2                 conda-forge/noarch       Cached
  + esmf                              8.2.0  mpi_mpich_h4975321_100       conda-forge/linux-64     Cached
  + esmpy                             8.2.0  mpi_mpich_py310hd9c82d4_101  conda-forge/linux-64     Cached
  + expat                             2.4.1  h9c3ff4c_0                   conda-forge/linux-64     Cached
  + flex                              2.6.4  h58526e2_1004                conda-forge/linux-64     Cached
  + font-ttf-dejavu-sans-mono          2.37  hab24e00_0                   conda-forge/noarch       Cached
  + font-ttf-inconsolata              3.000  h77eed37_0                   conda-forge/noarch       Cached
  + font-ttf-source-code-pro          2.038  h77eed37_0                   conda-forge/noarch       Cached
  + font-ttf-ubuntu                    0.83  hab24e00_0                   conda-forge/noarch       Cached
  + fontconfig                       2.13.1  hba837de_1005                conda-forge/linux-64     Cached
  + fonts-conda-ecosystem                 1  0                            conda-forge/noarch       Cached
  + fonts-conda-forge                     1  0                            conda-forge/noarch       Cached
  + fonttools                        4.28.3  py310h6acc77f_0              conda-forge/linux-64       2 MB
  + freeglut                          3.2.1  h9c3ff4c_2                   conda-forge/linux-64     Cached
  + freetype                         2.10.4  h0708190_1                   conda-forge/linux-64     Cached
  + freexl                            1.0.6  h7f98852_0                   conda-forge/linux-64     Cached
  + future                           0.18.2  py310hff52083_4              conda-forge/linux-64     Cached
  + g2clib                            1.6.3  heb9ad7a_1                   conda-forge/linux-64     Cached
  + gcc                              11.2.0  h702ea55_1                   conda-forge/linux-64      23 KB
  + gcc_impl_linux-64                11.2.0  h82a94d6_11                  conda-forge/linux-64      52 MB
  + gcc_linux-64                     11.2.0  h39a9532_1                   conda-forge/linux-64      24 KB
  + genutil                           8.2.1  py310h96516ba_2              conda-forge/linux-64     Cached
  + geos                             3.10.0  h9c3ff4c_0                   conda-forge/linux-64     Cached
  + geotiff                           1.7.0  hcfb7246_3                   conda-forge/linux-64     Cached
  + gettext                        0.19.8.1  h73d1719_1008                conda-forge/linux-64     Cached
  + gfortran                         11.2.0  h8811e0c_1                   conda-forge/linux-64      23 KB
  + gfortran_impl_linux-64           11.2.0  h7a446d4_11                  conda-forge/linux-64      15 MB
  + gfortran_linux-64                11.2.0  h777b47f_1                   conda-forge/linux-64      23 KB
  + giflib                            5.2.1  h36c2ea0_2                   conda-forge/linux-64     Cached
  + glib                             2.70.1  h780b84a_0                   conda-forge/linux-64     Cached
  + glib-tools                       2.70.1  h780b84a_0                   conda-forge/linux-64     Cached
  + gsl                                 2.7  he838d99_0                   conda-forge/linux-64     Cached
  + gst-plugins-base                 1.18.5  hf529b03_2                   conda-forge/linux-64     Cached
  + gstreamer                        1.18.5  h9f60fe5_2                   conda-forge/linux-64     Cached
  + hdf4                             4.2.15  h10796ff_3                   conda-forge/linux-64     Cached
  + hdf5                             1.12.1  mpi_mpich_h9c45103_2         conda-forge/linux-64     Cached
  + hdfeos2                            2.20  h64bfcee_1000                conda-forge/linux-64     Cached
  + hdfeos5                          5.1.16  h3fd3831_9                   conda-forge/linux-64     Cached
  + icu                                68.2  h9c3ff4c_0                   conda-forge/linux-64     Cached
  + idna                                3.1  pyhd3deb0d_0                 conda-forge/noarch       Cached
  + importlib-metadata                4.8.2  py310hff52083_0              conda-forge/linux-64     Cached
  + importlib_resources               5.4.0  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + ipython_genutils                  0.2.0  py_1                         conda-forge/noarch       Cached
  + jasper                           2.0.33  ha77e612_0                   conda-forge/linux-64     Cached
  + jbig                                2.1  h7f98852_2003                conda-forge/linux-64     Cached
  + joblib                            1.1.0  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + jpeg                                 9d  h36c2ea0_0                   conda-forge/linux-64     Cached
  + json-c                             0.15  h98cffda_0                   conda-forge/linux-64     Cached
  + jsonschema                        4.2.1  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + jupyter_core                      4.9.1  py310hff52083_1              conda-forge/linux-64     Cached
  + kealib                           1.4.14  h87e4c3c_3                   conda-forge/linux-64     Cached
  + kernel-headers_linux-64          2.6.32  he073ed8_15                  conda-forge/noarch       707 KB
  + kiwisolver                        1.3.2  py310h91b1402_1              conda-forge/linux-64     Cached
  + krb5                             1.19.2  hcc1bbae_3                   conda-forge/linux-64     Cached
  + lazy-object-proxy                 1.6.0  py310h6acc77f_1              conda-forge/linux-64     Cached
  + lcms2                              2.12  hddcbb42_0                   conda-forge/linux-64     Cached
  + ld_impl_linux-64                 2.36.1  hea4e1c9_2                   conda-forge/linux-64     Cached
  + lerc                                3.0  h9c3ff4c_0                   conda-forge/linux-64     Cached
  + libblas                           3.9.0  12_linux64_openblas          conda-forge/linux-64     Cached
  + libbrotlicommon                   1.0.9  h7f98852_6                   conda-forge/linux-64     Cached
  + libbrotlidec                      1.0.9  h7f98852_6                   conda-forge/linux-64     Cached
  + libbrotlienc                      1.0.9  h7f98852_6                   conda-forge/linux-64     Cached
  + libcblas                          3.9.0  12_linux64_openblas          conda-forge/linux-64     Cached
  + libcdms                           3.1.2  h3bdf4ba_117                 conda-forge/linux-64     Cached
  + libcf                             1.0.3  py310hdc0ccca_113            conda-forge/linux-64     Cached
  + libclang                         11.1.0  default_ha53f305_1           conda-forge/linux-64     Cached
  + libcurl                          7.80.0  h2574ce0_0                   conda-forge/linux-64     Cached
  + libdap4                          3.20.6  hd7c4107_2                   conda-forge/linux-64     Cached
  + libdeflate                          1.8  h7f98852_0                   conda-forge/linux-64     Cached
  + libdrs                            3.1.2  hf593df3_118                 conda-forge/linux-64     Cached
  + libdrs_f                          3.1.2  h7e76ec7_114                 conda-forge/linux-64     Cached
  + libedit                    3.1.20191231  he28a2e2_2                   conda-forge/linux-64     Cached
  + libev                              4.33  h516909a_1                   conda-forge/linux-64     Cached
  + libevent                         2.1.10  h9b69904_4                   conda-forge/linux-64     Cached
  + libffi                            3.4.2  h7f98852_5                   conda-forge/linux-64     Cached
  + libgcc-devel_linux-64            11.2.0  h0952999_11                  conda-forge/linux-64       3 MB
  + libgcc-ng                        11.2.0  h1d223b6_11                  conda-forge/linux-64     Cached
  + libgdal                           3.3.3  h18e3bf0_2                   conda-forge/linux-64     Cached
  + libgfortran-ng                   11.2.0  h69a702a_11                  conda-forge/linux-64     Cached
  + libgfortran5                     11.2.0  h5c6108e_11                  conda-forge/linux-64     Cached
  + libglib                          2.70.1  h174f98d_0                   conda-forge/linux-64     Cached
  + libglu                            9.0.0  he1b5a44_1001                conda-forge/linux-64     Cached
  + libgomp                          11.2.0  h1d223b6_11                  conda-forge/linux-64     Cached
  + libiconv                           1.16  h516909a_0                   conda-forge/linux-64     Cached
  + libkml                            1.3.0  h238a007_1014                conda-forge/linux-64     Cached
  + liblapack                         3.9.0  12_linux64_openblas          conda-forge/linux-64     Cached
  + libllvm11                        11.1.0  hf817b99_2                   conda-forge/linux-64     Cached
  + libnetcdf                         4.8.1  mpi_mpich_h319fa22_1         conda-forge/linux-64     Cached
  + libnghttp2                       1.43.0  h812cca2_1                   conda-forge/linux-64     Cached
  + libnsl                            2.0.0  h7f98852_0                   conda-forge/linux-64     Cached
  + libogg                            1.3.4  h7f98852_1                   conda-forge/linux-64     Cached
  + libopenblas                      0.3.18  pthreads_h8fe5266_0          conda-forge/linux-64     Cached
  + libopus                           1.3.1  h7f98852_1                   conda-forge/linux-64     Cached
  + libpng                           1.6.37  h21135ba_2                   conda-forge/linux-64     Cached
  + libpq                              13.5  hd57d9b9_0                   conda-forge/linux-64     Cached
  + librttopo                         1.1.0  h0ad649c_7                   conda-forge/linux-64     Cached
  + libsanitizer                     11.2.0  he4da1e4_11                  conda-forge/linux-64       6 MB
  + libspatialite                     5.0.1  h1d9e4f1_10                  conda-forge/linux-64     Cached
  + libssh2                          1.10.0  ha56f1ee_2                   conda-forge/linux-64     Cached
  + libstdcxx-ng                     11.2.0  he4da1e4_11                  conda-forge/linux-64     Cached
  + libtiff                           4.3.0  h6f004c6_2                   conda-forge/linux-64     Cached
  + libuuid                          2.32.1  h7f98852_1000                conda-forge/linux-64     Cached
  + libvorbis                         1.3.7  h9c3ff4c_0                   conda-forge/linux-64     Cached
  + libwebp-base                      1.2.1  h7f98852_0                   conda-forge/linux-64     Cached
  + libxcb                             1.13  h7f98852_1004                conda-forge/linux-64     Cached
  + libxkbcommon                      1.0.3  he3ba5ed_0                   conda-forge/linux-64     Cached
  + libxml2                          2.9.12  h72842e0_0                   conda-forge/linux-64     Cached
  + libzip                            1.8.0  h4de3113_1                   conda-forge/linux-64     Cached
  + libzlib                          1.2.11  h36c2ea0_1013                conda-forge/linux-64     Cached
  + lz4-c                             1.9.3  h9c3ff4c_1                   conda-forge/linux-64     Cached
  + m4                               1.4.18  h516909a_1001                conda-forge/linux-64     Cached
  + matplotlib                        3.5.0  py310hff52083_0              conda-forge/linux-64     Cached
  + matplotlib-base                   3.5.0  py310h23f4a51_0              conda-forge/linux-64     Cached
  + mpi                                 1.0  mpich                        conda-forge/linux-64     Cached
  + mpi4py                            3.1.3  py310h853ac07_0              conda-forge/linux-64     Cached
  + mpich                             3.4.2  h846660c_100                 conda-forge/linux-64     Cached
  + munkres                           1.1.4  pyh9f0ad1d_0                 conda-forge/noarch       Cached
  + mysql-common                     8.0.27  ha770c72_1                   conda-forge/linux-64     Cached
  + mysql-libs                       8.0.27  hfa10184_1                   conda-forge/linux-64     Cached
  + nbformat                          5.1.3  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + ncl                               6.6.2  h50d29c7_33                  conda-forge/linux-64     Cached
  + nco                               5.0.3  hcfc2ecc_0                   conda-forge/linux-64     Cached
  + ncurses                             6.2  h58526e2_4                   conda-forge/linux-64     Cached
  + netcdf-fortran                    4.5.3  mpi_mpich_h1364a43_6         conda-forge/linux-64     Cached
  + netcdf4                           1.5.8  nompi_py310hd7ca5b8_101      conda-forge/linux-64     Cached
  + nspr                               4.32  h9c3ff4c_1                   conda-forge/linux-64     Cached
  + nss                                3.73  hb5efdd6_0                   conda-forge/linux-64     Cached
  + numpy                            1.21.4  py310h57288b1_0              conda-forge/linux-64     Cached
  + olefile                            0.46  pyh9f0ad1d_1                 conda-forge/noarch       Cached
  + openblas                         0.3.18  pthreads_h4748800_0          conda-forge/linux-64     Cached
  + openjpeg                          2.4.0  hb52868f_1                   conda-forge/linux-64     Cached
  + openssl                          1.1.1l  h7f98852_0                   conda-forge/linux-64     Cached
  + packaging                          21.3  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + pcre                               8.45  h9c3ff4c_0                   conda-forge/linux-64     Cached
  + pillow                            8.4.0  py310h07f4688_0              conda-forge/linux-64     Cached
  + pip                              21.3.1  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + pixman                           0.40.0  h36c2ea0_0                   conda-forge/linux-64     Cached
  + poppler                         21.09.0  ha39eefc_3                   conda-forge/linux-64     Cached
  + poppler-data                     0.4.11  hd8ed1ab_0                   conda-forge/noarch       Cached
  + postgresql                         13.5  h2510834_0                   conda-forge/linux-64     Cached
  + proj                              8.1.1  h277dcde_2                   conda-forge/linux-64     Cached
  + pthread-stubs                       0.4  h36c2ea0_1001                conda-forge/linux-64     Cached
  + pycparser                          2.21  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + pyopenssl                        21.0.0  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + pyparsing                         3.0.6  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + pyqt                             5.12.3  py310hff52083_8              conda-forge/linux-64     Cached
  + pyqt-impl                        5.12.3  py310h1f8e252_8              conda-forge/linux-64     Cached
  + pyqt5-sip                       4.19.18  py310h122e73d_8              conda-forge/linux-64     Cached
  + pyqtchart                          5.12  py310hfcd6d55_8              conda-forge/linux-64     Cached
  + pyqtwebengine                    5.12.1  py310hfcd6d55_8              conda-forge/linux-64     Cached
  + pyrsistent                       0.18.0  py310h6acc77f_0              conda-forge/linux-64     Cached
  + pysocks                           1.7.1  py310hff52083_4              conda-forge/linux-64     Cached
  + python                           3.10.0  h62f1059_3_cpython           conda-forge/linux-64     Cached
  + python-dateutil                   2.8.2  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + python_abi                         3.10  2_cp310                      conda-forge/linux-64     Cached
  + pytz                             2021.3  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + qt                               5.12.9  hda022c4_4                   conda-forge/linux-64     Cached
  + readline                            8.1  h46c0cb4_0                   conda-forge/linux-64     Cached
  + requests                         2.26.0  pyhd8ed1ab_1                 conda-forge/noarch       Cached
  + setuptools                       59.4.0  py310hff52083_0              conda-forge/linux-64     Cached
  + six                              1.16.0  pyh6c4a22f_0                 conda-forge/noarch       Cached
  + sqlite                           3.37.0  h9cd32fc_0                   conda-forge/linux-64     Cached
  + sysroot_linux-64                   2.12  he073ed8_15                  conda-forge/noarch        31 MB
  + tempest-remap                     2.1.1  hfda0864_0                   conda-forge/linux-64     Cached
  + tiledb                            2.3.4  he87e0bf_0                   conda-forge/linux-64     Cached
  + tk                               8.6.11  h27826a3_1                   conda-forge/linux-64     Cached
  + tornado                             6.1  py310h6acc77f_2              conda-forge/linux-64     Cached
  + traitlets                         5.1.1  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + tzcode                            2021e  h7f98852_0                   conda-forge/linux-64     Cached
  + tzdata                            2021e  he74cb21_0                   conda-forge/noarch       Cached
  + udunits2                      2.2.27.27  hc3e0081_3                   conda-forge/linux-64     Cached
  + urllib3                          1.26.7  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + wheel                            0.37.0  pyhd8ed1ab_1                 conda-forge/noarch       Cached
  + xerces-c                          3.2.3  h9d8b166_3                   conda-forge/linux-64     Cached
  + xorg-fixesproto                     5.0  h7f98852_1002                conda-forge/linux-64     Cached
  + xorg-imake                        1.0.7  0                            conda-forge/linux-64     Cached
  + xorg-inputproto                   2.3.2  h7f98852_1002                conda-forge/linux-64     Cached
  + xorg-kbproto                      1.0.7  h7f98852_1002                conda-forge/linux-64     Cached
  + xorg-libice                      1.0.10  h7f98852_0                   conda-forge/linux-64     Cached
  + xorg-libsm                        1.2.3  hd9c2040_1000                conda-forge/linux-64     Cached
  + xorg-libx11                      1.6.12  h36c2ea0_0                   conda-forge/linux-64     Cached
  + xorg-libxau                       1.0.9  h7f98852_0                   conda-forge/linux-64     Cached
  + xorg-libxaw                      1.0.14  h7f98852_0                   conda-forge/linux-64     Cached
  + xorg-libxdmcp                     1.1.3  h7f98852_0                   conda-forge/linux-64     Cached
  + xorg-libxext                      1.3.4  h516909a_0                   conda-forge/linux-64     Cached
  + xorg-libxfixes                    5.0.3  h516909a_1004                conda-forge/linux-64     Cached
  + xorg-libxi                       1.7.10  h516909a_0                   conda-forge/linux-64     Cached
  + xorg-libxmu                       1.1.3  h516909a_0                   conda-forge/linux-64     Cached
  + xorg-libxpm                      3.5.13  h516909a_0                   conda-forge/linux-64     Cached
  + xorg-libxrender                  0.9.10  h516909a_1002                conda-forge/linux-64     Cached
  + xorg-libxt                        1.1.5  h516909a_1003                conda-forge/linux-64     Cached
  + xorg-makedepend                   1.0.6  he1b5a44_1                   conda-forge/linux-64     Cached
  + xorg-renderproto                 0.11.1  h7f98852_1002                conda-forge/linux-64     Cached
  + xorg-xextproto                    7.3.0  h7f98852_1002                conda-forge/linux-64     Cached
  + xorg-xproto                      7.0.31  h7f98852_1007                conda-forge/linux-64     Cached
  + xz                                5.2.5  h516909a_1                   conda-forge/linux-64     Cached
  + zipp                              3.6.0  pyhd8ed1ab_0                 conda-forge/noarch       Cached
  + zlib                             1.2.11  h36c2ea0_1013                conda-forge/linux-64     Cached
  + zstd                              1.5.0  ha95c52a_0                   conda-forge/linux-64     Cached

  Summary:

  Install: 227 packages

  Total download: 121 MB



2b. Validate hdf libraries are correctly installed

(amipbcs211203) bash-4.2$ ls -al /home/durack1/anaconda3/envs/amipbcs211203/lib/libmfhdf.*
-rw-rw-r-- 13 durack1 climate 347030 May 19  2021 /home/durack1/anaconda3/envs/amipbcs211203/lib/libmfhdf.a
lrwxrwxrwx  1 durack1 climate     17 Dec  3  2021 /home/durack1/anaconda3/envs/amipbcs211203/lib/libmfhdf.so -> libmfhdf.so.0.0.0
lrwxrwxrwx  1 durack1 climate     17 Dec  3  2021 /home/durack1/anaconda3/envs/amipbcs211203/lib/libmfhdf.so.0 -> libmfhdf.so.0.0.0
-rwxrwxr-x 13 durack1 climate 190992 May 19  2021 /home/durack1/anaconda3/envs/amipbcs211203/lib/libmfhdf.so.0.0.0



3. Copy previous directory contents to new directory and purge redundant files (and archive previous v1.2.0 version)
(amipbcs211203) bash-4.2$ date
Fri Dec  3 09:33:02 PST 2021
(base) bash-4.2$ pwd
/work/durack1/Shared/150219_AMIPForcingData
(amipbcs211203) bash-4.2$ mv SST_1-2-0 SST_1-2-0_old4
(amipbcs211203) bash-4.2$ mkdir SST_1-2-0
(amipbcs211203) bash-4.2$ cp -R SST_1-2-0_old4/* SST_1-2-0/ ; # Copy files from previous v1.2.0_old4 pre-release
(amipbcs211203) bash-4.2$ cd SST_1-2-0
(amipbcs211203) bash-4.2$ rm -f 211104* ; # Cleanup old files
(amipbcs211203) bash-4.2$ rm -f MODEL*.nc ; # Cleanup old files



4. Edit script files (there is a need to update PATH info to include the binaries installed into the conda env)
-> changes were made in the 1.2.0 pre-release, as the path was copied no updates needed, this time
-> did need to update month info, as an additional month (August) was included from the last run
-> download date needs to be updated in the *.csh file
> download.sh
-> Many edits required for paths, start end indexing etc
> SSTICE.Update.unf.csh
-> Paths require updating lines 86-88
# Environment variables set in download.sh, removed from NCL script
# Correct all indexing for new downloaded input data, and previous version files



5. Run
(amipbcs211203) bash-4.2$ date
Mon Nov 15 15:51:24 PST 2021
(amipbcs211203) bash-4.2$ ./download.sh



6. Edit sanitize and CMORize files
-> Update all files
> ~git/amipbcs/sanitize.py
-> update all version info, input filename
> ~git/amipbcs/CMOR/drive_input4MIPs_bcs.json
> ~git/amipbcs/CMOR/drive_input4MIPs_obs.json
-> update all version info



7. Run sanitize.py
(amipbcs211203) bash-4.2$ date
Fri Dec  3 09:42:45 PST 2021
-> ensure conda env includes f2py/gfortran, and compile in same env
(amipbcs211203) bash-4.2$ cd pcmdiAmipBcs/
(amipbcs211203) bash-4.2$ ./compile 
gfortran compiler on path
--fcompiler=gnu95 GNU Fortran 95 compiler (11.2.0)
Cleanup and compile will begin in 2 seconds..
..
-> run script that will process and CMORize data
(amipbcs211203) bash-4.2$ cd ..
(amipbcs211203) bash-4.2$ python sanitize.py



8. Validate data against previously published version - edit validation script
(amipbcs211203) bash-4.2$ conda activate cdms315vcsjoblib210915
(cdms315vcsjoblib210915) bash-4.2$ pwd
/home/durack1/git/amipbcs
-> edit validation script
> make_newVsOldDiffs.py ; # Be wary of any unit changes
(cdms315vcsjoblib210915) bash-4.2$ python make_newVsOldDiffs.py > /work/durack1/Shared/150219_AMIPForcingData/211202_1544_make_newVsOldDiffs.txt



9. Validate output movies
# Pay particular attention to end points (2018-)



10. Publish data to ESGF
> Alert Sasha for this step



11. Update github repo
> cd /home/durack1/git/amipbcs
# Delete old and update data contributed (*.nc files are not hosted in git - .gitignore)
(amipbcs210923) bash-4.2$ rm -rf SST_1-1-6 ; # Purge local dir (if offical release replaces it)
(amipbcs210923) bash-4.2$ rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/SST_1-2-0 . ; # Copy remote dir
(amipbcs210923) bash-4.2$ mv SST_1-2-0 SST_1-2-0_1-1-7 ; # Update to reflect current version
(amipbcs210923) bash-4.2$ chmod 755 SST_1-2-0_1-1-7 ; # Update file permissions
(amipbcs210923) bash-4.2$ chmod 755 SST_*/*.csh SST_*/*.dat SST_*/*.f SST_*/*.ncl SST_*/*.sh
(amipbcs210923) bash-4.2$ chmod 644 SST_*/*.bz2 SST_*/*.nc SST_*/*.txt
# Update this file once everything is updated
> rsync -vrut /work/durack1/Shared/150219_AMIPForcingData/README.txt .
# Update README.md
# Add new dir to repo
> git add SST_1-2-0_1-1-7
# Commit version
> git commit -am 'Updated for new release - 1.1.7; CMOR3.6.1'
> git push
