# amipbcs - AMIP dataset prepared for input4MIPs
Code to generate boundary condition data for the AMIP (Atmospheric Model Intercomparison Project) experiment as part of CMIP6 and the input4MIPs project

Generated data can be obtained from the [PCMDI ESGF data portal](https://esgf-node.llnl.gov/search/input4mips/?source_version=1.1.6)

The code depends upon the following packages:
- [**CDAT 8.1+ with cd77**](https://github.com/CDAT/cdat) (Available through [conda](https://anaconda.org/CDAT/repo))
- [**CMOR 3.5.0+**](https://github.com/PCMDI/cmor) (Available through [conda](https://anaconda.org/conda-forge/cmor/files))
- [**gfortran 7.3+**](https://gcc.gnu.org/wiki/GFortran)

The optional libraries of CDAT 8.1+ are required (these are bundled with cd77):
- [EzGet](https://github.com/UV-CDAT/EzGet)
- [lats](https://github.com/UV-CDAT/lats)
- [libdrs](https://github.com/UV-CDAT/libdrs)

These packages and libraries are available for linux-64 and osx-64 (x86_64) architectures
