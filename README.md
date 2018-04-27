# amipbcs - AMIP dataset prepared for input4MIPs
Code to generate boundary condition data for the AMIP (Atmospheric Model Intercomparison Project) experiment as part of CMIP6 and the input4MIPs project

Generated data can be obtained from the [PCMDI ESGF data portal](https://esgf-node.llnl.gov/search/input4mips/?dataset_version_number=1.1.4)

The code depends upon the following packages:
- [**CDAT 8.0 with cd77**](https://github.com/UV-CDAT/uvcdat) (Available through [conda](https://anaconda.org/uvcdat/uvcdat/files))
- [**CMOR 3.3.2+**](https://github.com/PCMDI/cmor) (Available through [conda](https://anaconda.org/PCMDI/cmor/files))
- [**gfortran 4.8.5+**](https://gcc.gnu.org/wiki/GFortran)

The optional libraries of CDAT 8.0 are required:
- [EzGet](https://github.com/UV-CDAT/EzGet)
- [lats](https://github.com/UV-CDAT/lats)
- [libdrs](https://github.com/UV-CDAT/libdrs)

These packages and libraries are available for linux-64 and osx-64 (x86_64) architectures
