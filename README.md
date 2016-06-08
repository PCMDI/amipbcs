# amipbcs - AMIP dataset prepared for input4MIPs
Code to generate boundary condition data for the AMIP (Atmospheric Model Intercomparison Project) experiment as part of CMIP6 and the input4MIPs project

Generated data can be obtained from the [PCMDI ESGF data portal](https://pcmdi.llnl.gov/search/input4mips/)

The code depends upon the following packages:
- [**UV-CDAT 2.6.0-prerelease**](https://github.com/UV-CDAT/uvcdat) (Available through [conda](https://anaconda.org/uvcdat/uvcdat/files) installation)
- [**CMOR 3**](https://github.com/PCMDI/cmor) (Available through [conda](https://anaconda.org/PCMDI/cmor/files) installation)
- [**gfortran 4.4.7+**](https://gcc.gnu.org/wiki/GFortran)

The sub-libraries of UV-CDAT 2.6.0+ are required:
- [EzGet](https://github.com/UV-CDAT/EzGet)
- [lats](https://github.com/UV-CDAT/lats)
- [libdrs](https://github.com/UV-CDAT/libdrs)

These libraries are available for linux-64 and osx-64 (x86_64) architectures
