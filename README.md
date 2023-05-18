# amipbcs - AMIP dataset prepared for input4MIPs
Code to generate boundary condition data for the AMIP (Atmospheric Model Intercomparison Project) experiment as part of CMIP6 and CMIP6Plus phases and the input4MIPs project

The latest generated data can be obtained from the [PCMDI ESGF data portal](https://esgf-node.llnl.gov/search/input4mips/?source_id=PCMDI-AMIP-1-1-9)

The code depends upon the following packages:
- [**CDMS2 3.1.5**](https://github.com/CDAT/cdms) (Available through [conda](https://anaconda.org/conda-forge/cdms2/files))
- [**CMOR 3.7.1+**](https://github.com/PCMDI/cmor) (Available through [conda](https://anaconda.org/conda-forge/cmor/files))
- [**gfortran 12.2.0+**](https://gcc.gnu.org/wiki/GFortran) (Available through [conda](https://anaconda.org/conda-forge/gfortran/files))

These packages and libraries are available for linux-64 and osx-64 (x86_64) architectures
