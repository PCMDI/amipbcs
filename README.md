# amipbcs
Code to generate boundary condition data for the AMIP experiment suite

Generated data can be obtained from the [PCMDI ESGF data portal](https://pcmdi.llnl.gov/search/input4mips/)

The fortran code included in the repo depends upon the following packages:
- [UV-CDAT 2.6.0-prerelease](https://github.com/UV-CDAT/uvcdat)
- [gfortran 4.4.7+](https://gcc.gnu.org/wiki/GFortran)

The following required libraries are included in UV-CDAT 2.6.0+
- [EzGet](https://github.com/UV-CDAT/EzGet)
- [lats](https://github.com/UV-CDAT/lats)
- [libdrs](https://github.com/UV-CDAT/libdrs)

These libraries are currently available on linux-64 and osx-64 (x86_64) architectures
