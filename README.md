# amipbcs - AMIP dataset prepared for input4MIPs

Code to generate boundary condition data for the AMIP (Atmospheric Model Intercomparison Project) experiment as part of CMIP6 and CMIP6Plus phases and the input4MIPs project

The latest generated data can be obtained from the [PCMDI ESGF data portal](https://esgf-node.llnl.gov/search/input4mips/?source_id=PCMDI-AMIP-1-1-9)

The code depends upon the following packages:
- [**CDMS2 3.1.5**](https://github.com/CDAT/cdms) (Available through [conda](https://anaconda.org/conda-forge/cdms2/files))
- [**CMOR 3.7.1+**](https://github.com/PCMDI/cmor) (Available through [conda](https://anaconda.org/conda-forge/cmor/files))
- [**gfortran 12.2.0+**](https://gcc.gnu.org/wiki/GFortran) (Available through [conda](https://anaconda.org/conda-forge/gfortran/files))

These packages and libraries are available for linux-64 and osx-64 (x86_64) architectures

## Contributors

[![Contributors](https://contrib.rocks/image?repo=PCMDI/amipbcs)](https://github.com/PCMDI/amipbcs/graphs/contributors)

Thanks to our contributors!

## Acknowledgement

Content in this repository is developed by climate and computer scientists from the Program for Climate Model Diagnosis and Intercomparison ([PCMDI][PCMDI]) at Lawrence Livermore National Laboratory ([LLNL][LLNL]). This work is sponsored by the Regional and Global Model Analysis ([RGMA][RGMA]) program, of the Earth and Environmental Systems Sciences Division ([EESSD][EESSD]) in the Office of Biological and Environmental Research ([BER][BER]) within the [Department of Energy][DOE]'s [Office of Science][OS]. The work is performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.

<p>
    <img src="https://pcmdi.github.io/assets/PCMDI/100px-PCMDI-Logo-NoText-square-png8.png"
         width="65"
         style="margin-right: 30px"
         title="Program for Climate Model Diagnosis and Intercomparison"
         alt="Program for Climate Model Diagnosis and Intercomparison"
    >&nbsp;
    <img src="https://pcmdi.github.io/assets/DOE/480px-DOE_Seal_Color.png"
         width="65"
         style="margin-right: 30px"
         title="United States Department of Energy"
         alt="United States Department of Energy"
    >&nbsp;
    <img src="https://pcmdi.github.io/assets/LLNL/212px-LLNLiconPMS286-WHITEBACKGROUND.png"
         width="65"
         style="margin-right: 30px"
         title="Lawrence Livermore National Laboratory"
         alt="Lawrence Livermore National Laboratory"
    >
</p>
