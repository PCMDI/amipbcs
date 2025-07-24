# amipbcs - AMIP dataset prepared for input4MIPs

Code to generate sea surface temperature (SST/tos) and sea ice concentration (siconc) boundary condition data for the AMIP (Atmospheric Model Intercomparison Project) experiment as part of CMIP6, CMIP6Plus and CMIP7 phases and the input4MIPs project.

The latest generated data can be obtained from the [ESGF2-US data portal](https://esgf-node.ornl.gov/search?project=input4MIPs&activeFacets={%22mip_era%22:%22CMIP7%22})

The code depends upon the following packages:
- [**xCDAT 0.9.1+**](https://github.com/xCDAT/xcdat) (Available through [conda-forge](https://anaconda.org/conda-forge/xcdat/files))
- [**CMOR 3.11.0+**](https://github.com/PCMDI/cmor) (Available through [conda-forge](https://anaconda.org/conda-forge/cmor/files))
- [**gfortran 14.2.0+**](https://gcc.gnu.org/wiki/GFortran) (Available through [conda-forge](https://anaconda.org/conda-forge/gfortran/files))

These packages and libraries are available for linux-64, osx-64 and osx-arm64 (x86_64) architectures

## Contributors

[![Contributors](https://contrib.rocks/image?repo=PCMDI/amipbcs)](https://github.com/PCMDI/amipbcs/graphs/contributors)

Thanks to our contributors!

## Acknowledgement

Content in this repository is developed by climate and computer scientists from the Program for Climate Model Diagnosis and Intercomparison ([PCMDI][PCMDI]) at Lawrence Livermore National Laboratory ([LLNL][LLNL]). This work is sponsored by the Regional and Global Model Analysis ([RGMA][RGMA]) program, of the Earth and Environmental Systems Modeling Division ([EESM][EESM]) in the Office of Biological and Environmental Research ([BER][BER]) within the [Department of Energy][DOE]'s [Office of Science][OS]. The work is performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.

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


[PCMDI]: https://pcmdi.llnl.gov
[LLNL]: https://www.llnl.gov
[RGMA]: https://eesm.science.energy.gov/program-area/regional-global-model-analysis
[EESM]: https://eesm.science.energy.gov
[BER]: https://www.energy.gov/science/ber/biological-and-environmental-research
[DOE]: https://www.energy.gov
[OS]: https://www.energy.gov/science/office-science
