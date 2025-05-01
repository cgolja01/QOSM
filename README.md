# QOSM

[Pronounced "cosm", it will make sense shortly]

When compelete, this will be a set of simple models of the QBO. These will couple the 
zonal momentum equation to radiation and ozone photochemistry. There will be
many such simple models (micro-QOSMs if you will) that will make
different approximations to the equations and together they will form some sort
of sensible hierarchy which we can use to understand the QBO.

This model is a fork of Edward Charlesworth SIRACHA model (Charlesworth et al., 2019).
More information about the original model can also be found in his thesis published 
[here](https://api.mountainscholar.org/server/api/core/bitstreams/e9e84e9f-2bbd-49e3-85a8-25bf547b5a46/content).

Sally Dacie has also extended the model chemistry (Dacie et al., 2019) but this fork
does not use these modifications. We make modifications to the photochemistry in the model
along similar lines to the Dacie version but they are not entirely equivalent (mostly
because I made them before becoming aware of the Dacie work).

This version of the model is as documented in Ming et al. (2025). The main changes from the
siracha version are:
- additions to the photochemistry
- the option to impose a background + pertuurbation that varies in height
- the option to impose a perturbation to the total NOx profile
- a fudge to the NO2 coefficient jno2 which really should be thought throught better.
- multiple calls to radiation and photochemistry where various quantities  (e.g., T, O3)
are held constant in turn and the heating rates and photochemical tendencies are stored.
This is used to diagnose the linearise coefficients.

The radiation code is a modified version of [RRTMG](http://rtweb.aer.com/rrtm_frame.html)
source files, so to distribute the files, we need to inform you of the license for RRTMG:

-----------------------------------------------------------------------
RRTM/RRTMG Copyright and Disclaimer

    Copyright © 2002-2010, Atmospheric and Environmental 
    Research, Inc. (AER, Inc.). This software may be used, 
    copied, or redistributed as long as it is not sold and 
    this copyright notice is reproduced on each copy made. 
    This model is provided as is without any express or 
    implied warranties. 
-----------------------------------------------------------------------

This code uses version 4.85 for the RRTM LW and 4.0 for RRTM SW component.

Dependencies:
- netCDF4
- numpy
- [pygeode](https://pygeode.github.io/index.html)

How to run the code:
- First compile RRTM LW and SW for your architecture and place the binaries
in the folder `rce_pce_code`
- Modify the `input_file.nc`
- Make any further changes to values and flags in `run_rce_pce.py`
- Run the model using `python run_rce_pce.py`
- The output is written to `output.nc`

# References

Charlesworth, E. J., Birner, T., & Albers, J. R. (2019). Ozone transport-radiation feedbacks in the tropical tropopause layer. Geophysical Research Letters, 46, 14195–14202. [doi:10.1029/2019GL084679](https://doi.org/10.1029/2019GL084679)

Dacie, S., and Coauthors, 2019: A 1D RCE Study of Factors Affecting the Tropical Tropopause Layer and Surface Climate. J. Climate, 32, 6769–6782, [doi:10.1175/JCLI-D-18-0778.1.](https://doi.org/10.1175/JCLI-D-18-0778.1.)

Ming, A., Hitchcock, P., Orbe, C., & Dubé, K. (2025). Phase and amplitude relationships between ozone, temperature, and circulation in the quasi-biennial oscillation. Journal of Geophysical Research: Atmospheres, 130, e2024JD042469. [doi:10.1029/2024JD042469](https://doi.org/10.1029/2024JD042469)




