# EcH2O-iso-nitrate: A tracer-aided ecohydrological modelling framework 

## Model introdution
EcH2O-iso-nitrate is an integrated modelling framework for simulations of hydrology, isotopes, and nitrogen. The model are fully distributed and physics-based, which consists of a hydrological model EcH2O-iso and a nitrate module.


### EcH2O-iso
EcH2O-iso was modified from Kuppel et al., (2018). It originally builds on the physically-based, spatially distributed ecohydrologic model EcH2O developed in C++ in the Regional Hydrology Lab at the University of Montana.

The specifity of EcH2O-iso is the implementation of stable
water isotopes (2H and 18O) and age tracking. It is mostly based on an immediate, complete implicit scheme for mixing
water transiting between compartments within time steps. Evaporative fractionation of isotopes is also included.

 For detailed documentation, please refer to [IGB wetsite](https://www.igb-berlin.de/en/ech2o-iso) or [Gitlab repository](https://gitlab.igb-berlin.de/ech2o-iso/ech2o-iso).
>*Kuppel, S., Tetzlaff, D., Maneta, M., & Soulsby, C. (2018). EcH2O-iso 1.0: Water isotopes and age tracking in a process-based, distributed ecohydrological model. Geoscientific Model Development, 11, 3045-3069. https://doi.org/10.5194/gmd-11-3045-2018*

### Nitrate module
Nitrate module was modified from Yang et al., (2024), which is originally based on the conceptualisation in the Hydrological Predictions for the Environment-HYPE model.

The key biogeochemical transformations are simualted for nitrogen, including denitrification, mineralisation, degradation, dissolution in terrestrail grids, as well as denitrification and assimilation in stream water.
For detailed documentation please refer to Yang et al., (2024).
>*Yang, X., Tetzlaff, D., Jin, J., Li, Q., Borchardt, D., & Soulsby, C. (2024). Linking terrestrial biogeochemical processes and water ages to catchment water quality: A new Damköhler analysis based on coupled modeling of isotope tracers and nitrate dynamics. Water Research, 262, 122118. https://doi.org/10.1016/j.watres.2024.122118*



## Data introduction
The data used for a recent modelling applications was also included (Wu et al., 2025).

The folder data_DMC contains the GIS setting of the catchment Demnitz Mill Creek, DMC, Germany, the climatic inputs, and the observations of discharge, in-stream isotopes, groundwater level, and in-stream nitrate concentrations.

>*Wu, S., Tetzlaff, D., Yang, X., Sauter, T., & Soulsby, C. Hydrological connectivity dominates NO3-N cycling in complex landscapes – insights from integration of isotopes and water quality modelling. Submitted*