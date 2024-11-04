# Seastate-Dependent Air-Sea Heat Fluxes with Sea Spray in High Winds

This repository contains subroutines to incorporate seastate-dependent sea spray heat flux physics into an existing bulk surface layer scheme in a regional or global Earth system model.  Parameterization of sea spray generation, which is required for the spray heat flux physics, is also included.  The model implemented herein originates from Barr, Chen, and Fairall (2023) (hereafter BCF23) and is actively being improved and updated with contributions from many collaborators addressing both model physics and multiscale interactions in model simulations.

The code in this repository is designed to implement spray heat flux physics into an existing Fortran bulk heat flux algorithm with minimal disruption to the existing bulk code.  All subroutines required to implement the spray model are contained in the Fortran module file `module_sprayHFs.F90`.  To use this code, a user will download `module_sprayHFs.F90` into their existing code base, modify their existing surface flux calculations to add a call to the top-level subroutine `sprayHFs()` and to manage its input/output (discussed fully in Section 2 below), and compile `module_sprayHFs.F90` with the existing model code.  Other than minor changes to integrate `module_sprayHFs.F90` into the existing code base (suppressing optional output, etc), the user should not have to change any code in `module_sprayHFs.F90` to implement it (i.e., it is ready to use "out of the box").

Please send any questions about model physics or this module's implementation to Ben Barr at benjamin.barr@whoi.edu.

Spray generation physics and subroutines are discussed in Section 1 below, and spray heat flux physics and subroutines are discussed in Section 2.  Other spray interactions including impacts on momentum fluxes and turbulence, as well as participation in cloud processes, are not addressed at this time.  

## 1. Sea Spray Generation

Breaking waves eject sea spray droplets into the air, and these droplets participate in physical processes that affect air-sea exchange in high winds.  There are many parameterizations for spray generation, with numerous unresolved issues.  We include several parameterizations for spray generation, based both on wind and on seastate.  We recommend using the seastate-based model from BCF23 (option `BCF23_Seastate` below) to directly address wave-based physical processes related to spray production, although there are not yet sufficient observations to determine whether any given wind or seastate based model is more quantitatively correct under a given set of environmental conditions.  Comparison of spray generation models across diverse air-sea-wave conditions in your modeling system is encouraged!  The available spray generation models are:

+ `BCF23_Seastate`: A seastate-dependent spray generation model presented in BCF23 Eq. 1.  In this model, both the total spray mass flux and droplet size distribution change with wind-wave conditions.
+ `F94_MOM80`: The widely-used Fairall et al. (1994) (hereafter F94) wind-dependent source function, with the droplet size distribution per Mueller and Veron (2014) and the original whitecap fraction per Monahan and O'Muircheartaigh (1980).
+ `F94_BCF23`: An update to the F94 wind-dependent source function given as BCF23 Eq. 3, with the droplet size distribution per Mueller and Veron (2014) and the whitecap fraction per BCF23 Eq. A2.

The spray generation model is selected by providing one of the string keys above to the subroutine `sprayHFs()`, as discussed in Section 2.

## 2. Air-Sea Heat Fluxes with Spray

The spray physics included in this parameterization produce changes to the surface sensible and latent heat fluxes, the Obukhov length determining surface layer stability, and the near-surface profiles of thermodynamic variables.  We quantify these changes for use in modifying an existing bulk surface layer scheme by computing a large number of quantities representing changes to affected variables, which are implemented as additive corrections to the same variables calculated by the existing algorith.  It is likely that not all provided variables will be required to update the existing scheme.  The variables included are the surface sensible and latent heat fluxes ($`H_{S,1}`$ and $`H_{L,1}`$ respectively), the turbulent flux scales for potential temperature, specific humidity, and virtual potential temperature ($`\theta_*`$, $`q_*`$, and $`\theta_{v*}`$ respectively), and values of potential temperature, temperature, specific humidity, and saturation ratio (i.e., fractional relative humidity or RH/100%) at a user-specified reference height ($`\theta_{ref}`$, $`T_{ref}`$, $`q_{ref}`$, and $`s_{ref}`$ respectively).  The changes (additive corrections) are defined below, with primed variables representing values without spray (i.e., those calculated by the existing bulk code), non-primed variables representing values with spray, and variables with the prefix $d$ and subscript $`_{spr}`$ (e.g., $`dH_{S,1,spr}`$) representing the additive corrections themselves.

```math
dH_{S,1,spr} = H_{S,1} - H^{\prime}_{S,1} \: \: (1a)
```
```math
dH_{L,1,spr} = H_{L,1} - H^{\prime}_{L,1} \: \: (1b)
```
```math
d\theta_{*,spr} = \theta_* - \theta^{\prime}_* \: \: (2a)
```
```math
dq_{*,spr} = q_* - q^{\prime}_* \: \: (2b)
```
```math
d\theta_{v*,spr} = \theta_{v*} - \theta^{\prime}_{v*} \: \: (2c)
```
```math
d\theta_{ref,spr} = \theta_{ref} - \theta^{\prime}_{ref} \: \: (3a)
```
```math
dT_{ref,spr} = T_{ref} - T^{\prime}_{ref} \: \: (3b)
```
```math
dq_{ref,spr} = q_{ref} - q^{\prime}_{ref} \: \: (3c)
```
```math
ds_{ref,spr} = s_{ref} - s^{\prime}_{ref} \: \: (3d)
```

A user incorporates spray heat flux physics into an existing bulk surface layer code by calling subroutine `sprayHFs()` within the existing code after calculation of the relevant primed quantities above.  This will typically occur near the end of the existing subroutine, unless there is an internal loop for calculating stability, in which case `sprayHFs()` should be called within this loop and likely near its end.  Implementing the call to `sprayHFs()` will likely involve passing additional fields (particularly surface wave properties) into the existing bulk model code.  The input parameter `z_ref` allows the user to specify the height for calculating corrections to reference values (e.g., 2m).  `sprayHFs()` returns the additive corrections listed above with the following variable names: `dHS1_spr`, `dHL1_spr`, `dthstar_spr`, `dqstar_spr`, `dthvstar_spr`, `dthref_spr`, `dtref_spr`, `dqref_spr`, and `dsref_spr`.  After calling `sprayHFs()`, a user adds the additive corrections to the relevant variables from the existing scheme.  In doing so, please note the following:
1. Sensible and latent heat fluxes are sometimes quantified using $`\theta_*`$ and $`q_*`$ respectively, so care must be taken to update these in addition to or instead of $`H_{S,1}`$ and $`H_{L,1}`$ as applicable.
2. Depending on the algorithm, the Obukhov length in the existing code may be calculated in terms of surface fluxes ($`H_{S,1}`$ and $`H_{L,1}`$), separate heat and moisture turbulent flux scales ($`\theta_*`$ and $`q_*`$), or the turbulent buoyancy flux scale ($`\theta_{v*}`$).  This should be updated consistently using the provided variables.  Note that `sprayHFs()` takes $L$ as an input.  The implicit relationship between $L$ and surface fluxes should be treated the same way in `sprayHFs()`as in the existing bulk code, i.e., the $L$ passed to `sprayHFs()` should come from 1) the previous model timestep if there is no internal loop for $L$ at each model timestep or 2) the previous internal loop iteration if there is an internal loop for $L$ at each model timestep.
3. Some models carry latent heat and moisture fluxes separately, so the moisture flux may need to be updated as well.
4. The spray parameterization does not produce any direct changes to the subgrid wind profile, so no correction is needed for the 10-m windspeed.

Subroutine `sprayHFs()` uses fields at the lowest atmospheric model mass level and at the surface, many of which are likely already used by the existing bulk algorithm code.  Additional information on input fields is given in the header of the `sprayHFs()` subroutine code in `module_sprayHFs.F90`.  `sprayHFs()` is set up for calling point-by-point, i.e., within the double DO loops over x and y horizontal gridpoints.  `sprayHFs()` automatically sets spray heat fluxes to zero if the (internally calculated) 10-m windspeed is less than a predetermined, hard-coded lower bound (currently 10 $`m \, s^{-1}`$).  The user is responsible for making sure that `sprayHFs()` is only called for water points, i.e., `sprayHFs()` has no internal check for land vs water.  In this code, heat fluxes from the ocean to the atmosphere are defined as positive, and momentum fluxes (stress) from the atmosphere to the ocean are defined as positive.

The user selects the spray generation function to use in `sprayHFs()` by providing one of the string keys in Section 1 for the input variable `whichSSGF`.  If using the seastate-based model `BCF23_Seastate`, all wave parameters are required (`eps`, `dcp`, `swh`, `mss`).  If using a wind-based model (`F94_MOM80` or `F94_BCF23`), `swh` is still required (for calculating the droplet settling timescale), but the remaining wave parameters (`eps`, `dcp`, `mss`) are not used (the user should pass dummy values, e.g., zeros).  Internal parameterizations for wave physics are in development to allow cautious use of this model in atmosphere-only situations.  The user may parameterize wave processes by setting the input variable `paramWaves` to `.TRUE.`.  In this case dummy values should be passed for all wave parameters.

Subroutine `sprayHFs()` makes internal calculations for several variables that may already appear in the existing bulk algorithm code, including $`H^{\prime}_{S,1}`$ (variable name `H_S0pr`), $`H^{\prime}_{L,1}`$ (`H_L0pr`), $`\theta^{\prime}_*`$ (`thstar_pr`), $`q^{\prime}_*`$ (`qstar_pr`), friction velocity $`u_*`$ (`ustar`), and the surface stress $\tau$ (`tau`).  `sprayHFs()` returns `tau`, `H_S0pr`, and `H_L0pr` by default in case the user wants to compare.  These internal values may not exactly match those in the existing bulk code if, for instance, different stability functions are used.  This is probably OK and should not affect the results very much.  Please contact Ben Barr if you are interested in discussing ways to reconcile any discrepancies.

`sprayHFs()` returns many `INTENT(OUT)` variables by default.  Some of these are required to implement the spray physics into the existing code as discussed above, but many are included only for debugging and interpretation of the model physics.  The user is free to suppress any unneeded `INTENT(OUT)` fields by removing the `INTENT(OUT)` tag from the variable definition and removing the variable name from the subroutine definition.

## References

The following publications/references provide additional information on the spray model, its implementation, and scientific findings coming from its use in coupled model simulations:

+ Barr, B. W. and S. S. Chen: Impacts of seastate-dependent sea spray heat fluxes on tropical cyclone structure and intensity in fully coupled atmosphere-wave-ocean model simulations. _Submitted to J. Adv. Model. Earth Systems_.

+ Barr, B. W., 2023: Seastate-dependent sea spray heat fluxes and impacts on tropical cyclone structure and intensity using fully coupled atmosphere-wave-ocean model simulations. Ph.D. Dissertation, University of Washington, Seattle, WA.

+ Barr, B. W., S. S. Chen, and C. W. Fairall, 2023: Sea-state-dependent sea spray and air-sea heat fluxes in tropical cyclones: A new parameterization for fully coupled atmosphere-wave-ocean models. _J. Atmos. Sci._, **80**, 933 - 960, https://doi.org/10.1175/JAS-D-22-0126.1.

The following additional publications are referenced in the documentation above:

+ Fairall, C. W., J. D. Kepert, and G. J. Holland, 1994: The effect of sea spray on surface energy transports over the ocean. _Glob. Atmos. Ocean System_, **2**, 121 - 142.

+ Monahan, E. C. and I. O'Muircheartaigh, 1980: Optimal power-law description of oceanic whitecap coverage dependence on wind speed. _J. Phys. Oceanogr._, **10**, 2094 - 2099, https://doi.org/10.1175/1520-0485(1980)010<2094:OPLDOO>2.0.CO;2.

+ Mueller, J. A. and F. Veron, 2014: Impact of sea spray on air-sea fluxes. Part II: feedback effects. _J. Phys. Oceanogr._, **44**, 2835 - 2853, https://doi.org/10.1175/JPO-D-13-0246.1.
