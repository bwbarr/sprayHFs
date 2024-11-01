# Seastate-Dependent Air-Sea Heat Fluxes with Sea Spray in High Winds

README EDITS IN PROGRESS

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

The spray physics modeled in this parameterization produces changes to the surface sensible and latent heat fluxes, the surface buoyancy flux (i.e., the Obukhov length) that determines surface layer stability, and the near-surface profiles of thermodynamic variables.  We quantify these changes for use in modifying an existing bulk surface layer scheme by computing a large number of quantities representing changes to affected variables, which should be implemented as additive corrections to the same variables calculated by the existing algorith.  The variables addressed are the surface sensible and latent heat fluxes ($`H_{S,1}`$ and $`H_{L,1}`$)It is likely that not all provided values will be required to update the existing scheme.  The which satisfy the following relationships:

```math
H_{S,1} = H^{\prime}_{S,1} + dH_{S,1,spr} \: \: (1a)
```
```math
H_{L,1} = H^{\prime}_{L,1} + dH_{L,1,spr} \: \: (1b)
```
```math
\theta_* = \theta^{\prime}_* + d\theta_{*,spr} \: \: (2a)
```
```math
q_* = q^{\prime}_* + dq_{*,spr} \: \: (2b)
```
```math
\theta_{v*} = \theta^{\prime}_{v*} + d\theta_{v*,spr} \: \: (2c)
```
```math
\theta_{ref} = \theta^{\prime}_{ref} + d\theta_{ref,spr} \: \: (3a)
```
```math
T_{ref} = T^{\prime}_{ref} + dT_{ref,spr} \: \: (3b)
```
```math
q_{ref} = q^{\prime}_{ref} + dq_{ref,spr} \: \: (3c)
```
```math
s_{ref} = s^{\prime}_{ref} + ds_{ref,spr} \: \: (3d)
```

Here $`H^{\prime}_{S,1}`$ and $`H^{\prime}_{L,1}`$ are the surface sensible and latent heat fluxes without spray, respectively, $`\theta^{\prime}_*`$, $`q^{\prime}_*`$, and $`\theta^{\prime}_{v*}`$ are turbulent flux scales for sensible heat, water vapor, and buoyancy without spray, respectively, and $`\theta^{\prime}_{ref}`$, $`T^{\prime}_{ref}`$, $`q^{\prime}_{ref}`$, and $`s^{\prime}_{ref}`$ are potential temperature $\theta$, temperature $T$, specific humidity $q$, and saturation ratio $s$ (i.e. fractional relative humidity or $RH$/100%) at a reference height without spray, respectively.  Non-primed versions of these variables represent the same quantities in the presence of spray, and $`dH_{S,1,spr}`$, $`dH_{L,1,spr}`$, $`d\theta_{*,spr}`$, $`dq_{*,spr}`$, $`d\theta_{v*,spr}`$, $`d\theta_{ref,spr}`$, $`dT_{ref,spr}`$, $`dq_{ref,spr}`$, and $`ds_{ref,spr}`$ represent the changes to these variables due to spray.hich are calculated by the existing surface layer scheme.  $H_{S,spr}$, $H_{R,spr}$, and $H_{L,spr}$ are spray heat fluxes, and $\gamma_S$ and $\gamma_L$ are feedback coefficients.  We define $dH_{S,1,spr} = \gamma_S \left( H_{S,spr} - H_{R,spr} \right)$ and $dH_{L,1,spr} = \gamma_L H_{L,spr}$ as the changes to the existing bulk heat fluxes due to spray.  Equation (1a,b) can be written equivalently in terms of turbulent flux scales as follows:

```math
```
```math
```

Here $\theta* = -H_{S,1}/(\rho_a c_{p,a} u*)$ and $q* = -H_{L,1}/(\rho_a L_v u*)$ are the turbulent flux scales for potential temperature $\theta$ and specific humidity $q$ with spray, and $`\theta*^{\prime} = H^{\prime}_S/(\rho_a c_{p,a} u*)`$ and $`q*^{\prime} = H^{\prime}_L/(\rho_a L_v u*)`$ are the same quantities without spray, with $\rho_a$ as the air density, $c_{p,a}$ as the air specific heat at constant pressure, $u*$ as the friction velocity, and $L_v$ as the latent heat of vaporization of water.  $`d\theta*_{spr} = -dH_{S,1,spr}/(\rho_a c_{p,a} u*)`$ and $`dq*_{spr} = -dH_{L,1,spr}/(\rho_a L_v u*)`$ are changes to $`\theta*^{\prime}`$ and $`q*^{\prime}`$ due to spray.

A user incorporates spray heat fluxes into an existing bulk surface layer code as follows:

1. Call subroutine `sprayHFs()` directly after the existing calculation of the bulk heat fluxes $`H^{\prime}_S`$ and $`H^{\prime}_L`$ _OR_ the bulk turbulent flux scales $`\theta*^{\prime}`$ and $`q*^{\prime}`$, depending on the existing code setup.  Implementing this call will likely involve passing additional fields (particularly surface wave properties) into the existing bulk model code.  `sprayHFs()` returns $dH_{S,1,spr}$ (variable name `dHS1_spr`, units of $`W \, m^{-2}`$), $dH_{L,1,spr}$ (variable name `dHL1_spr`, units of $`W \, m^{-2}`$), $`d\theta*_{spr}`$ (variable name `dthstar_spr`, units of $K$), and $`dq*_{spr}`$ (variable name `dqstar_spr`, units of $`kg \, kg^{-1}`$.
2. Directly after the `sprayHFs()` call, add $dH_{S,1,spr}$ to $`H^{\prime}_S`$ and $dH_{L,1,spr}$ to $`H^{\prime}_L`$ _OR_ add $`d\theta*_{spr}`$ to $`\theta*^{\prime}`$ and $`dq*_{spr}`$ to $`q*^{\prime}`$, producing the total surface heat fluxes or turbulent flux scales with spray.  Note that some models carry latent heat and moisture fluxes separately -- update moisture flux too if necessary.
3. Make any changes necessary to propagate the spray-induced changes to heat fluxes or turbulent flux scales into the calculation of the Obukhov length $L$ that is used in stability calculations.  `sprayHFs()` outputs an internal calculation of the spray-modified $`\theta_v*`$ (the turbulent flux scale for virtual potential temperature that goes into $L$) for reference.  Note that `sprayHFs()` takes $L$ as an input.  The implicit relationship between $L$ and $dH_{S,1,spr}$, $dH_{L,1,spr}$, $`d\theta*_{spr}`$, and $`dq*_{spr}`$ should be treated the same way as in the existing bulk code, i.e., the $L$ passed to `sprayHFs()` should come from 1) the previous model timestep if there is no internal loop for $L$ at each model timestep or 2) the previous internal loop iteration if there is an internal loop for $L$ at each model timestep.
4. Update any calculations for diagnosed reference surface layer air temperature $T$ and specific humidity $q$ (e.g., 2-m reference values) made by the existing code.  `sprayHFs()` outputs `dtref_spr` and `dqref_spr`, which are spray-induced changes to $T$ and $q$ at a user-specified height that is prescribed using the `sprayHFs()` input parameter `z_ref`.  Add these perturbations to the diagnosed reference $T$ and $q$ values coming from the existing code.  The spray parameterization does not produce any direct changes to the subgrid wind profile (i.e., no correction is needed for the 10-m windspeed).

Subroutine `sprayHFs()` uses fields at the lowest atmospheric model mass level and at the surface, many of which are likely already used by the existing bulk algorithm code.  Additional information on input fields is given in the header of the `sprayHFs()` subroutine code in `module_sprayHFs.F90`.  `sprayHFs()` is set up for calling point-by-point, i.e., within the double DO loops over x and y horizontal gridpoints.  `sprayHFs()` automatically sets spray heat fluxes to zero if the (internally calculated) 10-m windspeed is less than a predetermined, hard-coded lower bound (currently 10 $`m \, s^{-1}`$).  The user is responsible for making sure that `sprayHFs()` is only called for water points, i.e., `sprayHFs()` has no internal check for land vs water.  In this code, heat fluxes from the ocean to the atmosphere are defined as positive, and momentum fluxes (stress) from the atmosphere to the ocean are defined as positive.

The user selects the spray generation function to use in `sprayHFs()` by providing one of the string keys in Section 1 for the input variable `whichSSGF`.  If using the seastate-based model `BCF23_Seastate`, all wave parameters are required (`eps`, `dcp`, `swh`, `mss`).  If using a wind-based model (`F94_MOM80` or `F94_BCF23`), `swh` is still required (for calculating the droplet settling timescale), but the remaining wave parameters (`eps`, `dcp`, `mss`) are not used (the user should pass dummy values, e.g., zeros or NaNs).

Subroutine `sprayHFs()` makes internal calculations for several variables that may already appear in the existing bulk algorithm code, including $`H^{\prime}_S`$ (variable name `H_S0pr`), $`H^{\prime}_L`$ (`H_L0pr`), $`\theta*^{\prime}`$ (`thstar_pr`), $`q*^{\prime}`$ (`qstar_pr`), $`u*`$ (`ustar`), and the surface stress $\tau$ (`tau`).  `sprayHFs()` returns `tau`, `H_S0pr`, and `H_L0pr` by default in case the user wants to compare.  These internal values may not exactly match those in the existing bulk code if, for instance, different stability functions are used.  This is probably OK and should not affect the results very much.  Please contact Ben Barr if you are interested in discussing ways to reconcile any discrepancies.

`sprayHFs()` returns 21 `INTENT(OUT)` variables by default.  Some of these are required to implement the spray physics into the existing code (i.e., `dHS1_spr` and `dHL1_spr` _OR_ `dthstar_spr` and `dqstar_spr`, as well as `dtref_spr` and `dqref_spr` to update diagnosed reference conditions), but many are included only for debugging and interpretation of the model physics.  The user is free to suppress any unneeded `INTENT(OUT)` fields.

## References

The following publications/references provide additional information on the spray model, its implementation, and scientific findings coming from its use in coupled model simulations:

+ Barr, B. W. and S. S. Chen: Impacts of seastate-dependent sea spray heat fluxes on tropical cyclone structure and intensity in fully coupled atmosphere-wave-ocean model simulations. _Submitted to J. Adv. Model. Earth Systems_.

+ Barr, B. W., 2023: Seastate-dependent sea spray heat fluxes and impacts on tropical cyclone structure and intensity using fully coupled atmosphere-wave-ocean model simulations. Ph.D. Dissertation, University of Washington, Seattle, WA.

+ Barr, B. W., S. S. Chen, and C. W. Fairall, 2023: Sea-state-dependent sea spray and air-sea heat fluxes in tropical cyclones: A new parameterization for fully coupled atmosphere-wave-ocean models. _J. Atmos. Sci._, **80**, 933 - 960, https://doi.org/10.1175/JAS-D-22-0126.1.

The following additional publications are referenced in the documentation above:

+ Fairall, C. W., J. D. Kepert, and G. J. Holland, 1994: The effect of sea spray on surface energy transports over the ocean. _Glob. Atmos. Ocean System_, **2**, 121 - 142.

+ Monahan, E. C. and I. O'Muircheartaigh, 1980: Optimal power-law description of oceanic whitecap coverage dependence on wind speed. _J. Phys. Oceanogr._, **10**, 2094 - 2099, https://doi.org/10.1175/1520-0485(1980)010<2094:OPLDOO>2.0.CO;2.

+ Mueller, J. A. and F. Veron, 2014: Impact of sea spray on air-sea fluxes. Part II: feedback effects. _J. Phys. Oceanogr._, **44**, 2835 - 2853, https://doi.org/10.1175/JPO-D-13-0246.1.
