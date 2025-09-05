# Seastate-Dependent Air-Sea Heat Fluxes with Sea Spray in High Winds

This repository contains a Fortran module with subroutines to incorporate seastate-dependent sea spray heat flux physics into air-sea flux calculations in numerical forecast models.  Parameterization of sea spray generation, which is required for the spray heat flux physics, is also included.  The model originates from Barr, Chen, and Fairall (2023) (hereafter BCF23) and is in ongoing active development.

The code in this repository is designed to implement spray heat flux physics into an existing Fortran bulk surface heat flux algorithm with minimal disruption to the existing bulk code.  All required subroutines are contained in `module_sprayHFs.F90`.  To implement the code, a user will download `module_sprayHFs.F90` into their existing code base, modify the existing surface layer flux calculations to add a call to the top-level subroutine `sprayHFs()` and to manage its input/output (Section 2 below), and compile `module_sprayHFs.F90` with the existing model code.  The spray module is designed for use "right out of the box" with minimal changes to the module file itself by the user.

Please send any questions about model physics or this module's implementation to Ben Barr at benjamin.barr@whoi.edu.

Spray generation physics and subroutines are discussed in Section 1 below, and spray heat flux physics and subroutines are discussed in Section 2.  Other spray interactions including impacts on momentum fluxes, turbulence, and cloud processes are not addressed at this time.  

## 1. Sea Spray Generation

Breaking waves eject sea spray droplets into the air, and these droplets participate in physical processes that affect air-sea exchange in high winds.  There are many parameterizations for spray generation, with numerous unresolved issues.  We include two main classes of spray generation models herein, with several options for each class that are detailed in the header of the subroutine `sprayHFs()`.  The two classes are:

+ _Dissipation-Ejection-Based_ Spray Generation: Represents the processes of 1) droplet formation by fragmentation due to turbulent dissipation in breaking wave crests and 2) droplet ejection/entrainment by turbulent gusts over waves.  This model originates from Fairall et al. (2009) with current form described in BCF23.  Spray generation is parameterized by seastate, but a wind-only version is in development.
+ _Whitecap-Based_ Spray Generation: Represents spray generation by empirically defining the droplet flux from a unit surface whitecap area and scaling this by the whitecap fraction of the ocean surface.  This model originates from Fairall et al. (1994).  All options for this model are parameterized by winds only.

The spray generation model is selected by providing a string key to the variable `whichSSGF` in the subroutine `sprayHFs`.  For the dissipation-ejection (DE) generation model, the option `dissejec_SS_BCF23` is recommended.  This is the published form in BCF23.  All other DE options, including wind-only options, are in development and not yet ready for public use.  For the whitecap-based model, the option `whitecap_Wi_F94_BCF23_published` is recommended, as this is also the published version in BCF23.  We recommend using the DE model to directly address wave-based physical processes related to spray generation, although there are not yet sufficient observations to determine whether any given generation model is more quantitatively correct.  Comparison of spray generation models across diverse air-sea-wave conditions in your modeling system is encouraged!

## 2. Air-Sea Heat Fluxes with Spray

Ejected sea spray droplets warm near-surface air through direct heating by relatively warmer droplets, and they simulteneously moisten and cool near-surface air through enthalpy-neutral evaporation.  This produces changes to the surface sensible and latent heat fluxes (or, equivalently, the turbulent heat and moisture flux scales), the Obukhov length that determines surface layer stability, and the near-surface profiles of thermodynamic variables.  The subroutine `sprayHFs()` calculates these changes and provides them as additive corrections to variables calculated by the existing bulk surface layer code.  Many corrections are provided to accommodate the variety of approaches found in bulk surface flux models, and it is likely that not all provided variables are required to update any particular scheme.  Corrections are provided for the surface sensible and latent heat fluxes ($`H_{S,1}`$ and $`H_{L,1}`$ respectively), the turbulent flux scales for potential temperature, specific humidity, and virtual potential temperature ($`\theta_*`$, $`q_*`$, and $`\theta_{v*}`$ respectively), and values of potential temperature, temperature, specific humidity, and saturation ratio (i.e., fractional relative humidity or RH/100%) at a user-specified reference height `z_ref` ($`\theta_{ref}`$, $`T_{ref}`$, $`q_{ref}`$, and $`s_{ref}`$ respectively).  The changes (additive corrections) are defined below, with primed variables representing values without spray (i.e., those calculated by the existing bulk code), non-primed variables representing values with spray, and variables with the prefix $d$ and subscript $spr$ (e.g., $`dH_{S,1,spr}`$) representing the additive corrections themselves.  The corrections have corresponding variable names in the code of `dHS1_spr`, `dHL1_spr`, `dthstar_spr`, `dqstar_spr`, `dthvstar_spr`, `dthref_spr`, `dtref_spr`, `dqref_spr`, and `dsref_spr`.  

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

Subroutine `sprayHFs()` uses fields at the lowest atmospheric model mass level and at the surface, and it also requires several wave variables if using a seastate-dependent DE spray generation option.  The significant wave height is required even if using a wind-only spray generation option, but this may be parameterized internally using the `param_delspr_Wi` input variable.  The subroutine also requires inputs of momentum, heat, and moisture roughness lengths, the Obukhov stability length, and an optional gust factor for surface winds.  These are typically available from the existing code and are required as inputs to `sprayHFs()` to ensure consistency with the existing bulk model.  All inputs are described in the header of `sprayHFs()`.  A user incorporates the spray heat flux physics into their bulk surface layer code by calling `sprayHFs()` after all required inputs are determined.  Then, the returned corrections may be applied after each relevant bulk variable is calculated.  Due to the variety in code structure among bulk models, the user must be careful to understand their specific case and update surface heat fluxes, turbulent flux scales, the Obukhov length, and reference height variables correctly.  Note that some models carry latent heat and moisture fluxes separately, and both should be updated for spray.  Note also that the spray parameterization does not produce any direct changes to the subgrid wind profile, so no correction is needed for the 10-m windspeed.

`sprayHFs()` is set up for calling point-by-point, i.e., within the double DO loops over x and y horizontal gridpoints.  `sprayHFs()` automatically sets spray heat fluxes to zero if the (internally calculated) 10-m windspeed is less than a predetermined, hard-coded lower bound (currently 10 $`m \, s^{-1}`$).  The user is responsible for making sure that `sprayHFs()` is only called for water points, i.e., `sprayHFs()` has no internal check for land vs water.  In this code, heat fluxes from the ocean to the atmosphere are defined as positive, and momentum fluxes (stress) from the atmosphere to the ocean are defined as positive.

Subroutine `sprayHFs()` calculates and returns the bulk-only heat fluxes $`H^{\prime}_{S,1}`$ (variable name `H_S0pr`) and $`H^{\prime}_{L,1}`$ (`H_L0pr`) and the surface stress $\tau$ (`tau`), which are provided in case the user wants to compare to their existing bulk calculations.  These values may not exactly match those from the existing bulk code if, for instance, different stability functions are used.  This is probably OK and should not affect the results very much.  Please contact Ben Barr if you are interested in discussing or reconciling any discrepancies.  `sprayHFs()` also returns the spray mass flux `M_spr`, spray heat fluxes `H_Tspr`, `H_Rspr`, `H_SNspr`, and `H_Lspr`, and spray feedback coefficients `alpha_S`, `beta_S`, `beta_L`, `gamma_S`, and `gamma_L`, which are useful diagnostics for debugging, analysis, and interpretation.

## References

The following publications/references provide additional information on the spray model, its implementation, and scientific findings coming from its use in coupled model simulations:

+ Barr, B. W. and S. S. Chen: Impacts of seastate-dependent sea spray heat fluxes on tropical cyclone structure and intensity in fully coupled atmosphere-wave-ocean model simulations. _J. Adv. Model. Earth Syst._, **17**, e2024MS004550, https://doi.org/10.1029/2024MS004550.

+ Barr, B. W., 2023: Seastate-dependent sea spray heat fluxes and impacts on tropical cyclone structure and intensity using fully coupled atmosphere-wave-ocean model simulations. Ph.D. Dissertation, University of Washington, Seattle, WA.

+ Barr, B. W., S. S. Chen, and C. W. Fairall, 2023: Sea-state-dependent sea spray and air-sea heat fluxes in tropical cyclones: A new parameterization for fully coupled atmosphere-wave-ocean models. _J. Atmos. Sci._, **80**, 933 - 960, https://doi.org/10.1175/JAS-D-22-0126.1.

The following additional publications are referenced in the documentation above:

+ Fairall, C. W., M. L. Banner, W. L. Peirson, W. Asher, and R. P. Morison, 2009: Investigation of the physical scaling of sea spray spume droplet production. _J. Geophys. Res._, **114**, C10001, https://doi.org/10.1029/2008JC004918.

+ Fairall, C. W., J. D. Kepert, and G. J. Holland, 1994: The effect of sea spray on surface energy transports over the ocean. _Glob. Atmos. Ocean System_, **2**, 121 - 142.

