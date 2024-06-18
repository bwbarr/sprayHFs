# Seastate-Dependent Air-Sea Heat Fluxes with Sea Spray in High Winds

This repository contains subroutines to incorporate seastate-dependent sea spray heat flux physics into an existing bulk surface layer scheme in a coupled regional or global Earth system model.  Parameterization of sea spray generation, which is required for the spray heat flux physics, is also included.  The model implemented herein comes from Barr et al. (2023) (hereafter BCF23).  This model was originally developed by Ben Barr, Shuyi Chen, and Chris Fairall, and it is being actively updated and improved with additional contributions from others including Hyodae Seo, Cesar Sauvage, Jim Edson, and Carol Anne Clayson.

The code in this repository was designed to implement spray heat flux physics into an existing Fortran bulk heat flux algorithm with minimal disruption to the existing bulk code.  All subroutines required to implement the spray model are contained in the Fortran module file `module_sprayHFs.F`.  To use the code, a user will download `module_sprayHFs.F` into their existing code base, modify their existing surface flux calculations to add a call to the top-level subroutine `sprayHFs()` and to manage its input/output (discussed fully in Section 2 below), and compile `module_sprayHFs.F` with the existing model code.  Other than minor changes to integrate `module_sprayHFs.F` into the existing code base (suppressing optional output, etc), the user should not have to change any code in `module_sprayHFs.F` to implement it (i.e., it is ready to use "out of the box").

Please send any questions about model physics or this module's implementation to Ben Barr at benjamin.barr@whoi.edu.

Spray generation physics and subroutines are discussed in Section 1 below, and spray heat flux physics and subroutines are discussed in Section 2.

## 1. Sea Spray Generation

Breaking waves eject sea spray droplets into the air, and these droplets participate in physical processes that affect air-sea exchange in high winds.  There are many parameterizations for spray generation, with numerous unresolved issues.  We include several parameterizations for spray generation, based both on wind and on seastate.  We recommend using the seastate-based model from BCF23 (option `BCF23_Seastate` below) to directly address physical processes related to spray production, although there are not yet sufficient observations to determine whether any given wind or seastate based model is more quantitatively correct under a given set of environmental conditions.  Comparison of spray generation models across diverse air-sea-wave conditions in your modeling system is encouraged!  The available spray generation models are:

+ `BCF23_Seastate`: A seastate-dependent spray generation model presented in BCF23 Eq. 1.  In this model, both the total spray mass flux and droplet size distribution change with wind-wave conditions.
+ `F94_MOM80`: The widely-used Fairall et al. (1994) (hereafter F94) wind-dependent source function, with the droplet size distribution per Mueller and Veron (2014) and the original whitecap fraction per Monahan and O'Muircheartaigh (1980).
+ `F94_BCF23`: An update to the F94 wind-dependent source function given as BCF23 Eq. 3, with the droplet size distribution per Mueller and Veron (2014) and the whitecap fraction per BCF23 Eq. A2.

The spray generation model is selected by providing one of the string keys above to the subroutine `sprayHFs()`, as discussed in Section 2.

## 2. Air-Sea Heat Fluxes with Spray

Spray heat fluxes change the total surface sensible and latent heat fluxes and the buoyancy flux (i.e., the Obukhov length) used for stability calculations.  Spray heat flux physics are implemented according to BCF23.  Per BCF23 Eq. 16a and 16b, the total surface sensible and latent heat fluxes with spray, $H_{S,1}$ and $H_{L,1}$ respectively, are

```math
H_{S,1} = H^{\prime}_S + \gamma_S \left( H_{S,spr} - H_{R,spr} \right) = H^{\prime}_S + dH_{S,1,spr}
H_{L,1} = H^{\prime}_L + \gamma_L H_{L,spr} = H^{\prime}_L + dH_{L,1,spr}
```
Here $`H^{\prime}_S`$ and $`H^{\prime}_L`$ are the bulk sensible and latent heat fluxes without spray, which are calculated by the existing surface layer scheme.  $H_{S,spr}$, $H_{R,spr}$, and $H_{L,spr}$ are spray heat fluxes, and $\gamma_S$ and $\gamma_L$ are feedback coefficients.  We define $dH_{S,1,spr} = \gamma_S \left( H_{S,spr} - H_{R,spr} \right)$ and $dH_{L,1,spr} = \gamma_L H_{L,spr}$ as the changes to the existing bulk heat fluxes due to spray.

!     A user incorporates spray heat fluxes into an existing bulk surface layer 
! code by calling subroutine sprayHFs() directly after the existing calculation 
! of the bulk heat fluxes H_S0pr and H_L0pr.  sprayHFs() returns dHS1spr and 
! dHL1spr.  Then, the user adds dHS1spr to H_S0pr and dHL1spr to H_L0pr.  Note 
! that some models carry latent heat and moisture fluxes separately -- update 
! moisture flux here too if necessary.  Finally, if the existing code does not 
! calculate the Obukhov length L directly from the modified H_S0pr and H_L0pr 
! (for instance, L may be calculated from a "theta v star" (variable name 
! thvstar in this code) that is computed separately from H_S0pr and H_L0pr), the
! user must also make sure that dHS1spr and dHL1spr are incorporated correctly 
! into the calculation of L.  sprayHFs() outputs a calculation of thvstar for 
! reference.  Note that sprayHFs() takes L as an input, just like a standard 
! bulk flux algorithm.  The implicit relationship between L and the spray heat 
! fluxes should be treated the same as in the existing bulk model, i.e., the L 
! passed to sprayHFs() should come from the previous timestep or from the 
! previous iteration if there is an internal loop for L.
!     Subroutine sprayHFs() uses fields at the lowest model mass level and at 
! the surface, which should be the same fields that are used by the existing 
! bulk algorithm code.  Additional information on input fields can be found in 
! the header of the sprayHFs() subroutine code.  sprayHFs() is set up for
! calling point-by-point, i.e., within the double DO loops over x and y 
! horizontal gridpoints.  sprayHFs() automatically sets spray heat fluxes to 
! zero if the (internally calculated) 10-m windspeed is less than a 
! predetermined, hard-coded lower bound.  The user is responsible for making 
! sure that sprayHFs() is only called for water points, i.e., sprayHFs() has no 
! internal check for land vs water.
!     The user selects the spray generation function to use in sprayHFs() by
! specifying one of the string keys in section 1.0 for the input variable
! whichSSGF.  If using the seastate-based model 'BCF23_Seastate', all wave
! parameters are required (eps, dcp, swh, mss).  If using a wind-based model
! ('F94_MOM80' or 'F94_BCF23'), swh is still required (for calculating the
! droplet settling timescale), but the remaining wave parameters (eps, dcp, mss)
! are not used (the user should pass dummy values, e.g., zeros or NaNs).
!     Subroutine sprayHFs() makes internal calculations for H_S0pr and H_L0pr, 
! as well as for ustar and the wind profile.  These internal values may not 
! exactly match those in the existing bulk code if, for instance, different 
! stability functions are used.  This is probably OK and should not affect 
! dHS1spr and dHL1spr very much.  Please contact Ben Barr if you are interested 
! in discussing ways to reconcile any discrepancies.
!     All INTENT(OUT) fields returned by sprayHFs() other than dHS1spr and
! dHL1spr are for debugging/interpretation and can be suppressed by the user if
! desired.
!
!===============================================================================

