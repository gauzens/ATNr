# Changes in version 1.1.1:
  - Correction of the temperature effects on parameters depending on both predator and prey body masses 

# Changes in version 1.1.0:
  - hill exponent q is now a vector. Each value gives the q of a specific consumer.
  - correcting bug in unscaled model: plant competition is now properly defined.
  - Vignettes have been improved
  - Extinctions thresholds are removed for nutrient pools in the unscaled_nut model
  - R.rsp package is now in the imports
  
# Changes in version 1.0.2:
  - The initialisation() function of each model that was supposed to be called before each integration is no internal to the lsoda wrapper function. Like that, there is no need to call this function each time a change is made (it is done whatever happens in lsoda_wrapper). Should be much less error prone.
  - The initialisation() for the scaled model was assuming that smallest basal species comes first. This is not the case anymore (body masses of plant do not have to be sorted anymore).
  - initialisation() is now initialisations()
  - Jacobian() is now jacobian()
  - Correction of Boolean comparisons in the C files (use logical instead of bitwise)
  - Vignettes 

# Changes in version 1.0.1:
  - Changing the extension of the latex vignette from .tex to .ltx as it was missing in the previous release. 
  - Correcting some memory management issues for the model codes.