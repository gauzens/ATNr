Changes in version 1.0.2:
  The initialisation() function of each model that was supposed to be called before each integration is no internal to the lsoda wrapper function. Like that, there is no need to call this function each time a change is made (it is done whatever happens in lsoda_wrapper). Should be much less error prone.  

Changes in version 1.0.1:
  Changing the extension of the latex vignette from .tex to .ltx as it was missing in the previous release. 
  Correcting some memory management issues for the model codes.