<!--- CHANGELOG.md --->
# Changelog 
All notable changes to this project are documented in this file.

## [5.3.1] - 10/15/2025
### Added
* **Input Validation:** Added timestamp validation to rainfall and meteorological station input timeseries. Warning, older models will no longer run if there is missing data in the data files. PR #95

### Fixed
* Fixed bug in Rutter interception scheme that could result in small amount of negative wet canopy evaporation. (#94)
* Fixed bug when reading gridded landuse data that would result in the landuse table values being used instead under specific conditions.

---

## [5.3.0] - 08/16/2025
### Fixed (Scientific & Mass Balance)
* **Canopy Water Balance:** Refactored `InterceptRutter` to calculate evaporation from the wet canopy internally. This corrects a small discrepancy in the canopy water balance previously handled in `tEvapoTrans`. (#83)
* **Snow Interception:** Fixed vegetation fraction scaling; intercepted SWE now represents the actual state of the canopy rather than being incorrectly scaled. (#83)
* **Raster Resampling:** Improved numerical stability of raster resampling in `tResample.cpp`. Fixed a bug with specific Voronoi polygon geometry that resulted in `NaN` values and mass balance errors. (#86)
* **Mesh Geometry:** Fixed `ComputeVoronoiArea` function which was assigning incorrect area values to nodes connected to the outlet (node 0). (#83)
* **Land Use Logic:** Resolved a bug where dynamic land use grids reverted to table values prematurely after the final interpolation interval. (#85)
* **Dynamic Land Use:** Fixed incorrect model behavior when using dynamic land use grids with the interpolation option turned off (`luInterpOption = 0`). (#85)
* **Timing:** Corrected `nodeHour` misalignment issues in `tEvapoTrans.cpp` and synchronized `julianDay()` and `SetSunVariables()` with the central `tRunTimer`. (#84)

### Added
* **Soil Layer Control:** Added optional input parameters `SURFACESOILDEPTH` and `ROOTZONEDEPTH` (in mm) to allow user-defined layer depths (defaulting to 100mm and 1000mm for backward compatibility).
* **Gridded Parameters:** Added support for soil moisture stress thresholds for soil evaporation (`SE`) and plant transpiration (`ST`) in the gridded data file (`.gdf`).
* **New Output Variables:** Added `Qunsat` to the Mean Response File (MRF) and `shortRadSlope` to pixel files.
* **Error Handling:** Added a fatal error check for cases where `luInterpOption = 1` is selected but only a single raster is provided.

### Changed & Refactored
* **Output Standard:** Modified writing of soil water state variables in pixel, dynamic, and MRF files to convert from sloped state variables to vertical depths.
* **Solar Radiation:** Centralized slope, albedo, and vegetation corrections into `inShortWave()`, reducing redundancy and improving consistency across the energy balance module.
* **ET Partitioning:** Updated `tEvapoTrans` to prioritize potential evaporation partitioning: first to wet canopy, then transpiration, and lastly soil evaporation.
* **Snow Physics Refactor:** (#84)
    * Updated albedo decay function with a minimum albedo threshold.
    * Refactored latent and sensible heat flux for ground snowpack to prevent temperatures from dropping below zero when liquid water is present.
    * Incorporated Bulk Richardson Number stability correction for aerodynamic resistance.
* **Technical Cleanup:** 
    * Replaced non-standard Variable Length Arrays (VLAs) with `std::vector` for improved stability.
    * Removed hardcoded version numbers from source headers and updated copyright notices.
    * Cleaned up extraneous debug `cout` statements and fixed C17 compatibility standards.

---

## Version 5.2.1
### 6/21/2024
* Modified how tRIBS accounts for the condition when Nwt = Bedrock depth, through setting conditional cutoffs in the beta functions.
* Fixed some issues where meteorological variables were not correctly being initalized and set.
* Added black box testing feature for both Happy Jack and Big Spring benchmarks.
* Updated CMake with option to set a release or debug flag when building model.
###  5/12/2024
* Increased precision of time step in .qout files to two decimal points.
* Call to inShortWaveCan in tSnowPack::computeSub was commented out and needs to be fixed or removed
* Resolved issue related to Tso not properly being set.
* Setup proper initialization of skycover_flag
* Removed CNode::setrsrf(double value), was not being used.
* Removed unused variable BasAltitude, it was a red herring for debugging purposes, but otherwise useless
* Fixed bug, where vegetation fraction = 1, so that unloaded snow is not lost from the system. This was accomplished by conditionally checking if VF =1, then setting it to 0.99.
* Commented out 2012 modification of vegHeight for coeffH < 1 and resetting of coeffV = 0.1 in tSnowPack::resFactCalc. In specific cases (i.e. where coeffH < 1) may have caused some reduction in evaporation and sublimation from the snowpack.


## Version 5.2.0 — Summer/Fall 2023
The below information records some of the modifications leading to version 5.2.0. This documents the initial efforts to centralize the tRIBS code base with modernized  tools. Note this is not by any mean a complete record of modifications that occurred at this time. Additional information can be found in the Git logs.
### Added
- Added check to reservoir option in tKinemat.cpp
- Added catch tEvapoTrans to set evapSoil = 0 when bedrock depth = 0 and if coeffV cerr w/ exit(1) since this behavior is currently not represented in the model physics.
- Removed update of variable limit in tFlowResults object to avoid undefined behavior when restart function is used (tFlowNet.cpp L: 676).
- added tRIBSCodeReference.pdf to doc/ with doxygen
- updated build system with CMake functionality (CMakeLists.txt)
- Merged fixes from different versions of tRIBS code including from Josh Cederstrom, Ara Ko, Carlos Lizarraga, and Xiaoyang Tang
- added #include "tTimer.h" to tTimer.cpp
- markdown (md) subdirectory to display markdown files on github
- Catch for when groundwater == bedrock in tHydromodel
- Added docker file

### Fixed
- In tEvapoTrans::initialLUGridAssignment added num*Files>1 to conditional if statement to catch case where only one given landuse gird is available.
- Fixed issues with memory allocation related to reading files in tVariant and tEvapoTrans by replacing numeric values (i.e. 10) with the macro kMaxExt
- Merged/refactored tSnowIntercept.cpp with tSnowPack;
    - There was an issue of creating a separate instance of tEvapoTrans in tSnowIntercept as it this instance was never initialized and a probable source undefined behavior, including calls to read in meteorological data from station files.
    - To fix this issue I simplified both tSnowPack and tSnowIntercept by removing unused functions and variables, replaced code snippets that were pulled from tEvapoTrans functions (e.g., interpolatLUGrids, integratedLUVars, etc), created new function for self-contained sections of code (i.e checkShelter, updateRipeSnowPack, etc..) and moved the content of tSnowIntercept into tSnowPack.cpp. Note tRestart.cpp, tSimul.cpp, and CMakeLists.txt also needed to be updated to account for references to tSnowIntercept.
    - Other minor fixes include replacing 3.1416 with macro PI in sublimation functions for callSnowIntercept and commenting out define albedo in callSnowIntercept since that is updated in tSnowPack. Also added in Xiaoyang's fix for routing liquid for snow pack in the case of no precipitation heat flux.
- Updated tOutput.cpp.
  -In CreatAndOpenPixel removed numbers from the front of variables, replaced commas, and forward slashes with underscores. This was done to facilitate easier reading of .pixel files into python.
    - Increased precision in WritePixelInfo to standard of 7 for all variables. Precision varied from 1 to 7 prior to this change. At somepoint someone might want to evaluate if different levels of precision are warranted--but I updated to improve post model run calcuations on water balance estimates. Increased precision showed a minor but still noticeable change in values.
    - In CreateAndOpenOutLet and SetrInteriorOulet I put an >= to if statement checking outlets as nodeId =0  is valid, especially in single element runs.
    - - Fixed Compiler errors for Linux HPC
- Fixed multiple issues in tSnow classes
- Compiler errors related to assert statements with null pointers
- Compiler error for tPtrList.h (L 873) newlist.insertAtBack to newlist->insertAtBack
- Issue in tResample::convertToVoronoiFormat where L 1615-1617 InOrOut variable was not allocating enough memory
- If optres == 0, would define tReservoir/tResData as dangling pointer, fixed by making tReservoir public and setting as null if optres == 0.
- Commented out #include t*(parallel code).cpp in parallel header files because it led to redefinition
- Fixed memory errors related to parallel operation 

### Removed
- tSnowIntercept.cpp see above for details
- doc/doc/ (created by doxygen) will update in future with pdf of readthedocs verison
- removed register calls (no longer supported at c++ 17 or earlier)
- removed old make file

### Changed
- layout of folder structure all tRIBS source code is now in the [src](./../src) folder.

## Return to [README](../../README.md)

