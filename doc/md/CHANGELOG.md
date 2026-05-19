<!--- CHANGELOG.md --->
# Changelog 
All notable changes to this project are documented in this file.

## [6.0.0] - Unreleased
v6.0.0 is a major release focused on simplifying the user experience and streamlining the codebase. This version introduces improvements to snow physics, adds modern raster support via GDAL, and implements optimizations to improve simulation performance. **To view examples of the updated v6.0.0 input structure and the new snow parameter files (.spf), please visit our [Model Benchmark Repository](https://github.com/tRIBS-Model/tRIBS-benchmarks).**

### Added
* **GDAL Integration:** Added optional build configuration to link against GDAL, allowing tRIBS to read a wide variety of binary raster formats. Includes a new `WITH_GDAL` CMake flag. ([#103](https://github.com/tRIBS-Model/tRIBS/pull/103))
* **Static Land Use Grids:** Added `OPTLANDUSE = 2` to allow reading spatially variable but temporally constant land use parameters from non-timestamped rasters. ([#102](https://github.com/tRIBS-Model/tRIBS/pull/102))
* **Snow Parameter File:** Users can now provide a dedicated `.spf` file for snow physics constants, moving away from hardcoded defaults. ([#104](https://github.com/tRIBS-Model/tRIBS/pull/104))
* **New Snow Physics:** ([#104](https://github.com/tRIBS-Model/tRIBS/pull/104))
    * **Dynamic Density:** Implemented snow compaction and density evolution based on Jordan (1991) and Anderson (1976).
    * **Liquid Water Routing:** Added a physically-based routing scheme accounting for holding capacity and conductivity (Colbeck 1972).
    * **Phase Partitioning:** Added user-selectable thresholds for Wet-bulb or Air Temperature to determine precipitation phase.
* **Snow Outputs:** Added `Snow Depth` and `Snow Density` to standard pixel and dynamic output routines. ([#104](https://github.com/tRIBS-Model/tRIBS/pull/104))
* **Standardized Hydrologic Outputs:** All output files that report hydrologic variables have been standardized to CSV format with a single header line. ([#114](https://github.com/tRIBS-Model/tRIBS/pull/114))
* **Standardized Hydrologic Inputs:** All forcing and parameter input files were standardized to CSV format with a single header line. Centorid latitude, longitude, and UTC timezone offset are no longer specified in the station and gridded data files. ([#116](https://github.com/tRIBS-Model/tRIBS/pull/116))
* **Root Zone Depth:** Previously root zone depth used for determining soil mositure for calculating plant transpiration was hardcoded at 1m. Root zone depth is now a land use parameter that can be provided in the land use table or gridded. Users who do not want to use this feature can specify `9999.99` in the land use table to use the original default of 1m. ([#117](https://github.com/tRIBS-Model/tRIBS/pull/117))

### Fixed
* **ET Partitioning:** Fixed a scaling bug in `tEvapoTrans` where unscaled vegetation rates were subtracted from potential ET. ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
* **Snow Fluxes:** Removed the vegetation fraction scaling on sublimation and evaporation from the ground snowpack. ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
* **tIntercept Memory Leak:** Resolved a crash occurring at simulation termination due to improper deallocation of grid filenames when interception was disabled. ([#102](https://github.com/tRIBS-Model/tRIBS/pull/102))
* **Sub-hourly Precipitation:** Resolved a bug when using sub-hourly precipitation inputs with the snow module turned on caused by the snow module using the wrong timestep. ([#111](https://github.com/tRIBS-Model/tRIBS/pull/111))
* **Channel Transmission Losses:** Refactored and fixed multiple bugs related to the channel transmission losses. Previous versions capped transmission losses to the available volume of lateral fluxes into the stream node at each timestep. ([#112](https://github.com/tRIBS-Model/tRIBS/pull/112))
* **Surface Temperature Inputs:** Fixed multiple bugs related to point or gridded surface temperature input. This feature is valuable but was rarely used, now functioning as expected: eahc timestep surface temperature is caluclated using input surface temperature as the initial state. When surface temperature is not provided calculated surface temperature evolves freely from the energy balance (same as before). ([#116](https://github.com/tRIBS-Model/tRIBS/pull/116))

### Changed & Refactored
* **Input Simplification:** Streamlined the `.in` file by removing legacy or unused options including:
    * Stochastic storm generator. ([#97](https://github.com/tRIBS-Model/tRIBS/pull/97))
    * RIBS-output compatibility and alternative visualization options. ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
    * Conditional header writing (headers are now always written). ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
    * Simplified gridded rainfall options; the model now standardizes on mm/hr for both point and gridded data. ([#98](https://github.com/tRIBS-Model/tRIBS/pull/98))
    * Removed HydroMetConverter options for converting meteorologcial inputs from various sources internally. ([#96](https://github.com/tRIBS-Model/tRIBS/pull/98))
    * Simplified `OPTEVAPOTRANS` to support only the Penman-Monteith method for calculating ET. ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
    * Simplified `OPTMESHINPUT` to support only pre-generated 4-mesh files or points files. ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
    * Simplified `OPTINTERCEPT` to support only the Rutter method for calculating canopy interception. ([#109](https://github.com/tRIBS-Model/tRIBS/pull/109))
    * Removed multiple command line options that were either dead code or neaver used in practice. ([#113](https://github.com/tRIBS-Model/tRIBS/pull/113))
    * Removed legacy forecast module. ([#115](https://github.com/tRIBS-Model/tRIBS/pull/115))
* **C++ Performance Optimizations:** ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
    * **Memory Management:** Refactored geometry functions (`FindIntersectionCoords`, `PlaneFit`, and `setRVtx`) to pass `tArray<double>` by `const` reference, eliminating massive heap allocation overhead.
    * **Math:** Replaced `pow()` calls with `x*x` and `sqrt()` in core physics loops to reduce CPU cycles.
    * **Standardization:** Unified platform-specific headers and replaced `sprintf` with `snprintf` for modern compiler compatibility.
* **Error Warning Outputs:** While running the modle in parallel there were many error messages that would be duplicated across all processors or other issues resulting in massive log files. Many error messsages modified to only print 10 times before being silenced.
* **Solar Position Calculation Inputs:** Previously the input values for solar position calculations were required in multiple input files. They have been moved to the main input file under the keywords: `UTCOFFSET`, `CENTROIDLAT`, and `CENTROIDLONG`. ([#116](https://github.com/tRIBS-Model/tRIBS/pull/116))

### Removed
* Removed legacy code related to changes in main input file listed above.
* Removed unused platform-specific `#ifdef` blocks in headers. ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
* Removed legacy debug printing in `tTriangulator`. ([#107](https://github.com/tRIBS-Model/tRIBS/pull/107))
* Removed legacy command-line flags that were either never used or obsolete. ([#113](https://github.com/tRIBS-Model/tRIBS/pull/113))
* Removed optional humidity inputs. Relative humdity is the only accepted humidity forcing. ([#116](https://github.com/tRIBS-Model/tRIBS/pull/116))
* Removed optional `Kpan` parameter for `OPTEVPOTRANS = 2`. Any input ET forcing is now used directly as Potential ET. ([#116](https://github.com/tRIBS-Model/tRIBS/pull/116))

---

## [5.3.1] - 10/15/2025
### Added
* **Input Validation:** Added timestamp validation to rainfall and meteorological station input timeseries. Warning, older models will no longer run if there is missing data in the data files. ([#95](https://github.com/tRIBS-Model/tRIBS/pull/95))

### Fixed
* Fixed bug in Rutter interception scheme that could result in small amount of negative wet canopy evaporation. ([#94](https://github.com/tRIBS-Model/tRIBS/pull/94))
* Fixed bug when reading gridded landuse data that would result in the landuse table values being used instead under specific conditions.

---

## [5.3.0] - 08/16/2025
### Fixed
* **Canopy Water Balance:** Refactored `InterceptRutter` to calculate evaporation from the wet canopy internally. This corrects a small discrepancy in the canopy water balance previously handled in `tEvapoTrans`. ([#83](https://github.com/tRIBS-Model/tRIBS/pull/83))
* **Snow Interception:** Fixed vegetation fraction scaling; intercepted SWE now represents the actual state of the canopy rather than being incorrectly scaled. ([#83](https://github.com/tRIBS-Model/tRIBS/pull/83))
* **Raster Resampling:** Improved numerical stability of raster resampling in `tResample.cpp`. Fixed a bug with specific Voronoi polygon geometry that resulted in `NaN` values and mass balance errors. ([#86](https://github.com/tRIBS-Model/tRIBS/pull/86))
* **Mesh Geometry:** Fixed `ComputeVoronoiArea` function which was assigning incorrect area values to nodes connected to the outlet (node 0). ([#83](https://github.com/tRIBS-Model/tRIBS/pull/83))
* **Land Use Logic:** Resolved a bug where dynamic land use grids reverted to table values prematurely after the final interpolation interval. ([#85](https://github.com/tRIBS-Model/tRIBS/pull/85))
* **Dynamic Land Use:** Fixed incorrect model behavior when using dynamic land use grids with the interpolation option turned off (`luInterpOption = 0`). ([#85](https://github.com/tRIBS-Model/tRIBS/pull/85))
* **Timing:** Corrected `nodeHour` misalignment issues in `tEvapoTrans.cpp` and synchronized `julianDay()` and `SetSunVariables()` with the central `tRunTimer`. ([#84](https://github.com/tRIBS-Model/tRIBS/pull/84))

### Added
* **Soil Layer Control:** Added optional input parameters `SURFACESOILDEPTH` and `ROOTZONEDEPTH` (in mm) to allow user-defined layer depths (defaulting to 100mm and 1000mm for backward compatibility).
* **Gridded Parameters:** Added support for soil moisture stress thresholds for soil evaporation (`SE`) and plant transpiration (`ST`) in the gridded data file (`.gdf`). ([#85](https://github.com/tRIBS-Model/tRIBS/pull/85))
* **New Output Variables:** Added `Qunsat` to the Mean Response File (MRF) and `shortRadSlope` to pixel files. ([#83](https://github.com/tRIBS-Model/tRIBS/pull/83))
* **Error Handling:** Added a fatal error check for cases where `luInterpOption = 1` is selected but only a single raster is provided.

### Changed & Refactored
* **Output Standard:** Modified writing of soil water state variables in pixel, dynamic, and MRF files to convert from sloped state variables to vertical depths.
* **Solar Radiation:** Centralized slope, albedo, and vegetation corrections into `inShortWave()`, reducing redundancy and improving consistency across the energy balance module. ([#84](https://github.com/tRIBS-Model/tRIBS/pull/84))
* **ET Partitioning:** Updated `tEvapoTrans` to prioritize potential evaporation partitioning: first to wet canopy, then transpiration, and lastly soil evaporation.
* **Snow Physics Refactor:** ([#84](https://github.com/tRIBS-Model/tRIBS/pull/84))
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

