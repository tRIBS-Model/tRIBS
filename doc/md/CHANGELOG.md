<!--- CHANGELOG.md --->
# Changelog 
All notable changes to this project are documented in this file.

## Version 6.0.0
### 10/16/2025
* Removed multiple legacy or unused options in the input file (stochastic storm generator, hydrometeorological converter, RIBS output compatibility, alternative visualization options).
* The only available option for calculating ET from meteorological data is Penman-Monteith.
* There are no longer multiple methods to provided gridded rainfall data. The user can provide either gridded or point station data, both in units of mm/hr.
* The formats for providing a tin mesh to the model with `OPTMESHINPUT` have been simplified. Now there are only two options, a pre-generated mesh from the 4 mesh files or a point file.
* Removed conditional header writing in model outputs. Headers are now always written.
* Added new `OPTLANDUSE` option 2 for reading static landuse parameter rasters without needing dates in the filenames.
* Added capability to read in raster data using GDAL. To activate the GDAL features the user must have GDAL installed already. When compiling tRIBS there is a new flag, `WITH_GDAL`, for turning on GDAL support. This feature allows tRIBS to read any raster type that GDAL is capable above but note that tRIBS does not support multi-band rasters. If a user does not enable `WITH_GDAL` while compiling the model will function the exact same as previous versions.
## Version 5.3.1
### 10/15/2025
* Fixed bug in Rutter interception scheme that could result in small amount of negative wet canopy evaporation.
* Added timestamp validation to meteorological station input timeseries.
* Fixed bug when reading gridded landuse data that would result in the landuse table values being used instead under specific conditions.
## Version 5.3.0
### 8/16/2025
* Removed extraneous cout statement that prints out `OptRES` during initialization.
* Removed hardcoded version number from the top of all source code files to conform to modern standards and simplify future model updates.
* Replaced non-standard Variable Length Arrays (VLAs) with std::vector in tResample.cpp to resolve compiler warnings and improve code portability.
* Improved numerical stability of raster resampling in tResample.cpp. This change fixes a bug with specific voronoi polygon geometry that would result in a NaN values for gridded parameters.
### 8/11/2025
* Added new optional input for the gridded land use parameters. These parameters are the soil moisture stress thresholds for soil evaporation and plant transpiration, denoted by `SE` and `ST` in the gridded data file (.gdf), respectively.
* Resolved a bug that could cause incorrect model behavior when using dynamic land use grids with the interpolation option turned off (`luInterpOption = 0`).
Resolved a bug that would cause gridded land use parameters to revert back to the table values after the final interpolation interval. Now the model will hold the final grid's value until the end of the simulation.
* Added a fatal error check for when the interpolation option turned off (`luInterpOption = 1`) but the user only supplied a single raster. Previously, the model would crash without warning.
### 8/3/2025 
* Updated copyright notices to the current year.
* Added new optional input file parameters for specifying the depth that defines surface soil and root zone soil moisture. These can be specified with `SURFACESOILDEPTH` and `ROOTZONEDEPTH` in mm. If not specified, these values default to their original hardcoded values of 100mm and 1,000mm, respectively, ensuring backward compatibility.
* Updated how surface soil moisture is calculated in tEvapoTrans. Now the calculated potential evaporation is partitioned first to the wet canopy, then transpiration, and lastly soil evaporation.
### 7/28/2025
* Refactored snow albedo decay function in tSnowPack to have descriptive variable names and set a minimum albedo threshold for which albedo cannot go below.
* Refactor the calculation of latent and sensible heat flux for the ground snowpack in tSnowPack to prevent the snowpack temperature to be less than zero when there is liquid water present.
* Incorporate a stability correction factor based on the Bulk Richardson Number for the calculation of aerodynamic resistance in tSnowPack.
### 7/25/2025
* Refactor the function InterceptRutter to now calculate evaporation from the wet canopy rather than in tEvapoTrans. This change corrects a small error in the canopy water balance.
* Modify the function ComputeETComponents to handle the changes in tIntercept.cpp.
* Fix vegetation fraction scaling in snow interception scheme. Intercepted SWE now represents the actual state of the canopy rather than scaled by the vegetation fraction.
* Fix ComputeVoronoiArea function in meshElements.cpp that was setting the incorrect voronoi area value for the node connected to the outlet node in the mesh.
* Modify writing of soil water state variables to pixel, dynamic, and mrf outputs files to convert from sloped state variables to vertical depths.
* Add new variable, Qunsat, to the mean response file.
### 6/5/2025
* Refactored solar radiation handling for clarity and consistency:
* Centralized all slope, albedo, and vegetation corrections into inShortWave(), reducing redundancy in energyBalance().
* Reorganized computation of direct and diffuse radiation components using observed and computed values where applicable.
* Introduced new pixel file variable shortRadSlope to store terrain-corrected incoming shortwave radiation.
* Added setShortRadSlope() and getShortRadSlope() to the tCNode class.
* Corrected nodeHour misalignment issues in tEvapoTrans.cpp.
* Updated julianDay() and SetSunVariables() to use consistent time source from tRunTimer.
* Fixed declaration scope of RadGlobClr to clarify its limited usage within conditional blocks.
* Fixed multiple small bugs that were not meeting C17 compatibility standards.


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

