# tRIBS TODO

### Todo
- [ ] Refactor tResample::doIt, code is unclear and possible bottle neck.
- [ ] Update flow for reading in dynamic LU grids--this is a bottleneck in terms of speed. Potential speed improvements with caching resampling variables and others, tested but not implemented.
- [ ] Consider flexible approach for specifying outputs, I.E. could we make it so that you could pass in attribute related to a node and have them returned in the dynamic and integrated files. Could also be nice for input parameters, i.e. make list of hard coded parameters and expose as inputs.
- [ ] Resolve remaing problems in channel transmission losses routine for option 2 & 3, (transient & green-ampt methods) ([#112](https://github.com/tRIBS-Model/tRIBS/pull/112)).
- [ ] Update input timeseries and raster forcing files to use ISO 8601 timestamp format. No chnage to outputs for now, will still be elasped simulation hours. Example from `_MMddYYYYhh.asc` to `_20240615T1200.asc`
- [ ] The precip grid and gauge paths label intervals differently: gauge row `T` covers `[T, T+dt)` while grid file `T` covers `(T-dt, T]`, so identically stamped gauge and grid forcing are offset by one interval. Need to decide which, most radar data is end of period as implemented in the model already.
- [ ] `RAININTRVL` values that are not exact binary fractions (10, 5, 3, 2, 1 min) are still rejected by the `fmod(1.0, dtRain) != 0` check in `tRunTimer`.
- [ ] tRIBS `FillLakes` bug. When it looks for the lowest node around a pit, the check is if `(open boundary || z < lowestElev)`. Finding the outlet sets the running minimum to the outlet's elevation, so any later neighbour lower that that elev replaces it. The terminal pit therefore never drains to the outlet it's directly connected to.
- [ ] Potential bug, `vCell::convertToVoronoiFormat` seems to trust convexity of a polygon. Thus, sizes its buffers based on that which with bad triangles cause warning less problems. Look into how we can reject or repair self-intersecting cells with a warning.

### Finished
- [x] Finalize updated benchmarks
- [x] Remove invariant .pixel files--all relevant information can be written to the time integrated variable
- [x] Merged tSnowIntercept with tSnowPack
- [x] Fixed larger memory leaks when running in parallel
- [x] Update .in file and inputs, including moving command line flags and unnecessary options in .in files
- [x] Fix compiler warnings (sprintf to sprintnf, etc ...)
- [x] Update compiler flags and release version for optimized performance
- [x] Create parameter input file for snowpack.cpp
- [x] Check if GAUGEBASENAME: Rain Gauge data BASE name (*.mdf) is needed, this is superseded by the .sdf file.
- [x] Review/remove remaining user run flags
- [x] Remove legacy forecasting module
- [x] Simplify ET and interception schemes
- [x] Improve error output management
- [x] Clean up channel and hillslope routing parameters
- [x] Simplify humidity input. Only keep RH and remove others
- [x] Standardize input text files
- [x] Add root zone depth as land use parameter. If specified to -9999 then use original 1m.
- [x] Merge spatial outputs from parallel simulations within tRIBS
- [x] Generalize restart functionality
- [x] Integrate stomatal resistance changes from Becerra branch
- [x] Rerun tRIBS through a profiler for improving runtimes
- [x] Remove default parallelization option
- [x] Port MeshBuilder parallelization workflow into tRIBS
- [x] Package METIS from parallelization work with tRIBS



Return to [README](../../README.md)
