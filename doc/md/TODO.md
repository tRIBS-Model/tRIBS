# tRIBS TODO

### Todo
- [ ] Refactor tResample::doIt, code is unclear and possible bottle neck.
- [ ] Update flow for reading in dynamic LU grids--this is a bottleneck in terms of speed
- [ ] Consider flexible approach for specifying outputs, I.E. could we make it so that you could pass in attribute related to a node and have them returned in the dynamic and integrated files. Could also be nice for input parameters, i.e. make list of hard coded parameters and expose as inputs.
- [ ] Resolve remaing problems in channel transmission losses routine for option 2 & 3, (transient & green-ampt methods) ([#112](https://github.com/tRIBS-Model/tRIBS/pull/112)).
- [ ] Remove default parallelization option
- [ ] Port MeshBuilder parallelization workflow into tRIBS
- [ ] Package METIS from parallelization work with tRIBS
- [ ] Rerun tRIBS through a profiler for improving runtimes
- [ ] Integrate stomatal resistance changes from Becerra branch
- [ ] Generalize restart functionality
- [ ] Update input timeseries and raster forcing files to use ISO 8601 timestamp format. No chnage to outputs for now, will still be elasped simulation hours. Example from `_MMddYYYYhh.asc` to `_20240615T1200.asc`

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



Return to [README](../../README.md)
