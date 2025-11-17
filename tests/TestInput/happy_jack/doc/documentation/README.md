# TIN-Based Real-Time Integrated Basin Simulator point-scale benchmark
This repository hosts the setup for executing a point-scale model of the [Happy Jack SNOTEL Site](https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=969) in Northern Arizona, USA, using the TIN-Based Real-Time Integrated Basin Simulator ([tRIBS](https://tribshms.readthedocs.io/en/latest/)). tRIBS v5.2 (or later) uses CMake as a build system, instructions on downloading and building tRIBS can be found [here](https://tribshms.readthedocs.io/en/latest/man/Model_Execution.html#compilation-instructions). Or alternatively one may use the tRIBS docker image, more information on this can be found [here](https://tribshms.readthedocs.io/en/latest/man/Docker.html#docker).

## Model Execution
From the command line tRIBS can be executed as follows, assuming the executable is stored in the sub-directory `bin`:

```
bin/tRIBS <path/to/infile>
```

In this case the happy jack input files can be found in ```src/in_files/happy_jack.in``` and can be further modified to explore tRIBS functionality. Note: the current input file is setup with relative paths and results will be saved to ```results/test```. Also, any modification in ```data/model``` may require an update of the .in file.

## Benchmark 
The benchmark case results are stored as a zip file under ```results/reference```. These results can be visualized in comparison to SNOTEL[^1]<sup>,</sup>[^2] data for the Happy Jack site as demonstrated in the jupyter notebook in ```src/tRIBS_snotel_comparison.ipynb```.

## Directory Structure:
### data
Contains all the necessary data to run tRIBS at Happy Jack and includes Snotel and SWANN data as a calibration/validation set.
### doc 
Contains notebooks for running and analyzing this specific benchmark case, along side additional documentation.
### src
Is designed to contain source code for for the tRIBS executable, which can be obtained [here](https://github.com/tribshms/tRIBS).
### bin
Directory for building and storing tRIBS executable, with instructions [here](https://tribshms.readthedocs.io/en/latest/man/Model_Execution.html#compilation-instructions).
### results
Directory for results, with a _reference_ and _test_ sub-directories. The former contains reference outputs from the Happy Jack point tRIBS simulation, while the later is empty and intended to store additional model simulations. Note the main results also contains a zip of the _reference_ sub-directory.



## References

[^1]: Sun N, H Yan, M Wigmosta, R Skaggs, R Leung, and Z Hou. 2019. “Regional snow parameters estimation for large-domain hydrological applications in the western United States.” Journal of Geophysical Research: Atmospheres. doi: 10.1029/2018JD030140

[^2]: Yan H, N Sun, M Wigmosta, R Skaggs, Z Hou, and R Leung. 2018. “Next-generation intensity-duration-frequency curves for hydrologic design in snow-dominated environments.” Water Resources Research, 54(2), 1093–1108.
BCQC Data Format

