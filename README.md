| Operating System | Build Status |
|------------------|--------------|
| Linux            | ![Linux Build](https://img.shields.io/github/actions/workflow/status/tribshms/tRIBS/compile_and_test_linux.yml) |
| macOS            | ![macOS Build](https://img.shields.io/github/actions/workflow/status/tribshms/tRIBS/compile_and_test_macos.yml)|
| Windows          | *Not Supported* |

![](https://img.shields.io/readthedocs/tribshms)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06747/status.svg)](https://doi.org/10.21105/joss.06747)

# TIN-based Real-time Integrated Basin Simulator: Version 6.0.0
This repository contains source code for the fully distributed hydrological model: TIN-based Real-time Integrated Basin Simulator (tRIBS). We provide extensive documentation of the model [here](https://tribshms.readthedocs.io/en/latest/). Our documentation also contains useful information about working with the tRIBS code base through [GitHub](https://tribshms.readthedocs.io/en/latest/man/Using%20GitHub.html) and how to [contribute](https://tribshms.readthedocs.io/en/latest/man/Contributing.html). Finally we also provide [templates](https://tribshms.readthedocs.io/en/latest/man/Templates.html).

Licensing information can be found in [LICENSE.txt](./LICENSE.txt).

Release notes are provided [here](https://tribshms.readthedocs.io/en/latest/man/Release%20Notes.html#).

Model reference list, number of citations and H-index can be found [here](https://csdms.colorado.edu/wiki/Model:TIN-based_Real-time_Integrated_Basin_Simulator_(tRIBS)).

## Installation 
We provide four options for accessing the tRIBS model:

1) Compiled executables for latest [macOS and Ubuntu](https://tribshms.readthedocs.io/en/latest/man/Executables.html#executables) systems
2) [Docker](https://tribshms.readthedocs.io/en/latest/man/Docker.html)
3) [CMake build system](https://tribshms.readthedocs.io/en/latest/man/Model_Execution.html#cmake)
4) [Codespaces](https://github.com/tRIBS-Model/tRIBS-Workshop-Sandbox)

## Getting Started
We have prepared multiple example applications for users to apply and learn from:

1) [Benchmarks](https://github.com/tRIBS-Model/tRIBS-benchmarks): Two fully setup models that users can run themselves or explore various model inputs. Requires users to supply their own installation of the model.
2) [Codespaces](https://github.com/tRIBS-Model/tRIBS-Workshop-Sandbox): An example application of generating an input file and running the model in a cloud computing environment, free of charge. Requires nothing from the user other than a Github account and a web-browser.
3) [pytRIBS Examples](https://github.com/tRIBS-Model/pytRIBS-examples): Example applications of full model generation in a Jupyter notebook using the pytRIBS python package. Requires the users to have pytRIBS installed in a python environment and Docker installed.

## Release/Version Notes
tRIBS uses semantic versioning. We record updates of major, minor, and patch versions [here](./doc/md/CHANGELOG.md).

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/tribshms/tRIBS)

> **Important Disclaimer**
>
> The following link points to an unofficial, third-party AI assistant that is not maintained or endorsed by the tRIBS development team.
>
> Answers provided by the AI may be **inaccurate, or incomplete.** For all scientific and technical applications, the official [tRIBS Documentation](https://tribshms.readthedocs.io/en/latest/) and the source code itself must be considered the **only source of truth**.
>
> Use this tool at your own risk.
