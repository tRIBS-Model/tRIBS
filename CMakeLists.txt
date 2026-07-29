# Required CMake arguments
cmake_minimum_required(VERSION 3.20)
#set standard
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY") #needed for testing cross compilation see: https://stackoverflow.com/questions/53633705/cmake-the-c-compiler-is-not-able-to-compile-a-simple-test-program

#set(parallel ON) # currently set through -Dparallel, but can be manually set by uncommenting
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Debug build-specific options
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -g") # -fsanitize=address -fsanitize-address-use-after-scope

# Release build-specific options, -O2 is the fastest
# -fno-math-errno: skips the (unused) errno bookkeeping around libm calls so
# the compiler can use hardware sqrt and streamline pow/exp; ~2% faster, does
# not change results. Supported by GCC, Clang, and Intel icpx alike.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -fno-math-errno")

option(WITH_GDAL "Build with GDAL support for GeoTIFF/NetCDF/etc." OFF)
option(WITH_PROFILER "Build with gperftools CPU profiler" OFF)

################################################################################
# VERSION AND BUILD PROVENANCE
#
# Single source of truth for the version: both project() calls below consume
# TRIBS_VERSION, and packaging/utilities/tribsVersion.h.in is configured into the
# build tree so the C++ code reads these same values rather than hardcoding them.
# Bump the version here and nowhere else.
#
# Must precede project(), which cannot be given a variable it does not yet know.
################################################################################
set(TRIBS_VERSION 6.0.0)
set(TRIBS_RELEASE "Summer 2026")

set(TRIBS_VERSION_STAMP_ARGS
        "-DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
        "-DTEMPLATE_FILE=${CMAKE_CURRENT_SOURCE_DIR}/packaging/utilities/tribsVersion.h.in"
        "-DOUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/generated/tribsVersion.h"
        "-DTRIBS_VERSION=${TRIBS_VERSION}"
        "-DTRIBS_RELEASE=${TRIBS_RELEASE}")
message(STATUS "tRIBS version ${TRIBS_VERSION} (${TRIBS_RELEASE})")

#set compiler if necessary, this has to be done before project()
#set (CMAKE_CXX_COMPILER <specify compiler here>)
#if(CMAKE_CXX_COMPILER_LOADED)
#    message(STATUS "Configured for C++, with compiler: ")
#    message(STATUS "${CMAKE_CXX_COMPILER_VERSION}}")
#endif()

#################### BELOW SPECIFIES COMPILATION FOR SERIAL OR PARALLEL VERSION OF TRIBS ###############################
if(parallel)
    set(exe "tRIBSpar")
    project(${exe}
            DESCRIPTION "Parallel Distribution of the TIN-based Real-time Integrated Basin Simulator Distributed Hydrological Model"
            VERSION ${TRIBS_VERSION}
    )
    add_executable(${exe} src/main.cpp)

    message(WARNING "Compiling parallel version of tRIBS")

    target_compile_definitions(${exe}
            PRIVATE
            PARALLEL_TRIBS
            PRIVATE
            GRAPH_TRIBS
    )

    find_package(MPI REQUIRED) #Note: changes to MPI path etc may require reset of CMakeCache.txt
    if (MPI_FOUND) #Note:this should find whichever implementation of MPI is installed (i.e. Open MPI, Intel, mpich)
        include_directories(${MPI_INCLUDE_PATH})
        target_link_libraries(${exe} PUBLIC MPI::MPI_CXX)
        set(CMAKE_CXX_COMPILER mpicxx)
    else (MPI_FOUND)
        message(SEND_ERROR "This application cannot compile without MPI")
    endif (MPI_FOUND)

    target_sources(
            ${exe} PRIVATE src/main.cpp
            src/Headers/Classes.h
            src/Headers/Definitions.h
            src/Headers/Inclusions.h
            src/Headers/TemplDefinitions.h
            src/Headers/globalFns.cpp
            src/Headers/globalFns.h
            src/Headers/globalIO.h
            src/Mathutil/geometry.h
            src/Mathutil/mathutil.cpp
            src/Mathutil/mathutil.h
            src/Mathutil/predicates.cpp
            src/Mathutil/predicates.h
            src/tArray/tArray.h
            src/tArray/tMatrix.cpp
            src/tArray/tMatrix.h
            src/tCNode/tCNode.cpp
            src/tCNode/tCNode.h
            src/tFlowNet/tFlowNet.cpp
            src/tFlowNet/tFlowNet.h
            src/tFlowNet/tFlowResults.cpp
            src/tFlowNet/tFlowResults.h
            src/tFlowNet/tKinemat.cpp
            src/tFlowNet/tKinemat.h
            src/tFlowNet/tResData.cpp
            src/tFlowNet/tResData.h
            src/tFlowNet/tReservoir.cpp
            src/tFlowNet/tReservoir.h
            src/tGraph/tGraph.cpp
            src/tGraph/tGraph.h
            src/tGraph/tGraphNode.cpp
            src/tGraph/tGraphNode.h
            src/tHydro/tEvapoTrans.cpp
            src/tHydro/tEvapoTrans.h
            src/tHydro/tHydroMet.cpp
            src/tHydro/tHydroMet.h
            src/tHydro/tHydroModel.cpp
            src/tHydro/tHydroModel.h
            src/tHydro/tIntercept.cpp
            src/tHydro/tIntercept.h
            src/tHydro/tSnowPack.cpp
            src/tHydro/tSnowPack.h
            src/tHydro/tWaterBalance.cpp
            src/tHydro/tWaterBalance.h
            src/tInOut/tInputFile.cpp
            src/tInOut/tInputFile.h
            src/tInOut/tOstream.cpp
            src/tInOut/tOstream.h
            src/tInOut/tOutput.cpp
            src/tInOut/tOutput.h
            src/tList/tList.cpp
            src/tList/tList.h
            src/tListInputData/tListInputData.cpp
            src/tListInputData/tListInputData.h
            src/tMesh/heapsort.h
            src/tMesh/tMesh.cpp
            src/tMesh/tMesh.h
            src/tMesh/tTriangulator.cpp
            src/tMesh/tTriangulator.h
            src/tMeshElements/meshElements.cpp
            src/tMeshElements/meshElements.h
            src/tMeshList/tMeshList.h
            src/tParallel/tParallel.cpp
            src/tParallel/tParallel.h
            src/tParallel/tTimer.cpp
            src/tParallel/tTimer.h
            src/tParallel/tTimings.cpp
            src/tParallel/tTimings.h
            src/tPtrList/tPtrList.cpp
            src/tPtrList/tPtrList.h
            src/tRasTin/tInvariant.cpp
            src/tRasTin/tInvariant.h
            src/tRasTin/tRainGauge.cpp
            src/tRasTin/tRainGauge.h
            src/tRasTin/tRainfall.cpp
            src/tRasTin/tRainfall.h
            src/tRasTin/tResample.cpp
            src/tRasTin/tResample.h
            src/tRasTin/tShelter.cpp
            src/tRasTin/tShelter.h
            src/tRasTin/tVariant.cpp
            src/tRasTin/tVariant.h
            src/tSimulator/tControl.cpp
            src/tSimulator/tControl.h
            src/tSimulator/tPreProcess.cpp
            src/tSimulator/tPreProcess.h
            src/tSimulator/tRestart.cpp
            src/tSimulator/tRestart.h
            src/tSimulator/tRunTimer.cpp
            src/tSimulator/tRunTimer.h
            src/tSimulator/tSimul.cpp
            src/tSimulator/tSimul.h
    )

else()
    set(exe "tRIBS")
    project(${exe}
            DESCRIPTION "Serial Distribution of the TIN-based Real-time Integrated Basin Simulator Distributed Hydrological Model"
            VERSION ${TRIBS_VERSION})
    add_executable(${exe} src/main.cpp)
    message(WARNING "Compiling serial version of tRIBS")

    target_sources(
            ${exe} PRIVATE src/main.cpp
            src/Headers/Classes.h
            src/Headers/Definitions.h
            src/Headers/Inclusions.h
            src/Headers/TemplDefinitions.h
            src/Headers/globalFns.cpp
            src/Headers/globalFns.h
            src/Headers/globalIO.h
            src/Mathutil/geometry.h
            src/Mathutil/mathutil.cpp
            src/Mathutil/mathutil.h
            src/Mathutil/predicates.cpp
            src/Mathutil/predicates.h
            src/tArray/tArray.h
            src/tArray/tMatrix.cpp
            src/tArray/tMatrix.h
            src/tCNode/tCNode.cpp
            src/tCNode/tCNode.h
            src/tFlowNet/tFlowNet.cpp
            src/tFlowNet/tFlowNet.h
            src/tFlowNet/tFlowResults.cpp
            src/tFlowNet/tFlowResults.h
            src/tFlowNet/tKinemat.cpp
            src/tFlowNet/tKinemat.h
            src/tFlowNet/tResData.cpp
            src/tFlowNet/tResData.h
            src/tFlowNet/tReservoir.cpp
            src/tFlowNet/tReservoir.h
            src/tGraph/tGraph.cpp
            src/tGraph/tGraph.h
            src/tGraph/tGraphNode.cpp
            src/tGraph/tGraphNode.h
            src/tHydro/tEvapoTrans.cpp
            src/tHydro/tEvapoTrans.h
            src/tHydro/tHydroMet.cpp
            src/tHydro/tHydroMet.h
            src/tHydro/tHydroModel.cpp
            src/tHydro/tHydroModel.h
            src/tHydro/tIntercept.cpp
            src/tHydro/tIntercept.h
            src/tHydro/tSnowPack.cpp
            src/tHydro/tSnowPack.h
            src/tHydro/tWaterBalance.cpp
            src/tHydro/tWaterBalance.h
            src/tInOut/tInputFile.cpp
            src/tInOut/tInputFile.h
            src/tInOut/tOstream.cpp
            src/tInOut/tOstream.h
            src/tInOut/tOutput.cpp
            src/tInOut/tOutput.h
            src/tList/tList.cpp
            src/tList/tList.h
            src/tListInputData/tListInputData.cpp
            src/tListInputData/tListInputData.h
            src/tMesh/heapsort.h
            src/tMesh/tMesh.cpp
            src/tMesh/tMesh.h
            src/tMesh/tTriangulator.cpp
            src/tMesh/tTriangulator.h
            src/tMeshElements/meshElements.cpp
            src/tMeshElements/meshElements.h
            src/tMeshList/tMeshList.h
            src/tPtrList/tPtrList.cpp
            src/tPtrList/tPtrList.h
            src/tRasTin/tInvariant.cpp
            src/tRasTin/tInvariant.h
            src/tRasTin/tRainGauge.cpp
            src/tRasTin/tRainGauge.h
            src/tRasTin/tRainfall.cpp
            src/tRasTin/tRainfall.h
            src/tRasTin/tResample.cpp
            src/tRasTin/tResample.h
            src/tRasTin/tShelter.cpp
            src/tRasTin/tShelter.h
            src/tRasTin/tVariant.cpp
            src/tRasTin/tVariant.h
            src/tSimulator/tControl.cpp
            src/tSimulator/tControl.h
            src/tSimulator/tPreProcess.cpp
            src/tSimulator/tPreProcess.h
            src/tSimulator/tRestart.cpp
            src/tSimulator/tRestart.h
            src/tSimulator/tRunTimer.cpp
            src/tSimulator/tRunTimer.h
            src/tSimulator/tSimul.cpp
            src/tSimulator/tSimul.h
    )

endif()

################################################################################
# GENERATED VERSION HEADER
# Applies to both Serial and Parallel builds (target: ${exe}).
#
# tribsVersion.h is written into the build tree, never the source tree, so it never
# appears as a modified file. This section sits after the serial/parallel branch
# because both target_include_directories() and add_dependencies() need ${exe}
# to already exist.
#
# The stamp runs twice, for two different reasons:
#   1. once here at configure time, so the header exists for IDEs and clangd
#      immediately after configuring;
#   2. on every build via the version_stamp target, so the git hash and dirty
#      flag track the tree as it is when the binary is actually produced.
#
# Only tControl.cpp includes tribsVersion.h; keep it out of Inclusions.h so that a
# version bump does not rebuild every translation unit.
################################################################################
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/generated)

execute_process(
        COMMAND ${CMAKE_COMMAND} ${TRIBS_VERSION_STAMP_ARGS}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/packaging/utilities/GitVersion.cmake)

add_custom_target(version_stamp ALL
        COMMAND ${CMAKE_COMMAND} ${TRIBS_VERSION_STAMP_ARGS}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/packaging/utilities/GitVersion.cmake
        COMMENT "Stamping git provenance into tribsVersion.h"
        VERBATIM)
add_dependencies(${exe} version_stamp)

target_include_directories(${exe} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)

################################################################################
# METIS + in-process reach partitioning (tPartition)
# Applies to both Serial and Parallel builds (target: ${exe}).
#
# Vendored GKlib + libmetis are compiled into a static library and linked into
# tRIBS, so stream-reach graph partitioning is an in-process
# METIS_PartGraphKway() call. 
#
# IDXTYPEWIDTH/REALTYPEWIDTH are commented out in metis.h and MUST be supplied.
################################################################################
set(METIS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/metis_builds/METIS)
set(GKLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/metis_builds/GKlib)

file(GLOB GKLIB_SOURCES ${GKLIB_DIR}/*.c)
file(GLOB METIS_SOURCES ${METIS_DIR}/libmetis/*.c)

add_library(metis STATIC ${GKLIB_SOURCES} ${METIS_SOURCES})
target_include_directories(metis
        PUBLIC  ${METIS_DIR}/include
        PRIVATE ${GKLIB_DIR} ${METIS_DIR}/libmetis)
target_compile_definitions(metis PUBLIC  IDXTYPEWIDTH=32 REALTYPEWIDTH=32)
target_compile_definitions(metis PRIVATE USE_GKREGEX)
if(UNIX AND NOT APPLE)
    target_compile_definitions(metis PRIVATE LINUX _FILE_OFFSET_BITS=64)
endif()
if(UNIX)
    target_link_libraries(metis PRIVATE m)
endif()
target_compile_options(metis PRIVATE -w)   # vendored C is not warning-clean
set_target_properties(metis PROPERTIES C_STANDARD 99)

target_sources(${exe} PRIVATE
        src/tPartition/tPartition.cpp
        src/tPartition/tPartition.h)
target_link_libraries(${exe} PRIVATE metis)

################################################################################
# GDAL CONFIGURATION
# This applies to both Serial and Parallel builds (target: ${exe})
################################################################################
if(WITH_GDAL)
    find_package(GDAL)
    
    if(GDAL_FOUND)
        message(STATUS "GDAL Found: ${GDAL_VERSION}")
        message(STATUS "Enabling GDAL support.")

        # 1. Define the preprocessor macro so the C++ code knows to compile the new function
        target_compile_definitions(${exe} PRIVATE READ_GDAL)

        # 2. Include header files
        target_include_directories(${exe} PRIVATE ${GDAL_INCLUDE_DIRS})

        # 3. Link the library
        # We use PRIVATE because we don't want to expose GDAL to consumers of tRIBS
        target_link_libraries(${exe} PRIVATE ${GDAL_LIBRARIES})

    else()
        message(WARNING "GDAL was requested (-DWITH_GDAL=ON) but not found on the system.")
        message(WARNING "Compiling WITHOUT GDAL support (ASCII Grids only).")
    endif()
endif()

################################################################################
# GPERFTOOLS PROFILER CONFIGURATION
################################################################################
if(WITH_PROFILER)
    find_library(PROFILER_LIBRARY profiler)
    if(PROFILER_LIBRARY)
        message(STATUS "gperftools profiler found: ${PROFILER_LIBRARY}")
        target_link_libraries(${exe} PRIVATE ${PROFILER_LIBRARY})
        target_compile_options(${exe} PRIVATE -fno-omit-frame-pointer)
    else()
        message(WARNING "gperftools profiler requested (-DWITH_PROFILER=ON) but not found.")
    endif()
endif()

# can be passed via the command line with -Dpackaging=ON
#set(packaging ON)

if (packaging)

    if (NOT DEFINED os)
        set(os "macOS/Silicon") # options: macOS/Silicon, macOS/Intel, or Linux
    endif ()

    install(
        TARGETS ${PROJECT_NAME}
    )

    set(CPACK_VERBOSE ON)
    set(CPACK_PACKAGE_NAME "TIN-based Real-time Integrated Basin Simulator")
    set(CPACK_PACKAGE_FILE_NAME "TIN-based Real-time Integrated Basin Simulator-${PROJECT_VERSION}")
    set(CPACK_PACKAGE_VENDOR "tRIBS developers")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE.txt")
    # Set the custom script for the STGZ generator
    #list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/packaging/utilities/CPack.STGZ_Header.sh.in")
    #set(CPACK_STGZ_HEADER_FILE "${CMAKE_SOURCE_DIR}/packaging/utilities/CPack.STGZ_Header.sh.in")

    set(CPACK_GENERATOR "STGZ")
    set(CPACK_PACKAGE_DIRECTORY "${CMAKE_SOURCE_DIR}/packaging/${os}")
    set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_SOURCE_DIR}/packaging/${os}")

    include(CPack)
    cpack_add_component(CPACK_COMPONENTS_ALL)

endif ()

