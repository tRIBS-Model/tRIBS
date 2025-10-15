# Required CMake arguments
cmake_minimum_required(VERSION 3.20)
#set standard
set(CMAKE_CXX_STANDARD 17)
# Find MPI package with C++ support


set(exe "tRIBS_par")
set(cxx_flags "-O2")
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(${exe})
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_COMPILER mpicxx)

add_executable(${exe} src/main.cpp)
target_link_libraries(${exe} PUBLIC MPI::MPI_CXX)
# project files
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# used while testing and developing
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${cxx_flags}") # -O2 appears to be the fastest option

message(WARNING "Compiling parallel version of tRIBS")

target_compile_definitions(${exe}
        PRIVATE
        PARALLEL_TRIBS
        PRIVATE
        GRAPH_TRIBS
        PRIVATE
        LINUX_32
        )

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
        src/tHydro/tHydroMetStoch.cpp
        src/tHydro/tHydroMetStoch.h
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
        src/tStorm/tStorm.cpp
        src/tStorm/tStorm.h
    )
