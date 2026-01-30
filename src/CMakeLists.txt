cmake_minimum_required(VERSION 3.30)
project(rscmSolver)

set(CMAKE_CXX_STANDARD 26)

# The static OR-Tools package expects this imported target; provide it if missing.
if(NOT TARGET BZip2::bz2_static)
    add_library(BZip2::bz2_static STATIC IMPORTED)
    set_target_properties(BZip2::bz2_static PROPERTIES
        IMPORTED_LOCATION "/home/smith/Compiled/or-tools/build_static/install/lib/libbz2_static.a"
        INTERFACE_INCLUDE_DIRECTORIES "/home/smith/Compiled/or-tools/build_static/install/include"
    )
endif()

find_package(ortools CONFIG REQUIRED)

add_executable(rscmSolver main.cpp
        solver/Solver.cpp
        solver/Solver.h
        shared/VariableDefs.h
        shared/Variables.cpp
        shared/Variables.h
        shared/Adder.cpp
        shared/Adder.h
        shared/Layer.cpp
        shared/Layer.h
        cp/CPModel.cpp
        cp/CPModel.h
        solver/RSCM.cpp
        solver/RSCM.h
        shared/DAG.cpp
        shared/DAG.h
        writers/Verilog.cpp
        writers/Verilog.h
        writers/JSONDumper.cpp
        writers/JSONDumper.h
        writers/SnapshotIO.cpp
        writers/SnapshotIO.h
        solver/CostModel.h
        solver/CostComputer.cpp
        solver/CostComputer.h
        solver/ICostComputer.h
)

target_link_libraries(rscmSolver PRIVATE ortools::ortools)
