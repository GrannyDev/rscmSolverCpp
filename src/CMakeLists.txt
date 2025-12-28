cmake_minimum_required(VERSION 3.30)
project(rscmSolver)

set(CMAKE_CXX_STANDARD 26)

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

target_link_libraries(rscmSolver ortools::ortools)
