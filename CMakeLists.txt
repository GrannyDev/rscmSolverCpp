cmake_minimum_required(VERSION 3.30)
project(rmcmSolver)

set(CMAKE_CXX_STANDARD 26)

find_package(ortools CONFIG REQUIRED)

add_executable(rmcmSolver main.cpp
        Solver.cpp
        Solver.h
        VariableDefs.h
        Variables.cpp
        Variables.h
        Adder.cpp
        Adder.h
        Layer.cpp
        Layer.h
        CPModel.cpp
        CPModel.h
        RSCM.cpp
        RSCM.h
        DAG.cpp
        DAG.h
        Verilog.cpp
        Verilog.h
        JSONDumper.cpp
        JSONDumper.h
)

target_link_libraries(rmcmSolver ortools::ortools)
