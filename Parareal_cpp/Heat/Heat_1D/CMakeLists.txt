cmake_minimum_required(VERSION 3.10)
project(PararealHeatEquation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(parareal_heat_1d parareal_heat_1d.cpp)

# Add a custom target to execute the program
add_custom_target(run
    COMMAND parareal_heat_1d
    COMMENT "Running the Parareal Heat Equation simulation"
)

# Add a custom target to plot results using Python
add_custom_target(plot
	COMMAND cd .. && cp plot_tools.py ./build/ && cd build && python3 plot_tools.py
    #COMMAND python3 plot_tools.py
	COMMENT "Plotting results with Python3"
)

# Make "plot" depend on "run" to ensure simulation is executed before plotting
add_dependencies(plot run)
