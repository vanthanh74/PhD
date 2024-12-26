## Run manually
g++ parareal_heat_2D_BE.cpp -o parareal_heat_2D_BE && ./parareal_heat_2D_BE && python3 plot_tools_2D.py

## Run by CMake
mkdir build
cd build
cmake ..
make
make plot