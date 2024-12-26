## Run manually
g++ parareal_heat_1D.cpp -o parareal_heat_1D && ./parareal_heat_1D && python3 plot_tools.py
g++ parareal_heat_1D_BE.cpp -o parareal_heat_1D_BE && ./parareal_heat_1D_BE && python3 plot_tools.py

## Run by CMake
mkdir build
cd build
cmake ..
make
make plot