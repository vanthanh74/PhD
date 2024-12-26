## Run manually
g++ parareal_conv_react_diff_2D_BE.cpp -o parareal_conv_react_diff_2D_BE && ./parareal_conv_react_diff_2D_BE && python3 plot_tools_2D.py

## Run by CMake
mkdir build
cd build
cmake ..
make
make plot