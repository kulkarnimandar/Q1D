# Q1D
MATLAB code for computing quasi-1D nozzle flow and sensitivities using Ccontinuum Sensitivity Analysis

run_q1d.m is the main script for running the primary and sensitivitiy anlysis. 

Set CSA_flag = 1 for CSA

Grid_convergence.m is the script which runs the grid convergence study with 7 different grids. Data is stored in .MAT files for each grid level.

About primary analysis:
- First order van Leer upwind scheme 
- T0 = T0_BC and P0 = P0_BC inflow boundary conditions, extrapolation outflow BC for isentropic subsonic-supersonic case works in this version

About CSA:
- Uses Jacobian of steady state primary analysis
- CSA BCs derived based on primary analysis BCs

Things to update in this version are written in issue #1
