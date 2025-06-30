# Lateral Misalignment Example

This example demonstrates how to calculate the mutual inducatance of two circular coils with lateral misalignment using the `ivsolver` binary. The example is implemented in MATLAB and uses the C++ backend for high-performance calculations.

The code is self-explanatory and includes comments to guide you through the process. The main steps are:
1. Define the wire parameters and create a wire definition file (`wires.def`). Define the parametric variables for the offset of the second coil, both z-axis and x-axis offsets.
2. Iterate over the defined offsets to calculate the mutual inductance for each configuration.
3. Run the `ivsolver` binary with the wire definition file and output file.
4. Plot the coupling coefficient as a function of the lateral misalignment.

If you want to visualise the two coils, uncomment the section that begins with `%% UNCOMMENT TO PLOT WINDING TRAJECTORIES %%` and run the code again. This will generate a 3D plot of the two coils with the specified offsets.
