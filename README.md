# IVSolver
IPTVisual solver implemented in C++. A fast winding modeler and inductance calculator for arbitrary air-core windings. Used with parametric optimisation programmes in MATLAB, Python, etc..

## Features

- Efficient C++ implementation for high-performance calculations
- Supports modeling of arbitrary air-core winding geometries
- Calculates self and mutual inductance for complex coil arrangements
- Command-line interface for batch processing and integration with scripts
- Designed for use with optimization workflows in MATLAB, Python, and other environments
- Modular code structure for easy extension and customization
- Includes example usage and test cases in the `/examples` directory

## Installation
1. Clone the repository:
   ```bash
   git clone git@github.com:guszhang/ivsolver
   cd ivsolver
   ```

2. Build the project:
   on Linux:
   ```bash
   make
   ```

   on Windows:
   Use the provided Visual Studio solution file to build the project.

3. Test the installation:
   Run the example scripts in the `/examples` directory to verify that the solver works correctly.

## Usage
The solver takes text input files of definitions of wires. It has two modes of operation:

### **Impedance Calculation**: Computes the self and mutual inductance of the defined wires.

usage example:
   ```bash
   ivsolver wire_def output_file
   ```

   where `wire_def` contains the definitions of the wires and `output_file` will contain the calculated inductance values.

   The input file format is as follows:
   - First line: number of wires `n`. An positive interger number. The next `n` sections define each wire.
   - First line of each wire section: `cs_radius` (radius of the cross-section in mm, a positive float number), `m` number of commands (a positive integer number). The next `m` sections define the commands.
   - First line of each command section: `command_directive`. One of the following chars:
      - `l`: linear segment, followed by `k` number of points (a positive integer number), and then `k` lines of coordinates in the format `x y z` (three float numbers). It connects all key points with lines. It will automatically connect the previous last point (if exist) to the first point in this list.
      - `s`: spiral segment, followed by `k` number of points (a positive integer number), and then `k` lines of coordinates in the format `x y z` (three float numbers). It connects all key points with linearly interpolated spiral arcs along the z-axis. It will automatically connect the previous last point (if exist) to the first point in this list. 
      - `g`: saggy segment with $a$ parameter, followed by `k` number of points (a positive integer number) and `a` which is the sagging parameter as in $z=a \cosh(\frac{r-r_0}{a}+z_0)$. The next `k` lines are coordinates in the format `x y z` (three float numbers). It connects all key points with saggy arcs along the z-axis. It will automatically connect the previous last point (if exist) to the first point in this list. Please see the [saggy segment documentation]().
      - `c`: saggy segment with known length, followed by `k` number of points (a positive integer number). The next `k` lines are coordinates and length of wire in the format `x y z l` (three float numbers and one float number for length). It connects all key points with saggy arcs along the z-axis, where the sagging is determined by the length of the wire segment. It will automatically connect the previous last point (if exist) to the first point in this list. Please see the [saggy segment documentation](./docs/catenary_lines.md).
      - `x`: followed by `angle` (a float number in degrees), which rotates all previous wires by the specified angle around the x-axis. 
      - `y`: followed by `angle` (a float number in degrees), which rotates all previous wires by the specified angle around the y-axis.
      - `z`: followed by `angle` (a float number in degrees), which rotates all previous wires by the specified angle around the z-axis.
      - `o`: followed by `x y z` (three float numbers), which translates all previous wires by the specified vector in 3D space.
      - `n`: followed by `x y z` (three float numbers), which moves all previous wires with the original z-positive-direction to the specified point in 3D space, effectively setting the new normal direction for the wires if the wires were symmetric around the z-axis.
      - `b`: followed by `r` (a positive float number), which sets the bending radius of the wire to `r` m around the x-axis. This is used to bend a planar winding fitting a cylindrical surface. You may need to apply transformations several times if it is not bending around the x-axis.
      - `i`: followed by `real imaginary` (two float numbers), which sets the wire's current to the specified complex value. This is used for magnetic field calculations.
  
   All coordinates are in meteres, and angles are in degrees. The wires are defined in a 3D Cartesian coordinate system.
   
   After running the solver, `n+1` output files are generated with the designated path/name. The `.lmat` file contains the L matrix of size `n` x `n`, where the diagonal elements represent the self-inductance of each wire, and the off-diagonal elements represent the mutual inductance between pairs of wires. The unit is $\mathrm{H}$ (Henries). The last line of the `.lmat` file is the total length of each wire in meters. The `_x.trj` files where `x` is the wire index (0 to n-1) contain the trajectory points of each wire in the format `x y z` (three float numbers).

   The computation of the inductance matrix is based on the Neumann's formula, which is a well-known method for calculating the self and mutual inductance of wire configurations. The wire generation programme segments the wires into small segments of length equal to the wire diameter, and the inductance is calculated using the following formula given by Maxwell:
   $$
   L_{ij} = \frac{\mu_0}{4\pi} \int_{l_i} \int_{l_j} \frac{\mathrm{d}\mathbf{l}_i \cdot \mathrm{d}\mathbf{l}_j}{|\mathbf{r}_i - \mathbf{r}_j|}
   $$
   For self-inductance calculation where $\mathrm{d}\mathbf{l}_i$ and $\mathrm{d}\mathbf{l}_j$ refer to the same segment, a closed-form kernel equation is used to avoid divide-by-zero.
   $$
   \mathrm{d}L = \frac{\mu_0}{2\pi} \left( l\ln{\left(\frac{l + \sqrt{l^2 + a^2}}{a}\right)}  - \sqrt{l^2 + a^2} + \frac{l}{4} +a  \right)
   $$
   where $l$ is the length of the segment, $a$ is the radius of the wire, and $\mu_0$ is the permeability of free space. The solver uses numerical integration to compute the contributions from each segment of the wires.

   The accuracy of this method was verified against the analytical solutions for simple geometries and ANSYS simulations, reported in the appendix of the [IPTVisual paper](https://doi.org/10.3390/wevj13040063).

### **Field Calculation**: Computes the magnetic field at specified observation points.
   usage example:
   ```bash
   ivsolver wire_def output_file -b observation_def
   ```

   where `wire_def` contains the definitions of the wires, `output_file` will contain the calculated magnetic field values, and `observation_def` contains the definition of coordinates of the observation points.

   The input file format is the same to the impedance calculation. Notice that a current must be defined for each wire using the `i` command. The observation points are defined in a separate `observation_def` file, which has the following format:
   - First line: `x_min x_max x_n` (two float numbers and one positive integer number), which defines the range of x-coordinates and the number of points in that range, both inclusive.
   - Second line: `y_min y_max y_n` (two float numbers and one positive integer number), which defines the range of y-coordinates and the number of points in that range, both inclusive.
   - Third line: `z_min z_max z_n` (two float numbers and one positive integer number), which defines the range of z-coordinates and the number of points in that range, both inclusive.
  
   After running the solver, an output file is generated with the designated path/name. The file contains the magnetic field values at each observation point in the format `x y z Bx By Bz`, where `Bx`, `By`, and `Bz` are the components of the magnetic field vector in complex format at that point.

   **Complex Current and Field**: The solver supports complex currents and fields for single frequency ac analysis, so that currents of different phases can be defined. The calculation is based on the Biot-Savart law, which states that the magnetic field at a point in space is proportional to the current flowing through a wire and inversely proportional to the distance from the wire to the point.
   $$
   \mathbf{B}(\mathbf{r}) = \frac{\mu_0}{4\pi} \int_{l} \frac{I(\mathbf{r'}) \mathrm{d}\mathbf{l'}}{|\mathbf{r} - \mathbf{r'}|^3} \times (\mathbf{r} - \mathbf{r'})
   $$
   where $\mathbf{B}(\mathbf{r})$ is the magnetic field at point $\mathbf{r}$, $I(\mathbf{r'})$ is the current at point $\mathbf{r'}$, and $\mathrm{d}\mathbf{l'}$ is the differential length element along the wire. The solver uses numerical integration to compute the contributions from each segment of the wires.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss changes or improvements. Please contact me if you have any questions or suggestions.

**Author:** Gus Cheng Zhang, The University of Manchester, UK.

## License
This project is licensed under the GNU Affero General Public License v3.0. See the [LICENSE](LICENSE) file for details.


