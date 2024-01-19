# Energy Hub System Optimization

This project consists of MATLAB scripts designed for optimizing an energy hub system. The scripts handle various aspects of the energy system, including data reading, network calculation, and optimization (both deterministic and robust).
## Prerequisites
- MATLAB (version R2020a or later recommended)
- Optimization Toolbox
- Matpowers
- Symbolic Math Toolbox
- Optimisation Toolbox
- Global Optimization Toolbox

## Installation
Clone the repository or download the scripts directly into your MATLAB working directory.

## Usage
1. Compile with your yearly data the "Dati Hub.xlsx" file
2. Run `main_EH.m` to start optimization process.

## Scripts Overview

### `Dati Hub.xlsx`
- **Description**: contains various types of energy data (load, thermal, gas, cooling, irradiation) from an Excel file, scales the irradiation data, and reshapes it into matrices.

### `calc_network.m`
- **Description**: Calculates the electrical parameters like voltages and currents in an electrical network.
- **Key Functions**: `exp`, `abs`, `get_losses`, `makeYbus`.

### `idx_variables.m`
- **Description**: Defines and allocates indices for variables in an energy system optimization model.
- **Key Functions**: None (variable indexing).

### `deterministic_hub.m`
- **Description**: Sets up a deterministic optimization problem for an energy hub system. Defines parameters, constraints, and objective function.
- **Key Functions**: `linprog`, `zeros`, `Inf`.

### `robust_hub.m`
- **Description**: Similar to `deterministic_hub` but incorporates robust optimization elements to handle uncertainties.
- **Key Functions**: `gamma_def`, `zeros`, `Inf`.

## License
This software is licensed under MIT license. For more information read license.txt

## Contact
Project Link: https://github.com/RiccTrev/RobustOptimization
