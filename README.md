# Stochastic Rotation Dynamics (SRD) Simulation for Heterogeneous Catalysis (Code2)

## Project Overview

This project implements a particle-based simulation framework using Stochastic Rotation Dynamics (SRD) to model complex phenomena in heterogeneous catalytic reactors. It focuses on understanding multiscale processes, including convection, diffusion, and chemical reactions, particularly at catalytic surfaces. The simulation code is written in C/C++ and is complemented by Mathematica notebooks for data analysis and visualization.

This version of the code (`Code2`) is closely related to the research published in:
1.  A. Sengar et al., "Towards a particle based approach for multiscale modeling of heterogeneous catalytic reactors," *Chemical Engineering Science*, vol. 198, pp. 184-197, 2019.
2.  A. Sengar et al., "Particle-based modeling of heterogeneous chemical kinetics including mass transfer," *Physical Review E*, vol. 96, no. 2, p. 022115, 2017.

## Purpose and Functionality

The primary goal of this project is to bridge the gap between microscopic molecular events and macroscopic fluid dynamics in reactive systems. Key functionalities include:

*   **Stochastic Rotation Dynamics (SRD):** A mesoscale simulation technique that models fluid as coarse-grained particles undergoing streaming and collision steps.
*   **Multicomponent Diffusion:** Simulation of mixtures with different particle masses and concentrations, enabling the study of mutual diffusion.
*   **Heterogeneous Catalysis:** Detailed modeling of adsorption, desorption, and surface reaction kinetics (e.g., Langmuir-Hinshelwood) on catalytic surfaces.
*   **Boundary Conditions:** Support for both periodic boundary conditions and solid walls with bounce-back rules, including ghost particle methods for accurate hydrodynamics near walls. **Notably, this version (`Code2`) introduces an "open" boundary condition in the z-direction, allowing particles to leave the system, in addition to periodic boundaries.**
*   **Thermostat and Barostat:** Implementation of a Galilean invariant thermostat to maintain constant temperature and ensure physical realism. (While a barostat is mentioned in the context, the provided code primarily focuses on the thermostat).
*   **Flow Profiles:** Ability to impose and analyze flow profiles within the reactor.
*   **Data Analysis:** Functions to calculate key physical quantities such as:
    *   Mean Squared Displacement (MSD) for diffusion coefficients.
    *   Velocity Autocorrelation Functions (ACF) for self and mutual diffusion.
    *   Terminal velocity method for diffusion.
    *   Velocity and concentration profiles across the reactor.
*   **Axial Dispersion Method:** A specialized method to study dispersion in packed beds or similar geometries.

## Key Features and Contributions

This project represents a significant effort in developing a robust simulation tool for chemical engineering research. My contributions encompassed the entire development cycle, including:

*   **Core SRD Algorithm Implementation:** Developed the fundamental SRD simulation engine in C/C++.
*   **Custom Thermostat and Barostat:** Implemented a Galilean invariant thermostat to accurately control temperature during simulations. (The code primarily shows thermostat implementation, but the underlying principles can extend to barostats).
*   **Complex Interaction Potentials:** Designed and implemented particle-surface interaction potentials to model adsorption, desorption, and reaction events on catalytic sites.
*   **Advanced Boundary Conditions:** Developed and integrated sophisticated boundary conditions, including ghost particle methods for solid walls, to ensure accurate fluid-surface interactions. **The introduction of an open boundary condition in the z-direction in this version allows for simulating systems where particles can exit the simulation domain.**
*   **Multicomponent System Handling:** Extended the SRD framework to handle mixtures of particles with different masses and properties, crucial for studying multicomponent diffusion.
*   **Reaction Kinetics Integration:** Incorporated various reaction kinetics models, including Langmuir-Hinshelwood, directly into the particle-based framework.
*   **Extensive Data Analysis Tools:** Built in-house functions for calculating and analyzing various correlation functions and profiles, essential for extracting physical insights from simulation data.
*   **Large-Scale Data Generation:** The simulations are capable of generating **tens of gigabytes (10s of GBs)** of data, including particle trajectories, velocity fields, and concentration profiles, which are then processed for analysis.
*   **Mathematica Notebooks:** Developed Mathematica notebooks for post-processing, visualization, and further analytical treatment of the simulation results.

## Project Size

The codebase consists of several thousand lines of C/C++ code, meticulously developed and optimized for performance. Additionally, a significant number of lines of Mathematica code are used for data analysis, visualization, and generating publication-quality figures. The project involves handling and processing large datasets, typically in the range of **tens of gigabytes (10s of GBs)** per simulation run, reflecting the complexity and scale of the physical systems being modeled.

## Technologies Used

*   **Programming Languages:** C/C++
*   **Libraries/Frameworks:** Standard C libraries (`math.h`, `stdio.h`, `stdlib.h`, `time.h`, `string.h`)
*   **Analysis & Visualization:** Wolfram Mathematica (`.nb` notebooks)

## Getting Started

### Prerequisites

*   A C/C++ compiler (e.g., GCC)
*   `make` utility (for building the project)
*   (Optional) Wolfram Mathematica for analyzing notebooks.
*   (Optional) VMD (Visual Molecular Dynamics) or similar software for visualizing `.pdb` trajectory files.

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/srd-heterogeneous-catalysis-code2.git
    cd srd-heterogeneous-catalysis-code2
    ```
2.  **Build the project:**
    ```bash
    make all
    ```
    This will compile the source code and create an executable named `srd_simulation` in the `build/bin/` directory.

### Usage

1.  **Run a simulation:**
    ```bash
    ./build/bin/srd_simulation
    ```
    The simulation parameters can be adjusted by modifying `src/main.c`. Output data files (`.dat`, `.pdb`) will be generated in the `data/` directory.

2.  **Analyze results (using Mathematica):**
    Open the `.nb` files in the `notebooks/` directory with Wolfram Mathematica to analyze the generated data, plot profiles, and calculate derived quantities.

3.  **Visualize trajectories (using VMD):**
    Load the `trajectories.pdb` file (generated in the `data/` directory) into VMD to visualize particle movements and interactions.


\n## Compilation\n\nTo compile the code, run:\n```\nmake\n```
\n## Usage\n\nAfter compiling, run the simulation with:\n```\n./srd\n```
