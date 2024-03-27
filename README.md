# Mixed-quantum-classical-dynamics

This Python codebase is designed for solving mixed quantum-classical dynamics, focusing initially on the relaxation of a two-level quantum system influenced by vacuum fluctuations within a cavity. However, the code can be readily adapted to tackle various other systems, regardless of their quantum system's dimensionality.

Structure Overview:

Main File: "int.py"

This file serves as the central control point, responsible for setting parameters and choosing integration options for both quantum and classical systems.


Dynamics Solver: "EOM.py"

Contains functions for solving the equations of motion (EOM) for both the classical and quantum systems. Specifically, it defines the CEOM function for classical dynamics, which corresponds to a Brownian oscillator, and the QEOM function for quantum dynamics, described by the time-dependent Schrödinger equation.


Quantum Hamiltonian Definition: "Hamiltonian.py"

Defines the Hamiltonian for the quantum system and manages interactions with the classical system.


Integration Methods:

For Quantum System:

Runge-Kutta 4th Order: Integ_QEOM

SciPy Integrators: Integ_QEOM2

For Classical System:

RK4: Integ_CEOM

Runge–Kutta–Fehlberg: Integ_CEOM2

RK8: Integ_CEOM3


Parallel Execution:
The code is designed to run in parallel using the POOL library, enhancing computational efficiency.


Customization:
To adapt the code for different systems, modifications are primarily made in the following files:

"int.py": Adjust parameters and integration options.

"Hamiltonian.py": Define the quantum system's Hamiltonian and manage interactions.

"EOM.py": Modify classical system dynamics if it deviates from a 1D Brownian oscillator.
