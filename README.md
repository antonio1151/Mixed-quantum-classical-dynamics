# Mixed-quantum-classical-dynamics

This is a Python code for solving mixed quantum-classical dynamics.
This code has been done for solving the case of relaxation of a two-level system due to the vacuum fluctuations in a cavity.
However, the code can be easily modified to solve any other system regardless of the dimension of the quantum system.

The main file is "int.py". In this file is setting the parameters for the dynamics and the options for the intergrator used for the quantum and/or classical system. 
The dynamics is done by the file "EOM.py". The equation of motion for the classical system (CEOM function) corresponds to a Brownian oscillator. The dynamics of the quantum system 
is decribed by the time dependent Schridinger equation (QEOM function). 
The integration of the dynamics of the quantum system can be done either by Runge-Kutta of order 4th (Integ_QEOM function) or by using the integradors of scipy (Integ_QEOM2 function).
The integration of the dynamics of the classical system can be done by three options of different Runge-Kutta algorithms:
      1. Integ_CEOM:  is RK4
      2. Integ_CEOM2:  is the Runge–Kutta–Fehlberg method
      3. Integ_CEOM3: is RK8
The Hamiltonian for the Quantum system is set in the file "Hamiltonian.py". Also, in this file is the the interaction with the classical system.

The code run in parallel using the library POOL. 

In order to modify the code to adapt it for different problem, you should modify some parts on "int.py", set the Hamiltonian of the quantum system in "Hamiltonian.py", and the equation of motion for the classical system in "EOM.py" if the classical system is not a 1D Brownian oscillator. 
