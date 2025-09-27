# quantum-simulator

------------------------------------------------------------------------------------------------

QUANTUM SIMULATOR - Sept. 27 2025

------------------------------------------------------------------------------------------------

Overview:

This project is a toy quantum simulator which reads in a user-inputed initial 2-D particle state and an energetic potential function. It generates a corresponding Hamiltonian matrix and updates with discrete timesteps accordingly, outputing an mp4 video file. 

This is meant to be a tool to help visualize and conceptualize quantum mechanics.

------------------------------------------------------------------------------------------------

Prerequisites:

Python (numpy, scipy, matplotlib)

------------------------------------------------------------------------------------------------

To use:

1. Create python script with two defined functions: V(x,y), and wavefunction(x,y). 

2. Run "python quantum.py" from terminal.

3. Enter path to input python script.

Simulator will run script with a progress bar and output an mp4 video file to working directory.

