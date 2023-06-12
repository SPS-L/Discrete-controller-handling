# Modelica implementation of IBM
This project uses a novel interpolation based method to solve a state-space model with two differential equations controlled by an integral controller modeolled inside Modelica.

## How to use
The files of the project are made using Dymola. There are 4 separate folders named __*IBM Original*__, __*IBM Error*__, __*IBM Error 2*__, and __*IBM D Equation*__ for simulation of the classical integral controller example used in the paper. The first one produces the original results with a separate derivative block. The second one has the derivative explicitly as a separate equation that leads to an error. The third one is similar to the fourth one leading to the same error but having a different equation for calculation of derivative (proposed in the meeting). The fourth one contains an equation that tries to estimate the derivative however its accuracy is limited to the fixed time step size.
The file named __*IBM*__ in each folder has the custom Block the of controller of the specific folder. Please make sure you don't overwrite them on each other.
You don't need to open the files named __*IBM*__ in Dymola, separately.

### Interpolator settings
There are two intorpolation equations integrated into the custom made controller block. The first one that is commented out is the first order interpolator and the second one in use is second-order. Please comment in or out all 5 equations for all intermediate points based on your choice.

### Simulation settings

We have assumed a fixed time step size containing five intermediate points. Currently, it is assumed that the time step size is equal to 0.05 with the controller sampling time equal to 0.01. Therefore, please use a fixed time step solver like *RK2* with the exact size equal to 0.05.
