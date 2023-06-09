# Modelica implementation of IBM
This project uses a novel interpolation based method to solve discontinuous systems modeled by differential equations.
## How to use
The files of the project are made using Dymola. There are two files named __*IBM equation*__ and __*IBM error*__ for simulation of the classical integral controller example used in the paper. The first one has an estimate derivateive equation inside the IBM custom made controller block while the second one has it explicitely that leads to an error.

### Interpolator settings
There are two intorpolatation equations integrated in the custom made controller block. The first one that is commented out is the first order interpolator and the second one in use is second-order. Please comment in or out all 5 equations for all intermediate points based on your choice.

### Simulation settings

We have assumed a fixed time step size containing five intermediate points. Currently, it is assumed that the time step size is equal to 0.05 with the controller sampling time equal to 0.01. Therefore, please use a fixed time step solver like *RK2* with the exact size equal to 0.05.
