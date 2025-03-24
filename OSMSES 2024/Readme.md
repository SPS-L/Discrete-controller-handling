# OSMSES2024

This project uses a novel interpolation-based method to solve a set of DAEs alongside the difference equations of the digital controller. you may find the simulation files of the paper "A Modeilica IBM Implementation for Fast Simulation of Digital Controllers in Power Systems"

## How to use
There are 2 case studies named __*Integral Controller*__, and __*SMIB*__ for simulation of the classical integral controller, and the SMIB system, respectively, used in the paper.

### How to use

Each folder has 2 main files modeling the system in the SRM and IBM approaches. The files regarding the controller models can also be found in the folders for each case study. First, open the main file, then load the controller models as required. Please notice that the equations of the controllers are the same for both SRM and IBM methods since the goal is to provide a new method for solving power systems with digital controllers and not remodeling the components of the system.

Regarding the simulation settings, choice of solvers, and possible results refer to the paper published in OSMSES 2024 cited below:

