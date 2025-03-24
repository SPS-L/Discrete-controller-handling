# Light IBM, handling computationally heavy digital controllers
The codes for the case studies of the paper "Handling of Computationally Demanding Digital Controllers in Power System Dynamic Simulations" published in SyNERGY MED 2024 are included.
## How to use
The files of the project are paired together. For all the test systems, please run the SRM file first, then the IBM, and finally the Heavy IBM.

### Integration method settings

This method uses a predictor-corrector approach for integration. Currently, regarding the integration method, two options are available:
1. Predictor: Second-order Adams-Bashford  Corrector: Second-order Adams-Moulton
1. Predictor: Second-order Adams-Bashford  Corrector: Second-order BDF

In order to use the first option simply set the parameters: cp=5 & cc=3
In order to use the second option simply set the parameters: cp=5 & cc=2

### Quantization method settings

Two quantization methods are also available. In order to use the uniform quantization set cq=1, and for the nonuniform quantization set cq=2. Moreover, the number of bits used for the quantization can be set inside the function called *quz*.

### Interpolator settings
There are different interpolators available to use for the IBM method, two explicit interpolators and three implicit ones. To change the explicit and implicit interpolators inside the file __*IBM*__, you should comment the current and decomment the one you wish to choose inside the functions "*predictorg*" and "*funcg*", respectively. Please note that only one of the interpolators must be decomented.

### Simulation settings

The time of the simulation can be set by the parameter *tsim*. The parameters *h* *hmin*, and *hmax* determine the first time step size, the minimum limit of the time step size, and the maximum limit of the time step size.

### Controller settings

Each folder has the same test system with a different controller solved with 3 different methods. For the equation-based controller you just need to run the files in order. For the fuzzy controller, please first import the fuzzy controller to the workspace. For the machine-learning controller, first run the SRM file for the equation-based controller, then run the file "*ANNAVR*" to train the controller. Finally, run the rest of the files in the order.

## Results

The figures will be made automatically. Set the variable *save* to 1 in the *decoupled* files if you wish to save the figures as PDF files. Besides, the overall number of iterations is saved under the variable called *iterationsc*.
