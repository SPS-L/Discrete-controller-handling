# IBM, Interpolation-based method for solving digital controllers alongside power system DAEs
The codes for the case studies of the paper "Fast and Accurate Simulation of Smart Digital Controllers in Power System Dynamic Studies" published in Sustainable Energy, Grids, and Networks are included.
## How to use
The files of the project are paired together. For all the test systems, please run the SRM file first, then the SSM (simplified simulation method), and finally IBM. For the SMIB and Kundur systems, Please import the fuzzy controller files into the workspace before running the codes.

### Integration method settings

This method uses a predictor-corrector approach for integration. Currently, regarding the integration method, two options are available:
1. Predictor: Second-order Adams-Bashford  Corrector: Second-order Adams-Moulton
1. Predictor: Second-order Adams-Bashford  Corrector: Second-order BDF

In order to use the first option simply set the parameters: cp=5 & cc=3
In order to use the second option simply set the parameters: cp=5 & cc=2

### Quantization method settings

Two quantization methods are also available. In order to use the uniform quantization set cq=1, and for the nonuniform quantization set cq=2. Moreover, the number of bits used for the quantization can be set inside the function called *quz*.

### Interpolator settings
There are different interpolators available to use for the IBM method, two explicit interpolators and three implicit ones. To change the explicit and implicit interpolators inside the file __*IBM*__, you should comment the current and uncomment the one you wish to choose inside the functions "*predictorg*" and "*funcg*", respectively. Please note that only one of the interpolators must be uncommented.

### Simulation settings

The time of the simulation can be set by the parameter *tsim*. The parameters *h* *hmin*, and *hmax* determine the first time step size, the minimum limit of the time step size, and the maximum limit of the time step size.

### Controller settings

The parameters *G* and *T* control the gain and sampling rate of the integral controller. Please refer to the AVR and Governor controller functions to modify their parameters.

### System settings

The eigenvalues of the system (*a* and *b*) can be set inside of the function *basedmatrix*.

## Results

The figures will be made automatically. Set the variable *save* to 1 in the *IBM* files if you wish to save the figures as PDF files. Besides, the overall number of iterations is saved under the variable called *iterationsc*.
