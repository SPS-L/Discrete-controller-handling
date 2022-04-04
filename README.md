# Discrete-controller-handling
This project uses a novel interpolation based method to solve discontinuous systems modeled by differential equations.
## How to use
The files of the project are paired together. In the case of single controller, please use files __*Reference*__ and __*IBM*__. You can use the files __*ReferenceMultiController*__ & __*IBMMultiController*__ for the simulation of multiple controllers. Also, Please set the settings described below for both files of every pair. After the settings are set, first the reference code must be run then the IBM code.

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

The parameters *G* and *T* control the gain and sampling rate of the PI controller.

### System settings

The eigenvalues of the system (*a* and *b*) can be set inside of the function *basedmatrix*.

## Results

The figures will be made automatically. Besides, the overall number of iterations is saved under the variable called *iterationsc*.


