# Discrete-controller-handling

## How to use
The files of the project are paired together. In each test case, please run the __*Analog*__ file first, then the __*Digital*__ one. In the case of the 3-bus system, both the analog and digital files have 2 separate files solving the system with a constant and variable time-step method.

### Integration method settings

This method uses a predictor-corrector approach for integration. Currently, regarding the integration method, two options are available:
1. Predictor: Second-order Adams-Bashford  Corrector: Second-order Adams-Moulton
1. Predictor: Second-order Adams-Bashford  Corrector: Second-order BDF

In order to use the first option simply set the parameters: cp=5 & cc=3
In order to use the second option simply set the parameters: cp=5 & cc=2

### Quantization method settings

Two quantization methods are also available. In order to use the uniform quantization set cq=1, and for the nonuniform quantization set cq=2. Moreover, the number of bits used for the quantization can be set inside the function called *quz*.

### Simulation settings

The time of the simulation can be set by the parameter *tsim*. The parameters *h* *hmin*, and *hmax* determine the first time step size, the minimum limit of the time step size, and the maximum limit of the time step size.

### Controller settings

The parameters *G* and *T* control the gain and sampling rate of the PI controller.

### System settings

The eigenvalues of the system (*a* and *b*) can be set inside of the function *basedmatrix*.

## Results

The figures will be made automatically. Besides, the overall number of iterations is saved under the variable called *iterationsc*.
