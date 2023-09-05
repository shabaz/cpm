# Speedy CPM

This is a Cellular Potts Model implementation that is provided along with the manuscript

"A parallelized cellular Potts model that enables simulations at tissue scale"

by Shabaz Sultan, Sapna Devi, Scott Mueller, and Johannes Textor.

The framework has been developed for benchmarking purposes and is highly optimized for single cell and low occupancy simulation scenarios. It is not intended to be a fully-fledged CPM framework, and user friendliness has not been a major concern during development. We provide this code in the hope that some of the underlying ideas can be integrated into other CPM frameworks.

The framework is implemented in C++, but we mainly use it from Python through a wrapper. Therefore we provide this implementation as a Python module.  


## System requirements


To compile the framework, the following software is required:

 * Python, version 3.0 or higher
 * numpy, version 1.24 or higher
 * setuptools, version 45.1 or higher


If you compile the framework through pip (recommended), then the requirements will be automatically installed.

## Installation instructions

Clone the repository, change into its directory and run

```
python -m pip install .
```

We have tested compilation on an Ubuntu 20.04.6 LTS Linux system (Python 3.9.12) and an Apple M1 MacBook Pro using macOS 12.6 "Monterey" (Python 3.9.12).

## Demo
 
The file `examples/sorting_cpu.py` contains a more detailed example implementation of the classic cell sorting simulation of Graner and Glazier (https://doi.org/10.1103/PhysRevLett.69.2013). Running this simulation should take only a few seconds. The script will produce a png file showing the final state of the simulation.

The file loads a pre-computed initial state saved as a numpy array (`examples/initial_state.npy`), so you need to change into the `examples` directory and run it from there. 10,000 simulation steps are run, this should take less than 1 minute on a reasonably recent system. 
