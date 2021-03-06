# A myofibre model for the study of uterine excitation-contraction dynamics
<!--[![DOI]()]()-->
<!--This is the code developed for our paper ____:-->
<!--(link)-->
This is the code developed for our paper "[A myofibre model for the study of uterine excitation-contraction dynamics](http://www.nature.com/articles/s41598-020-72562-x)"
<!--(Manuscript has been accepted for publication in Scientific Reports and will be published online on 1st Oct):\\-->
Uri Goldsztejn, Arye Nehorai

Washington University in St. Louis, 2020

*If you find this code useful in your research, please consider citing.*

## Content
* Overview
* Files
* Software requirements
* How to use
* Citation
* Contact

## Overview

Uterine contraction disorders can contribute to preterm births, slow progression of labor, and failure to initiaite labor. We developed a model to study uterine excitation-contraction dynamics at the muscle fiber level.

## Files

* *myofibre_code/* - The source code for the model and a makefile for the project.
<!---* *myofiber_lite_code/* The source code for a simplyfied version of our model. The computations are much faster and don't require the GNU scientific library. This version produces limited and approximated results, see the readme in this folder for further info.!--->
* *myofibre_analysis* - Matlab scripts to read the simulation results and recreate the figures in the paper.

* *output format description* - A description of the files created in every simulation.

## Software requirements

* This model needs to be compiled using a C++11 compiler, or newer.
* GNU Scientific Library (GSL) 2.4 or above.
This library and it's documentation can be found [here](https://www.gnu.org/software/gsl/).
<!---*The lite version does not use the GSL and can be run without it.!--->


#### GSL installation on Unix based systems

Instructions to install the GSL on unix-based systems can be found [here](https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/).

#### GSL installation on Windows based systems

A precompiled version of the GSL libraries for windows can be found [here](https://www.bruot.org/hp/libraries/).
These libraries can be linked to the project.

## How to use

After the GNU scientific library is installed, compile and link the program (see makefile in the code folder).
If running on unix-based systems:
Note that the makefile includes the directory in which the library is installed, this needs to be changed to the relevant directory. Additionally, before running the program, the library needs to be made visible to the operating system. This can be done with the following lines in bash:

$ LD_LIBRARY_PATH=/usr/local/lib

$ export LD_LIBRARY_PATH

$ ./myofibreModel

See [this guide](https://www.gnu.org/software/gsl/doc/html/usage.html) for more information.

Before running the program "myofibreModel", there should be an empy directory called "simulations" in the same directory as the program. The outputs of the simulations will be stored in that directory.

#### Simulation parameters
When the program is executed, it prompts the user to enter values for the simulation parameters.
Default values appear in brakets.

+ Simulation name: identifier for the simulation, the output file names will include this name.
+ IR: the intercellular resistivity. See the original article for more information on this parameter and it's definition.
+ Ionic channel modulations: The modulation is a multiplicative factor that alters the conductivity of the given ionic channel.
A modulation of 1, will have no effect on the ionic current, a modulation of 0 will delete this current from the model.

+ Option to change more parameters. Entering Y will prompt the user to modifiy more simulation parameters.
The following parameters can be changed if this option is chosen:
+ Simulation time: the time that will be simulated.
+ Number of cells: The number of cells concatenated serially in the model.
+ Time to turn on the stimulator: the time during the simulation at which the stimulator will be turn on.
+ Time to turn off the stimulator: the time during the simulation at which the stimulator will be turn off.
+ Step current: the current that will be injected to the first cell in the fiber while the stimulator is on.
+ Notes: the user can add notes about this simulation for future reference. These notes will be saved in the metadata file.

## Citation

Goldsztejn, U., Nehorai, A. A myofibre model for the study of uterine excitation-contraction dynamics. Sci Rep 10, 16221 (2020). https://doi.org/10.1038/s41598-020-72562-x

## Contact
Uri Goldsztejn
uri.goldsztejn@wustl.edu
