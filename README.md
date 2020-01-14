# A myofiber model for the study of uterine excitation-contraction dynamics
[![DOI]()]()

This is the code developed for our paper in ____:
(link)

## Content
* Overview
* Files
* Software requirements
* How to use

## Overview

Uterine contraction disorders can contribute to preterm births, slow progression of labor, and failure to initiaite labor. We developed a model to study uterine excitation-contraction dynamics at the muscle fiber level.

## Files

* *myofiber_code/* - The source code for the model
* *myofiber_lite_code/* The source code for a simplyfied version of our model. The computations are much faster and don't require the GNU scientific library. This version produces limited and approximated results, see the readme in this folder for further info.
* *myofiber_analysis/* - Matlab scripts to read the simulation results and recreate the figures in the paper

## Software requirements

* This model needs to be compiled using a C++11 compiler, or newer.
* GNU Scientific Library (GSL) 2.4 or above.
This library and it's documentation can be found [here](https://www.gnu.org/software/gsl/).
The lite version does not use the GSL and can be run without it.


#### GSL installation on Unix based systems

Instructions to install the GSL on unix-based systems can be found [here](https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/).

#### GSL installation on Windows based systems



