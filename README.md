Created by Ann-Sofie P. Zinck (a.p.zinck@tudelft.nl | @aszinck)

# Intro to BURGEE
This repository contains the code for the Basal melt rates Using Rema and Google Earth Engine (BURGEE) method. This is a method to obtain basal melt rates at high spatial resolution and is accompanied with a paper to be submitted soon.

# About BURGEE
BURGEE creates high resolution co-registered DEMs using the Reference Elevation Model of Antarctica (REMA) in combination with CryoSat-2. The DEMs can be used to to calculate surface elevation changes in both a Eulerian and Lagrangian framework, from which the latter is used to obtain high resolution basal melt rates. 
So far BURGEE has been tested on the Dotson Ice Shelf.

# Structure of the repo
The Dotson/ folder contains the code used in [Zinck et al, 2023](https://doi.org/10.5194/tc-17-3785-2023).

The Dotson/ folder contains the scripts to
1. Perform the initial co-registration.
2. Perform the second co-registration.
3. Calculate the Eulerian elevation change.
4. Calculate the Lagrangian elevation change and the basal melt rate.

The Dotson/modules/ folder contains two scripts with functions used for the co-registration and the Eulerian elevation change, respectively.

# How to apply BURGEE
All code is provided in javascript which can be copy-pasted and run on the Google Earth Engine
