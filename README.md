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

The V2/ folder contains the code used in Zinck et al, 2024 (submitted)

# How to apply BURGEE V2
All code is provided in javascript which can be copy-pasted and run on the Google Earth Engine directly

1. Perform the intitial co-registration (V2/CoRegistration)
2. Perform the second co-registration (V2/overlapCoRegistration)
3. Add a band to all strips with firn air content (V2/addFACband)
4. Perform the lagrangian displacement (V2/lagrangianDisplacementOfStrips)
5. Calculate the basal melt rate (V2/basalMeltRate)

Local GEE paths should be adjusted for saving and loading assets in-between scripts. All scripts are set to focus on the area around Totten Ice Shelf as a demo and can be moved to other ice shelves as explained in comments in the code.

The computational time depends on the amount of REMA strips and the size of the ice shelf/area. It is not possible to run all years in one go in step 1,2 and 4 above due to the computational limit on the GEE, and the allowed time window can therefore be used to adjust this. Running all of the above scripts for Totten only will take around a week in total on the GEE, with the lagrangian displacement being the bottleneck. 


