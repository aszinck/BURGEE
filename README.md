Created by Ann-Sofie P. Zinck (a.p.zinck@tudelft.nl | @aszinck)

# Intro to BURGEE
This repository contains the code for the Basal melt rates Using Rema and Google Earth Engine (BURGEE) method. This is a method to obtain basal melt rates at high spatial resolution.

I highly recommend new uses to Google Earth Enigne to first have a look at [Google's Get Started with Earth Engine page](https://developers.google.com/earth-engine/guides/getstarted). It should, furthermore, be noted that using GEE, and thereby also BURGEE, requires an GEE account.

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
All code is provided in javascript which can be copy-pasted and run on the Google Earth Engine directly or just be clicking the links below

1. Perform the intitial co-registration (V2/CoRegistration) [GEE link](https://code.earthengine.google.com/15025af55c7e1bcb231e580d9825f239)
2. Perform the second co-registration (V2/overlapCoRegistration) [GEE link](https://code.earthengine.google.com/be1c68c298d2fe2db9f8cb3b9379205d)
3. Add a band to all strips with firn air content (V2/addFACband) [GEE link](https://code.earthengine.google.com/ac18b2889574c1919cb56b2b1a02674c)
4. Perform the lagrangian displacement (V2/lagrangianDisplacementOfStrips) [GEE link](https://code.earthengine.google.com/3a2fc8bd9f5cc26bd2fe5cad3abcadc8)
5. Calculate the basal melt rate (V2/basalMeltRate) [GEE link](https://code.earthengine.google.com/ba58d1b36b145060280a3422e91e9542)

Local GEE paths should be adjusted for saving and loading assets in-between scripts. All scripts are set to focus on the area around Totten Ice Shelf as a demo and can be moved to other ice shelves as explained in comments in the code.

The computational time depends on the amount of REMA strips and the size of the ice shelf/area. It is not possible to run all years in one go in step 1,2 and 4 above due to the computational limit on the GEE, and the allowed time window can therefore be used to adjust this. Running all of the above scripts for Totten only will take around a week in total on the GEE, with the lagrangian displacement being the bottleneck. 


