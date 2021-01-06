# bnchmrkhyd.py

## General Description
This code will loop through all "systems" (combination of solute and solvent) for one "level of theory" (combination of DFT-basis set-solvent model) to:
* Check for the job status
* Create directories and files
* Submit TeraChem jobs
* Calculate hydricities from output files
* Store results to a dat file
* Make scatter plots with a linear trendline

## Organization
For this code to work properly, txt files to be parsed should be stored in a folder called 'keys' and starting structures should be stored in a folder called 'startxyz'. All the TeraChem input and output files will be under a folder called 'molecules', and dat files will be under a folder called 'data'. Look at the example folder for guidance. The example folder contains:
* txt files that the code will read, parse, and make dictionaries from
* Starting structures for formate (donor) and CO2 (acceptor)
* TeraChem input/output files and data for formate calculated using ub3lyp/6-31gs_ldz/pcm
* Calculated vs experimental hydricities for molecules in acetonitrile

## Inputs
* Level of Theory ID: should start with "P"
* Job Type:
  * "minimize": Geometry optimization
  * "frequencies": Frequency Analysis
  * "hydricity": Calculate hydricities from output files, store to dat file, make scatter plots

## Dependencies
Can be installed by typing "conda install ______" on terminal
* os
* sys
* numpy
* shutil
* forcebalance

## Important Notes
* The "Calculated Hydricity" stored to the dat files do NOT include the free energy of the hydride anion in solvent, meaning it is technically the free energy of the "hydricity half reaction" (HHR), i.e. G(acceptor) - G(donor).
* If you have 2 or more data points for a given solvent, this code will calculate an "average y-intercept" (see the code for what this means), and shift the calculated hydricities accordingly. This shift corresponds to the free energy of the hydride anion in solvent. The "average y-intercepts" will be stored to a dat file as well.
