# bnchmrkhyd-v3.py

## General Description
This code will loop through all "systems" (combination of solute and solvent) for one "level of theory" (combination of DFT-basis set-solvent model) to:
* Check for the job status
* Create directories and files
* Submit jobs
* Calculate hydricities from output files and store the numbers to a dat file

## Organization
For this code to work properly, txt files to be parsed should be stored in a folder called 'keys' and starting structures should be stored in a folder called 'startxyz'. All the TeraChem input and output files will be under a folder called 'molecules', and dat files will be under a folder called 'data'. Look at the example folder to see how these files are formatted.

## Inputs
* Level of Theory ID: should start with "P"
* Job Type:
  * "minimize": Geometry optimization
  * "frequencies": Frequency Analysis
  * "hydricity": Calculate hydricities from output files and store to dat file

## Dependencies
Can be installed by typing "conda install ______" on terminal
* os
* sys
* numpy
* shutil
* forcebalance
