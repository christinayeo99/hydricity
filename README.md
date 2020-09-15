*****bnchmrkhyd-v3.py*****

General Description
-------------------
This code will loop through all "systems" (combination of solute and solvent) for one "level of theory" (combination of DFT-basis set-solvent model) to:
-check for the job status
-create directories and files
-submit jobs
-calculate hydricities from output files and store the numbers to a dat file

Organization
------------
For this code to work properly, txt files to be parsed should be stored in a folder called 'keys' (look at the example folder to see how these files are formatted) and starting
structures should be stored in a folder called 'startxyz'. All the TeraChem input and output files will be under a folder called 'molecules', and dat files will be under a
folder called 'data'.

Inputs
------
1) Level of Theory ID: should start with "P"
2) Job Type:
  -"minimize": Geometry optimization
  -"frequencies": Frequency Analysis
  -"hydricity": Calculate hydricities from output files and store to dat file

Dependencies
------------
Can be installed by typing "conda install ______" on terminal
-os
-sys
-numpy
-shutil
-forcebalance