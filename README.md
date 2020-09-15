*****bnchmrkhyd-v3.py*****

General Description
-------------------
This code will loop through all "systems" (combination of solute and solvent) for one "level of theory" (combination of DFT-basis set-solvent model) to:
1) Check for the job status
2) Create directories and files
3) Submit jobs
4) Calculate hydricities from output files and store the numbers to a dat file

Organization
------------
For this code to work properly, txt files to be parsed should be stored in a folder called 'keys' (look at the example folder to see how these files are formatted) and starting
structures should be stored in a folder called 'startxyz'. All the TeraChem input and output files will be under a folder called 'molecules', and dat files will be under a
folder called 'data'.

Inputs
------
1) Level of Theory ID: should start with "P"
2) Job Type:
  a. "minimize": Geometry optimization
  b. "frequencies": Frequency Analysis
  c. "hydricity": Calculate hydricities from output files and store to dat file

Dependencies
------------
Can be installed by typing "conda install ______" on terminal
1) os
2) sys
3) numpy
4) shutil
5) forcebalance
