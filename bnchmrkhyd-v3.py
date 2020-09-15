#!/usr/bin/env python

import os, sys
import numpy as np
import shutil
from forcebalance.nifty import _exec

def parser():
    """
    Parse out info from system.txt, molecule.txt, lvlthry.txt, and solvent.txt

    Returns
    -------
    mold: dict
          Nested dictionary containing information about molecule. Each dictionary labeled by molecule ID
          donchg: donor charge
          accchg: acceptor charge
          donspn: donor spin mult
          accspn: acceptor spin mult
    pard: dict
          Nested dictionary containing information about levels of theory. Each dictionary labelled by level of theory ID.
          dft: DFT functional
          basis: basis set
          solvmod: solvent model
    sold: dict
          Nested dictionary containing information about solvents. Each dictionary labelled by solvent ID.
          sname: solvent name
          epsilon: dielectric constant of solvent
    sysd: dict
          Nested dictionary containing information about systems (combination of molecule and solvent). Each dictionary labelled by system ID
          molecule: molecule ID
          solvent: solvent ID
          hyd: experimental hydricity
          cit: citation number
    """  
    mold = {}

    for line in open(os.path.join('keys', 'molecules.txt')):
        #Each line that user wants parsed will start with molecule ID (M####), otherwise will start with # and will be skipped.
        if line.startswith('M'):
            #Make dictonaries inside mold, each labelled with their respective molecule IDs.
            info=line.split()
            label=info[0]
            mold[label]={}
            #Add donor and acceptor charges to dictionary
            mold[label]['donchg']=info[3]
            mold[label]['accchg']=info[5]
            #Add donor and acceptor spin mults to dictionary
            dspns0=info[4]
            aspns0=info[6]
            #Reformat spin mults into a list
            mold[label]['donspn']=dspns0.split(",")
            mold[label]['accspn']=aspns0.split(",")

    lvld = {}

    for line in open(os.path.join('keys', 'lvlthry.txt')):
        #Level of theory IDs start with P for parameters (P##)
        if line.startswith('P'):
            #Make dictionaries inside lvld, each labelled with their respective level of theory IDs.
            info=line.split()
            label=info[0]
            lvld[label]={}
            #Add information about the density functional, basis set, and solvent model to dictionary
            lvld[label]['dft']=info[1]
            lvld[label]['basis']=info[2]
            lvld[label]['solvmod']=info[3]

    sold = {}
            
    for line in open(os.path.join('keys', 'solvent.txt')):
        #Solvent IDs start with S (S##)
        if line.startswith('S'):
            #Make dictionaries inside sold, each labelled with their respective solvent IDs.
            info=line.split()
            label=info[0]
            #Add solvent name and dielectric constant to dictionary.
            sold[label]={}
            sold[label]['sname']=info[1]
            sold[label]['epsilon']=info[2]

    sysd = {}
    
    for line in open(os.path.join('keys', 'system.txt')):
        #System IDs start with Y (Y####)
        if line.startswith('Y'):
            #Make dictionaries inside sysd, each labelled with their respective system IDs.
            info=line.split()
            label=info[0]
            #Add molecule iD, solvent ID, experimental hydricity, and citation number to dictionary.
            sysd[label]={}
            sysd[label]['molecule']=info[1]
            sysd[label]['solvent']=info[3]
            sysd[label]['hyd']=info[5]
            sysd[label]['cit']=info[6]
            
    return mold, lvld, sold, sysd

def checkjobstatus(MID,SID,PID,sold,job,donacc,spin,jobtype):
    """
    Checks for prerequisites that need to be met before a job can be submitted.
    Parameters
    ----------
    MID, SID, PID: str
                   molecule, solvent, and level of theory ID
    sold: dict
          Nested dictionary created by parser, contains information about solvents
    job: str
         type of job to submit, e.g. minimize (optimization), frequencies (frequency analysis), etc.
    donacc: str
            type of structure (donor or acceptor)
    spin: str
          spin multiplicity
    jobtype: str
             type of job (geometry optimization or frequency analysis)
    Returns
    -------
    status: str
            Status of the requested job.
            -DNE: no jobs were submitted yet
            -queued: job is queued (running or pending)
            -done: job finished successfully
            -failed: job didn't finish correctly
            -optinc: user requested frequency analysis when optimization for that molecule is incomplete
    """
    solname = sold[SID]['sname']
    status = 'DNE'
    username = os.environ.get('USER')
    
    if jobtype == 'minimize':
        optfreq = 'opt'
    if jobtype == 'frequencies':
        optfreq = 'freq'
    
    #String containing the folder path
    jobdir = os.path.join('molecules', PID, MID, solname, spin, donacc, optfreq)
    
    #Check if folder contains running or queued job
    if os.path.exists(os.path.join(jobdir, 'submittedjob.txt')):
        #Read submittedjob.txt for job number
        subjobtxt = open(os.path.join(jobdir, 'submittedjob.txt'), "r")
        #Turn the first word of the first line into an integer
        jobID = int(list(subjobtxt.readlines())[0].split()[0])
        subjobtxt.close()
        #Execute "squeue -u username"
        jobqueue = _exec('squeue -u %s' % username)
        #Parse each line then check if there's a queued job that matches the job number in submittedjob.txt. If the job is queued, do nothing and move on, if not, check whether job crashed or finished successfully
        jqlist = []
        for i in range(len(jobqueue)):
            if i == 0:
                continue
            else:
                line = jobqueue[i]
                info = line.split()
                jobnum = int(info[0])
                jqlist.append(jobnum)
        for i in range(len(jqlist)):
            if jqlist[i] == jobID:
                status = 'queued'
            
        if status == 'DNE':
            if optfreq == 'opt':
                if os.path.exists(os.path.join(jobdir, 'end.xyz')):
                    status = 'done'
                else:
                    status = 'failed'
            if optfreq == 'freq':
                output = _exec('grep Gibbs run.out', cwd=jobdir)
                line = output[0]
                if line.startswith('Gibbs'):
                    status = 'done'
                else:
                    status = 'failed'
                    
    else:
        if optfreq == 'freq':
            optdir = os.path.join('molecules', PID, MID, solname, spin, donacc, 'opt')
            if not os.path.exists(os.path.join(optdir, 'end.xyz')):
                status = 'optinc'
        
    return status

                
def submit(MID,SID,PID,mold,lvld,sold,job,donacc,spin):
    """
    Creates necessary files then submits a TeraChem calculation job.
    
    Parameters
    ----------
    MID, SID, PID: str
                   molecule, solvent, and level of theory ID
    mold, lvld, sold: dict
                      Nested dictionaries created by parser
    job: str
         type of job to submit, e.g. minimize (optimization), frequencies (frequency analysis), etc.
    donacc: str
            type of structure (donor or acceptor)
    spin: str
          spin multiplicity
    """
    if donacc == "donor":
        name = "don"
        key = "donchg"
    else:
        name = "acc"
        key = "accchg"
    
    solname = sold[SID]['sname']
    
    #Submit an optimization job
    if job == "minimize":
        optdir = os.path.join('molecules', PID, MID, solname, spin, donacc, 'opt')
        #Create a folder
        if not os.path.exists(optdir): os.makedirs(optdir)
    
        #Write run.in file
        with open(os.path.join(optdir, 'run.in'), 'w') as f:
            f.write("run minimize\n")
            
            dft = lvld[PID]['dft']
            basis = lvld[PID]['basis']
            f.write("method %s\ndispersion yes\nbasis %s\n" % (dft, basis))
    
            #What to write if solvent model is PCM
            if lvld[PID]['solvmod'] == 'pcm':
                chg = mold[MID][key]
                f.write("charge %s\nspinmult %s\ncoordinates start.xyz\npcm cosmo\n" % (chg, spin))
    
                dec = sold[SID]['epsilon']
                f.write("epsilon %s\nconvthre 1e-6\nthreall 1e-18\ndftgrid 3\nthermotemp 298.15\nscf diis+a\n" % dec)
                f.write("\nmin_converge_e    1.0e-7\nmin_converge_grms 3.0e-5\nmin_converge_gmax 4.5e-5\nmin_converge_drms 1.2e-4\nmin_converge_dmax 1.8e-4")
    
            #What to write if using "gas phase"
            if lvld[PID]['solvmod'] == 'gas':
                chg = mold[MID][key]
                f.write("charge %s\nspinmult %s\ncoordinates start.xyz\n" % (chg, spin))
    
                f.write("convthre 1e-6\nthreall 1e-18\ndftgrid 3\nthermotemp 298.15\nscf diis+a\n")
                f.write("\nmin_converge_e    1.0e-7\nmin_converge_grms 3.0e-5\nmin_converge_gmax 4.5e-5\nmin_converge_drms 1.2e-4\nmin_converge_dmax 1.8e-4")
                
        #Copy over starting structure and rename as start.xyz        
        shutil.copy(os.path.join('startxyz', "%s-%s.xyz" % (MID, name)), optdir)
        os.rename(os.path.join(optdir, "%s-%s.xyz" % (MID, name)), os.path.join(optdir, "start.xyz"))

        #Write batch.script (doing this so that time limit can be set to 7 days, not 2 days (the default)
        #Read start.xyz file to check for number of atoms
        startcoord = open(os.path.join(optdir, "start.xyz"), "r")
        #Turn the first word of the first line into an integer
        atomnum = int(list(startcoord.readlines())[0].split()[0])
        startcoord.close()
        with open(os.path.join(optdir, 'batch.script'), 'w') as g:
            #Depending on the number of atoms, set the number of GPUs to either 1, 2, or 4.
            if atomnum < 30:
                g.write("#!/bin/bash -l\n#SBATCH -p gpu\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 2\n#SBATCH --gres=gpu:1\n#SBATCH --mem=8000\n#SBATCH -J tera\n#SBATCH -t 7-00:00:00\n")
            if 30 <= atomnum < 50:
                g.write("#!/bin/bash -l\n#SBATCH -p gpu\n#SBATCH -N 1\n#SBATCH -n 2\n#SBATCH -c 2\n#SBATCH --gres=gpu:2\n#SBATCH --mem=16000\n#SBATCH -J tera\n#SBATCH -t 7-00:00:00\n")
            if atomnum >= 50:
                g.write("#!/bin/bash -l\n#SBATCH -p gpu\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH -c 2\n#SBATCH --gres=gpu:4\n#SBATCH --mem=32000\n#SBATCH -J tera\n#SBATCH -t 7-00:00:00\n")
            g.write("\n#SBATCH --no-requeue\n")
            g.write('\n# Record job info\necho -e "$SLURM_JOB_ID  $HOSTNAME  $(pwd)" >> ~/.myslurmhistory\n')
            g.write("\nmodule load intel cuda/9.0\nexport CUDA_CACHE_PATH=/scratch/$USER/.nv/ComputeCache\nexport TeraChem=/home/leeping/opt/terachem/current\nexport PATH=$TeraChem/bin:$PATH\nexport LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH\n\n")
            g.write("\nterachem  run.in &> run.out\n")
            g.write("\nna=$(head -1 scr/optim.xyz) && tail -$((na+2)) scr/optim.xyz > end.xyz")
        
        #Submit optimization job
        submit = _exec('sbatch batch.script', cwd=optdir)
        submitstring = submit[0]
        submitsplt = submitstring.split()
        jobnum = submitsplt[3]
        #Save job number to txt file
        with open(os.path.join(optdir, 'submittedjob.txt'), 'w') as h:
            h.write(jobnum)
    
    #Submit a frequency job
    else:  
        #String containing the folder path
        freqdir = os.path.join('molecules', PID, MID, solname, spin, donacc, 'freq')
        #Create a folder
        if not os.path.exists(freqdir): os.makedirs(freqdir)
        
        #Copy over run.in file from opt, change "minimize" to "frequencies"
        optdir = os.path.join('molecules', PID, MID, solname, spin, donacc, 'opt')
        shutil.copy(os.path.join(optdir, 'run.in'), freqdir)
        
        runin = open(os.path.join(freqdir, 'run.in'), "r")
        lines = runin.readlines()
        lines[0] = "run frequencies\n"

        runin = open(os.path.join(freqdir, 'run.in'), "w")
        runin.writelines(lines)
        runin.close()
        
        #Copy over end.xyz file from opt, rename as start.xyz
        shutil.copy(os.path.join(optdir, 'end.xyz'), freqdir)
        os.rename(os.path.join(freqdir, 'end.xyz'), os.path.join(freqdir, 'start.xyz'))
        
        #Write batch.script (doing this so that time limit can be set to 7 days, not 2 days (the default)
        #Read start.xyz file to check for number of atoms
        startcoord = open(os.path.join(freqdir, "start.xyz"), "r")
        #Turn the first word of the first line into an integer
        atomnum = int(list(startcoord.readlines())[0].split()[0])
        startcoord.close()
        with open(os.path.join(freqdir, 'batch.script'), 'w') as g:
            #Depending on the number of atoms, set the number of GPUs to either 1, 2, or 4.
            if atomnum < 30:
                g.write("#!/bin/bash -l\n#SBATCH -p gpu\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 2\n#SBATCH --gres=gpu:1\n#SBATCH --mem=8000\n#SBATCH -J tera\n#SBATCH -t 7-00:00:00\n")
            if 30 <= atomnum < 50:
                g.write("#!/bin/bash -l\n#SBATCH -p gpu\n#SBATCH -N 1\n#SBATCH -n 2\n#SBATCH -c 2\n#SBATCH --gres=gpu:2\n#SBATCH --mem=16000\n#SBATCH -J tera\n#SBATCH -t 7-00:00:00\n")
            if atomnum >= 50:
                g.write("#!/bin/bash -l\n#SBATCH -p gpu\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH -c 2\n#SBATCH --gres=gpu:4\n#SBATCH --mem=32000\n#SBATCH -J tera\n#SBATCH -t 7-00:00:00\n")
            g.write("\n#SBATCH --no-requeue\n")
            g.write('\n# Record job info\necho -e "$SLURM_JOB_ID  $HOSTNAME  $(pwd)" >> ~/.myslurmhistory\n')
            g.write("\nmodule load intel cuda/9.0\nexport CUDA_CACHE_PATH=/scratch/$USER/.nv/ComputeCache\nexport TeraChem=/home/leeping/opt/terachem/current\nexport PATH=$TeraChem/bin:$PATH\nexport LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH\n\n")
            g.write("\nterachem  run.in &> run.out")
        
        #Submit frequency analysis job
        submit = _exec('sbatch batch.script', cwd=freqdir)
        submitstring = submit[0]
        submitsplt = submitstring.split()
        jobnum = submitsplt[3]
        #Save job number to txt file
        with open(os.path.join(freqdir, 'submittedjob.txt'), 'w') as h:
            h.write(jobnum)

def gethyd(MID,SID,PID,mold,sold,donacc):
    """
    Reads run.out file to get the Gibbs free energy of the optimal structure for a molecule.
    
    Parameters
    ----------
    MID, SID, PID: str
                   molecule, solvent, and level of theory ID
    mold, sold: dict
                nested dictionaries created by parser
    donacc: str
            type of structure (donor or acceptor)
    
    Returns
    -------
    minsm: str
           spin multiplicity of the structure with the minimum free energy
    minG: str
          minimum Gibbs free energy
    """
    if donacc == "donor":
        key = "donchg"
    else:
        key = "accchg"
    solname = sold[SID]['sname']
    Edict= {}

    #Loop over all spin multiplicities for the molecule of interest
    for sm in mold[MID][key]:
        freqdir = os.path.join('molecules', PID, MID, solname, sm, donacc, 'freq')
        energy = _exec('grep Gibbs run.out', cwd=freqdir)
        #Take the last line and parse to get energy
        Esplit = energy.split()
        Gibbs = Esplit[-2]
        #Store to dictionary
        Edict[sm] = Gibbs
    #Search dictionary for spin mult with lowest energy
    minsm = min(Edict, key=Edict.get)
    minG = Edict[minsm]
    
    return minsm, minG

def main():
    #Create dictionaries
    mold, lvld, sold, sysd = parser()
    #Ask user to provide level of theory and job type
    PID = str(input('Enter level of theory ID: \n'))
    jobtype = str(input('Enter job type: \n'))
    
    if jobtype == "minimize":
        for YID in sysd:
            MID = sysd[YID]['molecule']
            SID = sysd[YID]['solvent']
            #Using spin mults as key since number of files to be created depends on number of possible spin mults
            for key in mold[MID]['donspn']:
                stat = checkjobstatus(MID,SID,PID,sold,jobtype,'donor',key,jobtype)
                if stat == "DNE":
                    submit(MID,SID,PID,mold,lvld,sold,jobtype,'donor',key)
                if stat == "queued":
                    print("Geometry optimization for %s donor already queued." % MID)
                if stat == "failed":
                    print("Geometry optimization for %s donor failed." % MID)
                if stat == "done":
                    print("Geometry optimization done for %s donor, time to submit frequency jobs." % MID)
            #Same thing for acceptor structures
            for key in mold[MID]['accspn']:
                stat = checkjobstatus(MID,SID,PID,sold,jobtype,'acceptor',key,jobtype)
                if stat == "DNE":
                    submit(MID,SID,PID,mold,lvld,sold,jobtype,'acceptor',key)
                if stat == "queued":
                    print("Geometry optimization for %s acceptor already queued." % MID)
                if stat == "failed":
                    print("Geometry optimization for %s acceptor failed." % MID)
                if stat == "done":
                    print("Geometry optimization done for %s acceptor, time to submit frequency jobs." % MID)
                    
    if jobtype == "frequencies":
        for YID in sysd:
            MID = sysd[YID]['molecule']
            SID = sysd[YID]['solvent']
            #Using spin mults as key since number of files to be created depends on number of possible spin mults
            for key in mold[MID]['donspn']:
                stat = checkjobstatus(MID,SID,PID,sold,jobtype,'donor',key,jobtype)
                if stat == "DNE":
                    submit(MID,SID,PID,mold,lvld,sold,jobtype,'donor',key)
                if stat == "queued":
                    print("Frequency analysis for %s donor already queued." % MID)
                if stat == "failed":
                    print("Frequency analysis for %s donor failed." % MID)
                if stat == "done":
                    print("Frequency analysis for %s donor done, time to calculate hydricity." % MID)
                if stat == "optinc":
                    print("Geometry optimization for %s donor incomplete." % MID)

            #Same thing for acceptor structures
            for key in mold[MID]['accspn']:
                stat = checkjobstatus(MID,SID,PID,sold,jobtype,'acceptor',key,jobtype)
                if stat == "DNE":
                    submit(MID,SID,PID,mold,lvld,sold,jobtype,'acceptor',key)
                if stat == "queued":
                    print("Frequency analysis for %s acceptor already queued." % MID)
                if stat == "failed":
                    print("Frequency analysis for %s acceptor failed." % MID)
                if stat == "done":
                    print("Frequency analysis for %s acceptor done, time to calculate hydricity." % MID)
                if stat == "optinc":
                    print("Geometry optimization for %s acceptor incomplete." % MID)
                    
    #If freq analysis done for all spin mults, calculate hydricity
    if jobtype == "hydricity":
        for YID in sysd:
            MID = sysd[YID]['molecule']
            SID = sysd[YID]['solvent']
            donminspn, donfreeE = gethyd(MID,SID,PID,mold,sold,'donor')
            accminspn, accfreeE = gethyd(MID,SID,PID,mold,sold,'acceptor')
            delG_HHR = accfreeE - donfreeE
            if SID == 'S03': #water
                hydricity = delG_HHR - 420.7
            else: #acetonitrile, benzonitrile, DMSO
                hydricity = delG_HHR - 419.34
            #Save results to dat file (YID, donor charge, donor spin, acceptor charge, acceptor spin, hydricity)
            with open(os.path.join('data', 'data-%s.dat' % PID), 'a+') as f:
                f.write(YID + "\t" + mold[MID]['donchg'] + "\t" + donminspn + "\t" + mold[MID]['accchg'] + "\t" + accminspn + "\t" + hydricity + "\n")
    
if __name__ == "__main__":
    main()
