#!/usr/bin/env python

import os#, sys
import numpy as np
import IPython
import shutil
from forcebalance.nifty import _exec
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from forcebalance.molecule import Molecule

def parser(project, molfile):
    """
    Parse out info from system.txt, molecule.txt, lvlthry.txt, and solvent.txt

    Parameters
    ----------
    project: str
             Determines which class of molecules to create dictionaries for (just organic,
             just organometallic, or both)
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

    for line in open(molfile):
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
            #Add molecule type (organic or organometallic) to dictionary
            mold[label]['type']=info[7]

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
            
    modd = {}
    
    for line in open(os.path.join('keys', 'models.txt')):
        if line.startswith('D'):
            info=line.split()
            label=info[0]
            modd[label]={}
            modd[label]['optfreq']=info[1]
            if len(info)==3:
                modd[label]['sglpt']=info[2]

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
    
    if project == 'organic':
        sysfile = 'system-org.txt'
    if project == 'organometallic':
        sysfile = 'system-met.txt'
    if project == 'both':
        sysfile = 'system.txt'
    for line in open(os.path.join('keys', sysfile)):
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
            
    return mold, lvld, sold, sysd, modd

def checkjobstatus(MID,SID,PIDof,PIDsp,sold,job,donacc,spin):
    """
    Checks for prerequisites that need to be met before a job can be submitted.
    Parameters
    ----------
    MID, SID, PID: str
                   molecule, solvent, and level of theory ID
    sold: dict
          Nested dictionary created by parser, contains information about solvents
    job: str
         type of job to submit
    donacc: str
            type of structure (donor or acceptor)
    spin: str
          spin multiplicity
    Returns
    -------
    status: str
            Status of the requested job.
            -DNE: no jobs were submitted yet
            -queued: job is queued (running or pending)
            -done: job finished successfully
            -failed: job didn't finish correctly
            -optinc: user requested frequency analysis or single point calculation when
            optimization for that molecule is incomplete
    """
    solname = sold[SID]['sname']
    status = 'DNE'
    username = os.environ.get('USER')
    
    if job == 'minimize':
        jobword = 'opt'
        PID=PIDof
    if job == 'frequencies':
        jobword = 'freq'
        PID=PIDof
    if job == 'energy':
        jobword = 'sglpt_fr_%s' % PIDof
        PID=PIDsp
    
    #String containing the folder path
    jobdir = os.path.join('molecules', PID, MID, solname, spin, donacc, jobword)
    
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
        
        #Job not queued but was submitted
        if status == 'DNE':
            if jobword == 'opt':
                if os.path.exists(os.path.join(jobdir, 'end.xyz')):
                    status = 'done'
                else:
                    status = 'failed'
            if jobword == 'freq':
                output = _exec('grep Gibbs run.out', cwd=jobdir)
                line = output[0]
                if line.startswith('Gibbs'):
                    status = 'done'
                else:
                    status = 'failed'
            if job == 'energy':
                output = _exec('grep FINAL run.out', cwd = jobdir)
                line = output[0]
                if line.startswith('FINAL'):
                    status = 'done'
                else:
                    status = 'failed'
          
    #Nothing's been submitted yet
    else:
        optdir = os.path.join('molecules', PIDof, MID, solname, spin, donacc, 'opt')
        if jobword != 'opt':
            if not os.path.exists(os.path.join(optdir, 'end.xyz')):
                status = 'optinc'
        
    return status

                
def submit(MID,SID,PIDof,PIDsp,mold,lvld,sold,job,donacc,spin):
    """
    Creates necessary files then submits a TeraChem calculation job.
    
    Parameters
    ----------
    MID, SID, PID: str
                   molecule, solvent, and level of theory ID
    mold, lvld, sold: dict
                      Nested dictionaries created by parser
    job: str
         type of job to submit
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
        
    if job == "minimize":
        jobword = "opt"
        ofs = "na=$(head -1 scr/optim.xyz) && tail -$((na+2)) scr/optim.xyz > end.xyz"
        PID=PIDof
    if job == 'frequencies':
        jobword = "freq"
        ofs = ""
        PID=PIDof
    if job == 'energy':
        jobword = 'sglpt_fr_%s' % PIDof
        ofs = ""
        PID=PIDsp
    
    solname = sold[SID]['sname']
    dft = lvld[PID]['dft']
    if dft == "ub3lyp":
        dft = "b3lyp"
        unr = "unrestricted          true"
    if dft == "ubp86":
        dft = "bp86"
        unr = "unrestricted          true"
    else:
        unr = ""
    basis = lvld[PID]['basis']
    chg = mold[MID][key]
    dec = sold[SID]['epsilon']
    
    dirname = os.path.join('molecules', PID, MID, solname, spin, donacc, jobword)
    if not os.path.exists(dirname): os.makedirs(dirname)
    optdir = os.path.join('molecules', PIDof, MID, solname, spin, donacc, 'opt')
    if not os.path.exists(optdir): os.makedirs(optdir)
    
    #QChem
    if lvld[PID]['solvmod'] == 'smd':
      
        #String that will be written in to the run.in file
        qcinstr = """\
$molecule
{chg} {spin}
{coords}
$end

$rem
jobtype               sp
method                {dft}
basis                 mixed
symmetry              off
sym_ignore            true
{unr}
incdft                false
incfock               0
scf_convergence       8
thresh                14
max_scf_cycles        200
solvent_method        smd
$end

$smx
solvent               {solname}
$end

$basis
{coeffs}
$end
"""
            
        #String that will be written in to the batch.script file
        batchscriptstr = """\
#!/bin/bash -l
#SBATCH -p med
#SBATCH -N 1
#SBATCH -n {cpus}
#SBATCH -c 2
#SBATCH --gres=gpu:0
#SBATCH --mem={mem}
#SBATCH -J {MID}
#SBATCH -t 7-00:00:00

#SBATCH --no-requeue


# Record job info
echo -e "$SLURM_JOB_ID  $HOSTNAME  $(pwd)" >> ~/.myslurmhistory

export QC=/home/leeping/opt/qchem/5.1.1-openmp
export PATH=$QC/bin:$PATH
mkdir -p /scratch/cyeo99/$SLURM_JOB_ID
export QCLOCALSCR=/scratch/cyeo99/$SLURM_JOB_ID
export QCSCRATCH=.


qchem  -nt 4  qc.in qc.out

rm -r /scratch/cyeo99/$SLURM_JOB_ID
"""     
        #Write qc.in file
        if job == "minimize":
            startcoord = open(os.path.join('startxyz', "%s-%s.xyz" % (MID, name)), "r")
        else:
            startcoord = open(os.path.join(optdir, "end.xyz"), "r")
        coordlist = list(startcoord.readlines())
        atomnum = int(coordlist[0].replace("\n", ""))
        coords = coordlist[2:]
        coords = ''.join(coords)
        startcoord.close()
        
        #GET BASIS SET COEFFS
        if basis == "6-31gs_ldz":
            basisloc = os.path.join('molecules', "P01", MID, solname, spin, donacc, "opt", "scr", "start.basis")
        if basis == "tzvp_ltz":
            basisloc = os.path.join('molecules', "P02", MID, solname, spin, donacc, "opt", "scr", "start.basis")
        if basis == "6-31gs_ldz":
            basisloc = os.path.join('molecules', "P08", MID, solname, spin, donacc, "sglpt_fr_P02", "scr", "start.basis")
        basisfile = open(basisloc, "r")
        startbasis = list(basisfile.readlines())
        del startbasis[:3]
        for i in range(len(startbasis)):
            if startbasis[i].startswith("ATOM"):
                atomstr = startbasis[i].split()
                startbasis[i] = atomstr[-1] + "  " + atomstr[1].replace(":", "") + "\n"
            if startbasis[i].startswith("S"):
                startbasis[i] = startbasis[i].replace("\n", " 1.0\n")
            if startbasis[i].startswith("P"):
                startbasis[i] = startbasis[i].replace("\n", " 1.0\n")
            if startbasis[i].startswith("D"):
                startbasis[i] = startbasis[i].replace("\n", " 1.0\n")
            if startbasis[i].startswith("F"):
                startbasis[i] = startbasis[i].replace("\n", " 1.0\n")
            if startbasis[i] == "\n":
                startbasis[i] = "****\n"
        coeffs = ''.join(startbasis)
        basisfile.close()
        
        with open(os.path.join(dirname, 'qc.in'), 'w') as f:    
            f.write(qcinstr.format(unr=unr,chg=chg,spin=spin,coords=coords,dft=dft,solname=solname,coeffs=coeffs))

        #Depending on the number of atoms, set the number of CPUs and mem
        if atomnum < 30:
            cpus = 1
            mem = 8000
        if 30 <= atomnum < 50:
            cpus = 2
            mem = 16000
        if atomnum >= 50:
            cpus = 4
            mem = 32000
        #Write batch.script (doing this so that time limit can be set to 7 days, not 2 days (the default)
        with open(os.path.join(dirname, 'batch.script'), 'w') as g:
            g.write(batchscriptstr.format(mem=mem,MID=MID,cpus=cpus))
        
    #TeraChem
    else:
        #What to write if solvent model is PCM
        if lvld[PID]['solvmod'] == 'pcm':
            pcmorgas = "pcm cosmo\nepsilon %s" % dec
        #What to write if using "gas phase"
        if lvld[PID]['solvmod'] == 'gas':
            pcmorgas = ""
      
        #String that will be written in to the run.in file
        runinstr = """\
run {job}
method {dft}
dispersion yes
basis {basis}
charge {chg}
spinmult {spin}
coordinates start.xyz
{pcmorgas}
convthre 1e-6
threall  1e-18
dftgrid  3
thermotemp 298.15
scf diis+a

min_converge_e    1.0e-7
min_converge_grms 3.0e-5
min_converge_gmax 4.5e-5
min_converge_drms 1.2e-4
min_converge_dmax 1.8e-4
"""
            
        #String that will be written in to the batch.script file
        batchscriptstr = """\
#!/bin/bash -l
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH -n {gpus}
#SBATCH -c 2
#SBATCH --gres=gpu:{gpus}
#SBATCH --mem={mem}
#SBATCH -J {MID}
#SBATCH -t 7-00:00:00

#SBATCH --no-requeue

# Record job info
echo -e "$SLURM_JOB_ID  $HOSTNAME  $(pwd)" >> ~/.myslurmhistory

module load intel cuda/9.0
export CUDA_CACHE_PATH=/scratch/$USER/.nv/ComputeCache
export TeraChem=/home/leeping/opt/terachem/current
export PATH=$TeraChem/bin:$PATH
export LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH
    
    
terachem  run.in &> run.out

{ofs}
"""     
        #Write run.in file
        with open(os.path.join(dirname, 'run.in'), 'w') as f:    
            f.write(runinstr.format(job=job, dft=dft, basis=basis, chg=chg, spin=spin, pcmorgas=pcmorgas))
            
        if job == "minimize":
            #Copy over starting structure and rename as start.xyz     
            shutil.copy(os.path.join('startxyz', "%s-%s.xyz" % (MID, name)), dirname)
            os.rename(os.path.join(optdir, "%s-%s.xyz" % (MID, name)), os.path.join(dirname, "start.xyz"))
        else:
            #Copy over last frame of optim.xyz from opt, rename as start.xyz
            M = Molecule(os.path.join(optdir, 'scr', 'optim.xyz'))
            M[-1].write(os.path.join(dirname, 'start.xyz'))
    
        #Read start.xyz file to check for number of atoms
        startcoord = open(os.path.join(dirname, "start.xyz"), "r")
        #Turn the first word of the first line into an integer
        atomnum = int(list(startcoord.readlines())[0].split()[0])
        startcoord.close()
        #Depending on the number of atoms, set the number of GPUs to either 1, 2, or 4.
        if atomnum < 30:
            gpus = 1
            mem = 8000
        if 30 <= atomnum < 50:
            gpus = 2
            mem = 16000
        if atomnum >= 50:
            gpus = 4
            mem = 32000
        #Write batch.script (doing this so that time limit can be set to 7 days, not 2 days (the default)
        with open(os.path.join(dirname, 'batch.script'), 'w') as g:
            g.write(batchscriptstr.format(gpus=gpus, mem=mem, MID=MID, ofs=ofs))
      
    #Submit job
    submit = _exec('sbatch batch.script', cwd=dirname)
    submitstring = submit[0]
    submitsplt = submitstring.split()
    jobnum = submitsplt[3]
    #Save job number to txt file
    with open(os.path.join(dirname, 'submittedjob.txt'), 'w') as h:
        h.write(jobnum)


def gethyd(MID,SID,DID,mold,sold,modd,lvld,donacc):
    #---get sgl pt seperately---
    """
    Reads run.out file to get the Gibbs free energy of the optimal structure for a molecule.
    
    Parameters
    ----------
    MID, SID, DID: str
                   molecule, solvent, and model ID
    mold, sold, modd: dict
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
        key = "donspn"
    else:
        key = "accspn"
    solname = sold[SID]['sname']
    Edict= {}

    if len(modd[DID])==1:
        PID = modd[DID]['optfreq']
        #Loop over all spin multiplicities for the molecule of interest
        for sm in mold[MID][key]:
            freqdir = os.path.join('molecules', PID, MID, solname, sm, donacc, 'freq')
            energy = _exec('grep Gibbs run.out', cwd=freqdir)
            Gibbs = float(energy[0].split()[-2])
            #Store to dictionary
            Edict[sm] = Gibbs
    else:
        PIDof = modd[DID]['optfreq']
        PIDsp = modd[DID]['sglpt']
        sglptname = 'sglpt_fr_%s' % PIDof
        for sm in mold[MID][key]:
            freqdir = os.path.join('molecules', PIDof, MID, solname, sm, donacc, 'freq')
            sglptdir = os.path.join('molecules', PIDsp, MID, solname, sm, donacc, sglptname)
            freqout = _exec("grep 'Free Energy Correction' run.out", cwd=freqdir)
            if lvld[PIDsp]['solvmod'] == 'smd':
                sglptout = _exec("grep '(6)  G-S(liq) free energy of system' qc.out", cwd=sglptdir)
            else:
                sglptout = _exec('grep FINAL run.out', cwd=sglptdir)
            freeEcorr = float(freqout[0].split()[-2])
            finalE = float(sglptout[0].split()[-2]) * 627.509 #convert from AU to kcal/mol
            Gibbs = freeEcorr + finalE
            Edict[sm] = Gibbs
            
    #Search dictionary for spin mult with lowest energy
    minsm = min(Edict, key=Edict.get)
    minG = Edict[minsm]
    
    return minsm, minG

def bigbraintime(MID,SID,PIDof,PIDsp,mold,lvld,sold,jobtype,donacc,spin):
    """
    Decides what to do based on job status
    """
    if jobtype == "minimize":
        name = "Geometry optimization"
        nextstep = "submit frequency or single point jobs."
    if jobtype == "frequencies":
        name = "Frequency analysis"
        nextstep = "calculate hydricity."
    if jobtype == "energy":
        name = "Single point calculation"
        nextstep = "calculate hydricity."
        
    stat = checkjobstatus(MID,SID,PIDof,PIDsp,sold,jobtype,donacc,spin)
    
    if stat == "DNE":
        submit(MID,SID,PIDof,PIDsp,mold,lvld,sold,jobtype,donacc,spin)
    if stat == "queued":
        print("%s for %s %s already queued." % (name, MID, donacc))
    if stat == "failed":
        print("%s for %s %s failed." % (name, MID, donacc))
    if stat == "done":
        print("%s done for %s %s, time to %s" % (name, MID, donacc, nextstep))
    if stat == "optinc":
        print("Geometry optimization for %s %s incomplete." % (MID, donacc))
        
def dataanalysis(xlist, ylist, fixedslope):
    #Add trendline and calculate R^2 & RMSE
    #Get trendline with slope fixed
    yintlist = []
    for i in range(len(xlist)):
        yintlist.append(ylist[i] - fixedslope * xlist[i])
    yint = np.mean(yintlist)
    
    #Get R^2
    #residual
    predictedy = fixedslope*xlist+yint
    residual = []
    for i in range(len(ylist)):
        residual.append(ylist[i] - predictedy[i]) #actual y - predicted y
    ressq = []
    for i in range(len(residual)):
        ressq.append(residual[i]**2)
    #average y
    aveexp = np.sum(ylist)/len(ylist)
    diffmean = []
    for i in range(len(ylist)):
        diffmean.append(ylist[i] - aveexp) #actual y - average y
    diffmeansq = []
    for i in range(len(diffmean)):
        diffmeansq.append(diffmean[i]**2)
    #R^2
    ressqsum = np.sum(ressq)
    diffmeansqsum = np.sum(diffmeansq)
    Rsquared = 1 - ressqsum/diffmeansqsum
    
    #Get RMSE
    difflist=[]
    for i in range(len(xlist)):
        diffexpcalsq = (predictedy[i]-ylist[i])**2
        difflist.append(diffexpcalsq)
    diffsqsum = np.sum(difflist)
    RMSE = np.sqrt(diffsqsum/len(xlist))
    
    return Rsquared, RMSE, yint

def savedata(project, DID, sysd, mold, sold, modd,lvld):
    """
    Saves free energy of hydricity half reaction to dat file, vertically shifts data points according to solvent,
    which gives hydricity, and makes plots for each solvent and one plot that includes points for all solvents.
    """
    if project == 'organic':
        doorg = True
        doorm = False
    if project == 'organometallic':
        doorg = False
        doorm = True
    if project == 'both':
        doorg = True
        doorm = True
    #Lists to store calculated delta G_HHRs and experimental hydricities, according to solvent and molecule type
    HHRdict = {}
    #Lists to store vertically shifted calculated values
    hyddict = {}
        
    for solvent in sold:
        solname = sold[solvent]['sname']
        
        
        calOR = "cal_%s_OR" % solname
        expOR = "exp_%s_OR" % solname
        calOM = "cal_%s_OM" % solname
        expOM = "exp_%s_OM" % solname
        
        if doorg == True:
            HHRdict[calOR] = []
            
            hyddict[calOR] = []
            hyddict[expOR] = []
        if doorm == True:
            HHRdict[calOM] = []
            
            hyddict[calOM] = []
            hyddict[expOM] = []
        
    #Make directories for data
    datadir = os.path.join('data', '%s' %DID)
    if not os.path.exists(datadir): os.makedirs(datadir)

    #Add column titles
    with open(os.path.join(datadir, 'delG_HHR-%s.dat' % DID), 'w') as f:
        f.write("YID       don chg   don spin  acc chg   acc spin  Calculated Hydricity     Experimental Hydricity\n")
        f.write("---       -------   --------  -------   --------  --------------------     ----------------------\n")

    for YID in sysd:
        try:
            #Calculate delta G_HHR
            MID = sysd[YID]['molecule']
            SID = sysd[YID]['solvent']
            donminspn, donfreeE = gethyd(MID,SID,DID,mold,sold,modd,lvld,'donor')
            accminspn, accfreeE = gethyd(MID,SID,DID,mold,sold,modd,lvld,'acceptor')
            delG_HHR = accfreeE - donfreeE
                
            #Save results to dat file (YID, donor charge, donor spin, acceptor charge, acceptor spin, free energy of hydricity half reaction)
            with open(os.path.join(datadir, 'delG_HHR-%s.dat' % DID), 'a+') as f:
                f.write('{0:10}{1:10}{2:10}{3:10}{4:10}{5:25}{6:25}\n'.format(YID, mold[MID]['donchg'], donminspn, mold[MID]['accchg'], accminspn, str(delG_HHR), sysd[YID]['hyd']))
            
            #Store points to lists accordingly        
            if doorg == True:
                if mold[MID]['type'] == 'OR':
                    calOR = "cal_%s_OR" % sold[SID]['sname']
                    expOR = "exp_%s_OR" % sold[SID]['sname']
                    HHRdict[calOR].append(delG_HHR)
                    hyddict[expOR].append(float(sysd[YID]['hyd']))
            if doorm == True:
                if mold[MID]['type'] == 'OM':
                    calOM = "cal_%s_OM" % sold[SID]['sname']
                    expOM = "exp_%s_OM" % sold[SID]['sname']
                    HHRdict[calOM].append(delG_HHR)
                    hyddict[expOM].append(float(sysd[YID]['hyd']))
        except:
            print("Calculation failed, check %s" %MID)


    #Make directories for plots
    plotdir = os.path.join('plots', '%s' %DID)
    if not os.path.exists(plotdir): os.makedirs(plotdir)
    
    #For all molecules: This plot is PURELY for visual purposes
    figall, axall = plt.subplots()

    #For each solvent
    for solvent in sold:
        solname = sold[solvent]['sname']
        
        calOR = "cal_%s_OR" % solname
        expOR = "exp_%s_OR" % solname
        calOM = "cal_%s_OM" % solname
        expOM = "exp_%s_OM" % solname
        
        #Make scatter plot
        figsolv, axsolv = plt.subplots()

        #Combine HHR data to single list
        if project == 'both':
            solvcal = HHRdict[calOR] + HHRdict[calOM]
            solvexp = hyddict[expOR] + hyddict[expOM]
        
        #Calculate vertical shift that corresponds to free energy of hydride anion in solvent
        #and get R^2 and RMSE for all, organic, and organometallic molecules
        fixedslope = 1
        
        if project == 'both':
            rsqall, RMSEall, yintall = dataanalysis(solvexp, solvcal, fixedslope)
            rsqall = str(round(rsqall, 3))
            RMSEall = str(round(RMSEall, 3))
        
        if doorg == True:
            if len(hyddict[expOR])>1:
                rsqorg, RMSEorg, yintorg = dataanalysis(hyddict[expOR], HHRdict[calOR], fixedslope)
                rsqorg = str(round(rsqorg, 3))
                RMSEorg = str(round(RMSEorg, 3))
            else:
                rsqorg = "N/A"
                RMSEorg = "N/A"
                yintorg = "N/A"
        
        if doorm == True:
            if len(hyddict[expOM])>1:
                rsqorm, RMSEorm, yintorm = dataanalysis(hyddict[expOM], HHRdict[calOM], fixedslope)
                rsqorm = str(round(rsqorm, 3))
                RMSEorm = str(round(RMSEorm, 3))
            else:
                rsqorm = "N/A"
                RMSEorm = "N/A"
                yintorm = "N/A"
            
        if project == 'both':
            #Store y-intercepts to dat file.
            with open(os.path.join(datadir, 'yint-%s-%s.dat' % (DID, solname)), 'w') as g:
                g.write('{:20}{:20}{:20}\n{:20}{:20}{:20}'.format("All", "Organic", "Organometallic",\
                        str(yintall), str(yintorg), str(yintorm)))
            
            #Apply vertical shift and save to hydricity dictionary
            for i in range(len(HHRdict[calOR])):
                hyddict[calOR].append(HHRdict[calOR][i] - yintall)
                
            for i in range(len(HHRdict[calOM])):
                hyddict[calOM].append(HHRdict[calOM][i] - yintall)
        if project == 'organic':
            #Store y-intercepts to dat file.
            with open(os.path.join(datadir, 'yint-%s-%s.dat' % (DID, solname)), 'w') as g:
                g.write(str(yintorg))
            
            #Apply vertical shift and save to hydricity dictionary
            if type(yintorg)==str:
                for i in range(len(HHRdict[calOR])):
                    hyddict[calOR].append(HHRdict[calOR][i])
            else:
                for i in range(len(HHRdict[calOR])):
                    hyddict[calOR].append(HHRdict[calOR][i] - yintorg)
        if project == 'organometallic':
            #Store y-intercepts to dat file.
            with open(os.path.join(datadir, 'yint-%s-%s.dat' % (DID, solname)), 'w') as g:
                g.write(str(yintorm))
            
            #Apply vertical shift and save to hydricity dictionary
            if type(yintorm)==str:
                for i in range(len(HHRdict[calOM])):
                    hyddict[calOM].append(HHRdict[calOM][i])
            else:
                for i in range(len(HHRdict[calOM])):
                    hyddict[calOM].append(HHRdict[calOM][i] - yintorm)

        #Add trendline
        x = np.linspace(0, 140, 1000)
        axsolv.plot(x, fixedslope*x, linestyle = 'dashed', color = 'black')
        
        #Plot points
        if doorg == True:
            organic = axsolv.scatter(hyddict[expOR], hyddict[calOR], c='turquoise')
        if doorm == True:
            organometallic = axsolv.scatter(hyddict[expOM], hyddict[calOM], c='darkgray')

        #Add labels and legend
        axsolv.set_title('Calculated vs Experimental Hydricities in %s' % solname.capitalize())
        axsolv.set_xlabel('Experimental Hydricity (kcal/mol)')
        axsolv.set_ylabel('Calculated Hydricity (kcal/mol)')
        axsolv.set_xlim([0,140])
        axsolv.set_ylim([0,140])
        axsolv.set_aspect('equal', adjustable='box')
        extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        
        color = next(axall._get_lines.prop_cycler)['color']
        
        if project == 'both':
            RsqRMSE='{:4}{:8}{:8}\n{:4}{:8}{:8}\n{:4}{:8}{:8}\n{:4}{:8}{:8}'.format(' ', 'R^2', 'RMSE',\
                     'All', rsqall, RMSEall, 'OR', rsqorg, RMSEorg, 'OM', rsqorm, RMSEorm)
            axsolv.legend([organic, organometallic, extra], ('Organic (OR)', 'Organometallic (OM)', RsqRMSE), prop={'family': 'monospace'})
            
            labelOR = solname + ', organic'
            labelOM = solname + ', organometallic'
            
            axall.scatter(hyddict[expOR], hyddict[calOR], color=color, marker='x', label=labelOR)
            axall.scatter(hyddict[expOM], hyddict[calOM], color=color, marker='.', label=labelOM)
        if project == 'organic':
            RsqRMSE='{:8}{:8}\n{:8}{:8}'.format('R^2', 'RMSE', rsqorg, RMSEorg)
            axsolv.legend([organic, extra], ('Organic', RsqRMSE), prop={'family': 'monospace'})
            
            axall.scatter(hyddict[expOR], hyddict[calOR], color=color, label=solname)
        if project == 'organometallic':
            RsqRMSE='{:8}{:8}\n{:8}{:8}'.format('R^2', 'RMSE', rsqorm, RMSEorm)
            axsolv.legend([organometallic, extra], ('Organometallic', RsqRMSE), prop={'family': 'monospace'})
            
            axall.scatter(hyddict[expOM], hyddict[calOM], color=color, label=solname)
        
        #Save plot to pdf
        figsolv.savefig(os.path.join(plotdir, '%s_%s.pdf' % (DID, solname)))
    
    #Add labels and legend
    axall.set_title('All Solvents')
    axall.set_xlabel('Experimental hydricity (kcal/mol)')
    axall.set_ylabel('Calculated Hydricity (kcal/mol)')
    axall.set_aspect('equal', adjustable='box')
    axall.set_xlim([0,140])
    axall.set_ylim([0,140])
    axall.legend(fontsize=9.8)
    
    #Save plot to pdf
    figall.savefig(os.path.join(plotdir, '%s_all.pdf' % DID))
    

def main():
    #Ask user to provide level of theory, job type, and name of file containing molecule info
    project = str(input('Enter molecule class: \n'))
    DID = str(input('Enter model ID: \n'))
    jobtype = str(input('Enter job type: \n'))
    molfile = os.path.join("keys", "molecules.txt")
    
    #Create dictionaries
    mold, lvld, sold, sysd, modd = parser(project, molfile)
    
    #If freq analysis done for all spin mults, calculate hydricity, then store data to dat file and plots
    if jobtype == "hydricity":
        savedata(project, DID, sysd, mold, sold, modd, lvld)
    
    else:
        PIDof = modd[DID]['optfreq']
        if len(modd[DID])>1:
            PIDsp = modd[DID]['sglpt']
        else:
            PIDsp = PIDof
            
        for YID in sysd:
            MID = sysd[YID]['molecule']
            SID = sysd[YID]['solvent']
            #Using spin mults as key since number of files to be created depends on number of possible spin mults
            for key in mold[MID]['donspn']:
                bigbraintime(MID,SID,PIDof,PIDsp,mold,lvld,sold,jobtype,'donor',key)
            #Same thing for acceptor structures
            for key in mold[MID]['accspn']:
                bigbraintime(MID,SID,PIDof,PIDsp,mold,lvld,sold,jobtype,'acceptor',key)

if __name__ == "__main__":
    main()
