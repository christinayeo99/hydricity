#!/usr/bin/env python

import os, sys
import numpy as np
import IPython
import shutil
from forcebalance.nifty import _exec
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def parser(molfile):
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

def checkjobstatus(MID,SID,PID,sold,job,donacc,spin):
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
    
    if job == 'minimize':
        optfreq = 'opt'
    if job == 'frequencies':
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
        
    if job == "minimize":
        optfreq = "opt"
        optorfreq = "na=$(head -1 scr/optim.xyz) && tail -$((na+2)) scr/optim.xyz > end.xyz"
    else:
        optfreq = "freq"
        optorfreq = ""
    
    solname = sold[SID]['sname']
    dft = lvld[PID]['dft']
    basis = lvld[PID]['basis']
    chg = mold[MID][key]
    dec = sold[SID]['epsilon']
    
    dirname = os.path.join('molecules', PID, MID, solname, spin, donacc, optfreq)
    if not os.path.exists(dirname): os.makedirs(dirname)
    optdir = os.path.join('molecules', PID, MID, solname, spin, donacc, 'opt')
    if not os.path.exists(optdir): os.makedirs(optdir)
    
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

{optorfreq}
"""     
    #Write run.in file
    with open(os.path.join(dirname, 'run.in'), 'w') as f:    
        f.write(runinstr.format(job=job, dft=dft, basis=basis, chg=chg, spin=spin, pcmorgas=pcmorgas))
        
    if job == "minimize":
        #Copy over starting structure and rename as start.xyz     
        shutil.copy(os.path.join('startxyz', "%s-%s.xyz" % (MID, name)), dirname)
        os.rename(os.path.join(optdir, "%s-%s.xyz" % (MID, name)), os.path.join(dirname, "start.xyz"))
    else:
        #Copy over end.xyz file from opt, rename as start.xyz
        shutil.copy(os.path.join(optdir, 'end.xyz'), dirname)
        os.rename(os.path.join(dirname, 'end.xyz'), os.path.join(dirname, 'start.xyz'))

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
        g.write(batchscriptstr.format(gpus=gpus, mem=mem, MID=MID, optorfreq=optorfreq))
        
    #Submit job
    submit = _exec('sbatch batch.script', cwd=dirname)
    submitstring = submit[0]
    submitsplt = submitstring.split()
    jobnum = submitsplt[3]
    #Save job number to txt file
    with open(os.path.join(dirname, 'submittedjob.txt'), 'w') as h:
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
        key = "donspn"
    else:
        key = "accspn"
    solname = sold[SID]['sname']
    Edict= {}

    #Loop over all spin multiplicities for the molecule of interest
    for sm in mold[MID][key]:
        freqdir = os.path.join('molecules', PID, MID, solname, sm, donacc, 'freq')
        energy = _exec('grep Gibbs run.out', cwd=freqdir)
        Gibbs = float(energy[0].split()[-2])
        #Store to dictionary
        Edict[sm] = Gibbs
    #Search dictionary for spin mult with lowest energy
    minsm = min(Edict, key=Edict.get)
    minG = Edict[minsm]
    
    return minsm, minG

def bigbraintime(MID,SID,PID,mold,lvld,sold,jobtype,donacc,spin):
    """
    Decides what to do based on job status
    """
    if jobtype == "minimize":
        name = "Geometry optimization"
        nextstep = "submit frequency jobs."
    if jobtype == "frequencies":
        name = "Frequency analysis"
        nextstep = "calculate hydricity."
        
    stat = checkjobstatus(MID,SID,PID,sold,jobtype,donacc,spin)
    
    if stat == "DNE":
        submit(MID,SID,PID,mold,lvld,sold,jobtype,donacc,spin)
    if stat == "queued":
        print("%s for %s %s already queued." % (name, MID, donacc))
    if stat == "failed":
        print("%s for %s %s failed." % (name, MID, donacc))
    if stat == "done":
        print("%s done for %s %s, time to %s" % (name, MID, donacc, nextstep))
    if stat == "optinc":
        print("Geometry optimization for %s %s incomplete." % (MID, donacc))
        
def dataanalysis(xlist, ylist, fixedslope):
    #Add trendline and calculate R^2 & RMSD
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
        residual.append(ylist[i] - predictedy[i])
    ressq = []
    for i in range(len(residual)):
        ressq.append(residual[i]**2)
    #average y
    aveexp = np.sum(ylist)/len(ylist)
    diffmean = []
    for i in range(len(ylist)):
        diffmean.append(ylist[i] - aveexp)
    diffmeansq = []
    for i in range(len(diffmean)):
        diffmeansq.append(diffmean[i]**2)
    #R^2
    ressqsum = np.sum(ressq)
    diffmeansqsum = np.sum(diffmeansq)
    Rsquared = 1 - ressqsum/diffmeansqsum
    
    #Get RMSD
    difflist=[]
    for i in range(len(xlist)):
        diffexpcalsq = (xlist[i]-ylist[i])**2
        difflist.append(diffexpcalsq)
    diffsqsum = np.sum(difflist)
    RMSD = np.sqrt(diffsqsum/len(xlist))
    
    return Rsquared, RMSD, yint

def savedata(PID, sysd, mold, sold):
    """
    Saves data to dat file, makes and saves plots as pdfs.
    """
    #Lists to store calculated and experimental hydricities, according to solvent and molecule type
    listdict = {}
    for solvent in sold:
        calOR = "cal_%s_OR" % sold[solvent]['sname']
        expOR = "exp_%s_OR" % sold[solvent]['sname']
        calOM = "cal_%s_OM" % sold[solvent]['sname']
        expOM = "exp_%s_OM" % sold[solvent]['sname']
        
        listdict[calOR] = []
        listdict[expOR] = []
        listdict[calOM] = []
        listdict[expOM] = []
    
    if not os.path.exists('data'): os.makedirs('data')
    with open(os.path.join('data', 'data-%s.dat' % PID), 'w') as f:
        f.write("YID       don chg   don spin  acc chg   acc spin  Calculated Hydricity     Experimental Hydricity\n")
        f.write("---       -------   --------  -------   --------  --------------------     ----------------------\n")
    for YID in sysd:
        #Calculate hydricity
        MID = sysd[YID]['molecule']
        SID = sysd[YID]['solvent']
        donminspn, donfreeE = gethyd(MID,SID,PID,mold,sold,'donor')
        accminspn, accfreeE = gethyd(MID,SID,PID,mold,sold,'acceptor')
        delG_HHR = accfreeE - donfreeE
        if sold[SID]['sname'] == 'water':
            hydricity = delG_HHR - 420.7
        else: #acetonitrile, benzonitrile, DMSO
            hydricity = delG_HHR - 419.34
            
        #Save results to dat file (YID, donor charge, donor spin, acceptor charge, acceptor spin, hydricity)
        with open(os.path.join('data', 'data-%s.dat' % PID), 'a+') as f:
            f.write('{0:10}{1:10}{2:10}{3:10}{4:10}{5:25}{6:25}\n'.format(YID, mold[MID]['donchg'], donminspn, mold[MID]['accchg'], accminspn, str(hydricity), sysd[YID]['hyd']))
        
        #Store points to lists accordingly
        if mold[MID]['type'] == 'OR':
            calOR = "cal_%s_OR" % sold[SID]['sname']
            expOR = "exp_%s_OR" % sold[SID]['sname']
            listdict[calOR].append(hydricity)
            listdict[expOR].append(float(sysd[YID]['hyd']))
        else:
            calOM = "cal_%s_OM" % sold[SID]['sname']
            expOM = "exp_%s_OM" % sold[SID]['sname']
            listdict[calOM].append(hydricity)
            listdict[expOM].append(float(sysd[YID]['hyd']))


    #Make Plots
    plotdir = os.path.join('plots', '%s' %PID)
    if not os.path.exists(plotdir): os.makedirs(plotdir)
    
    #All
    #Make scatter plot
    figall, axall = plt.subplots()

    for solvent in sold:
        solname = sold[solvent]['sname']
        
        calOR = "cal_%s_OR" % solname
        expOR = "exp_%s_OR" % solname
        calOM = "cal_%s_OM" % solname
        expOM = "exp_%s_OM" % solname
        
        color = next(axall._get_lines.prop_cycler)['color']
        
        labelOR = solname + ', organic'
        labelOM = solname + ', organometallic'
        
        axall.scatter(listdict[expOR], listdict[calOR], color=color, marker='x', label=labelOR)
        axall.scatter(listdict[expOM], listdict[calOM], color=color, marker='.', label=labelOM)
    
    #Add labels and legend
    axall.set_title('All Solvents')
    axall.set_xlabel('Experimental hydricity (kcal/mol)')
    axall.set_ylabel('Calculated Hydricity (kcal/mol)')
    axall.set_aspect('equal', adjustable='box')
    axall.set_xlim([0,140])
    axall.set_ylim([0,140])
    axall.legend(fontsize=9.8, prop={'family': 'monospace'})
    
    #Save plot to pdf
    figall.savefig(os.path.join(plotdir, '%s_all.pdf' % PID))
    
    
    #For each solvent
    for solvent in sold:
        solname = sold[solvent]['sname']
        calOR = "cal_%s_OR" % solname
        expOR = "exp_%s_OR" % solname
        calOM = "cal_%s_OM" % solname
        expOM = "exp_%s_OM" % solname
        
        #Make scatter plot
        figsolv, axsolv = plt.subplots()
        organic = axsolv.scatter(listdict[expOR], listdict[calOR], c='turquoise')
        organometallic = axsolv.scatter(listdict[expOM], listdict[calOM], c='darkgray')
        
        #Combine data to single list
        solvcal = listdict[calOR] + listdict[calOM]
        solvexp = listdict[expOR] + listdict[expOM]
        
        #Get R^2 and RMSD for all, organic, and organometallic molecules, then change type from float to str
        fixedslope = 1
        
        rsqall, rmsdall, yintall = dataanalysis(solvexp, solvcal, fixedslope)
        rsqall = str(round(rsqall, 3))
        rmsdall = str(round(rmsdall, 3))
        
        if len(listdict[expOR])>1:
            rsqorg, rmsdorg, yintorg = dataanalysis(listdict[expOR], listdict[calOR], fixedslope)
            rsqorg = str(round(rsqorg, 3))
            rmsdorg = str(round(rmsdorg, 3))
        else:
            rsqorg = "N/A"
            rmsdorg = "N/A"
        if len(listdict[expOM])>1:
            rsqorm, rmsdorm, yintorm = dataanalysis(listdict[expOM], listdict[calOM], fixedslope)
            rsqorm = str(round(rsqorm, 3))
            rmsdorm = str(round(rmsdorm, 3))
        else:
            rsqorm = "N/A"
            rmsdorm = "N/A"
        
        #Add trendline
        x = np.linspace(0, 140, 1000)
        axsolv.plot(x, fixedslope*x+yintall, linestyle = 'dashed', color = 'black')

        #Add labels and legend
        axsolv.set_title('%s' % solname.capitalize())
        axsolv.set_xlabel('Experimental hydricity (kcal/mol)')
        axsolv.set_ylabel('Calculated Hydricity (kcal/mol)')
        axsolv.set_xlim([0,140])
        axsolv.set_ylim([0,140])
        axsolv.set_aspect('equal', adjustable='box')
        extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        
        RsqRMSD='{:4}{:8}{:8}\n{:4}{:8}{:8}\n{:4}{:8}{:8}\n{:4}{:8}{:8}'.format(' ', 'R^2', 'RMSD',\
                 'All', rsqall, rmsdall, 'OR', rsqorg, rmsdorg, 'OM', rsqorm, rmsdorm)
        
        axsolv.legend([organic, organometallic, extra], ('Organic (OR)', 'Organometallic (OM)',RsqRMSD), prop={'family': 'monospace'})
        
        #Save plot to pdf
        figsolv.savefig(os.path.join(plotdir, '%s_%s.pdf' % (PID, solname)))
    

def main():
    #Ask user to provide level of theory, job type, and name of file containing molecule info
    PID = str(input('Enter level of theory ID: \n'))
    jobtype = str(input('Enter job type: \n'))
    molfile = str(input('Enter name of file containing molecule info: \n'))
    
    #Create dictionaries
    mold, lvld, sold, sysd = parser(molfile)
    
    #If freq analysis done for all spin mults, calculate hydricity, then store data to dat file and plots
    if jobtype == "hydricity":
        savedata(PID, sysd, mold, sold)
    
    else:
        for YID in sysd:
            MID = sysd[YID]['molecule']
            SID = sysd[YID]['solvent']
            #Using spin mults as key since number of files to be created depends on number of possible spin mults
            for key in mold[MID]['donspn']:
                bigbraintime(MID,SID,PID,mold,lvld,sold,jobtype,'donor',key)
            #Same thing for acceptor structures
            for key in mold[MID]['accspn']:
                bigbraintime(MID,SID,PID,mold,lvld,sold,jobtype,'acceptor',key)

if __name__ == "__main__":
    main()
