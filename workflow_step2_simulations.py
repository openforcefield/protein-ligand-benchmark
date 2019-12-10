#!/usr/bin/env python
# coding: utf-8

# In[21]:


import sys
import numpy as np
import os, shutil
import re
import subprocess
import glob
import argparse
from PLBenchmarks import targets, edges


# # Overview
# The workflow is divided into 3 steps
# 
# `Step 1`: preparation of hybrid structures and topologies.
# Here, ligand maps from the literature are used: this ensures equivalent comparison to the earlier calculations.
# Atom mapping between ligand pairs is generated with the script `atoms_to_morph.py` script.
# Hybrid structure and topology is generated with `make_hybrid.py` script.
# Further standard steps: placing system in the box, solvating, adding ions, creating energy minimization input.
# This steps concludes with the energy minimization step (best to be run on the cluster).
# 
# `Step 2`: equilibrium runs.
# 
# `Step 3`: non-equilibrium transitions.

# In[22]:


# here are some functions that are used in this workflow
def printInfo(runtype='em', run=1, target='xxx', edge='yyy_zzz', wp='water', state='stateA'):
    print('=' * 80)
    print(f'=== {runtype}{run:<5d} {target:8s} {edge:29s} {wp:12s} {state:12s} ===' )
    print('=' * 80)

# read edges from a file
def read_edges( target ):
    df = edges.getEdgesSet(target)
    ed = {}
    for i, e in df.iterrows():
        print(e[0])
        ed[f'edge_{e[0]}_{e[1]}'] = [f'lig_{e[0]}', f'lig_{e[1]}']
    return(ed)


# create a folder
def create_folder( path, fname ):
    if not os.path.exists(path+'/'+fname):
        os.makedirs(path+'/'+fname)


####################################
# functions to generate jobscripts #
####################################

def decide_on_resources( wp, simType ):
    simtime = 1
    simcpu = 4
    if simType=='eq':
        simtime = 2
        simcpu = 8
        if wp=='protein':
            simtime = 7
            simcpu = 8
    return(simtime,simcpu)


def create_SGE_jobscript( fname, simpath, jobname, runtype, run, simtime=4, simcpu=1 ):
    fp = open(fname,'w')
    fullsimpath = os.path.abspath(simpath)

    jobline = f'#! /usr/bin/bash\n\
#$ -N {jobname}\n\
#$ -cwd\n
#$ -q all.q\n
#$ -j yes\n
#$ -l slot_type=gromacs,affinity_group=default\n\n
source /home/dhahn3/bin/gromacs-patched/bin/GMXRC\n\n\

gmx mdrun -ntmpi 1 -ntomp ${simcpu} -pme gpu -pin on -s {runtype}{run}.tpr -deffnm {runtype}{run}\n\n'
    fp.write(jobline)    
    fp.close()
    
    
    
def create_SLURM_jobscript( fname, simpath, jobname, runtype, run, simtime=4, simcpu=1 ):
    fp = open(fname,'w')
    fullsimpath = os.path.abspath(simpath)

    jobline = f'#! /usr/bin/bash\n\
#SBATCH --job-name={jobname}\n\
#SBATCH --error={jobname}.e\n\
#SBATCH --partition=nes2.8\n\n\
#SBATCH --nodes=1\n\
#SBATCH --constraint=wes2.8\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --cpus-per-task={simcpu}\n\
#SBATCH --time={simtime}-00:00:00\n\
\n\
#Function to call to run the actual code\n\
slurm_startjob(){{\n\
#----------------- Actual calculation command goes here: ---------------------------\n\
\n\
module purge\n\
module load gnu/8.2.0  \n\
module load openmpi/3.1.2\n\
\n\
source /export/home/dfhahn/bin/gromacs_2019.4/bin/GMXRC\n\
PATH=/export/home/dfhahn/bin/gromacs_2019.4/bin/:$PATH\n\
\n\
WORKDIR=$(pwd)\n\
mkdir -p $TMPDIR/$SLURM_JOB_NAME\n\
cd $TMPDIR/$SLURM_JOB_NAME\n\
\n\
cp $WORKDIR/{runtype}{run}.tpr .\n\
gmx mdrun -ntomp 8 -ntmpi 1 -s {runtype}{run}.tpr -deffnm {runtype}{run}\n\
rsync -avzui ./ $WORKDIR/\n\
}}\n\
\n\
slurm_info_out(){{\n\
\n\
echo "=================================== SLURM JOB ==================================="\n\
date\n\
echo\n\
echo "The job will be started on the following node(s):"\n\
echo $SLURM_JOB_NODELIST\n\
echo\n\
echo "Slurm User:         $SLURM_JOB_USER"\n\
echo "Run Directory:      $(pwd)"\n\
echo "Job ID:             ${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"\n\
echo "Job Name:           $SLURM_JOB_NAME"\n\
echo "Partition:          $SLURM_JOB_PARTITION"\n\
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"\n\
echo "Number of tasks:    $SLURM_NTASKS"\n\
echo "Submitted From:     $SLURM_SUBMIT_HOST:$SLURM_SUBMIT_DIR"\n\
echo "=================================== SLURM JOB ==================================="\n\
echo\n\
echo "--- SLURM job-script output ---"\n\
                    }}\n\
\n\
slurm_info_out\n\
\n\
slurm_startjob\n'

    
    fp.write(jobline)    
    fp.close()


def prepareSimulations(target, forcefield, runtype, queueType='sge', verbose=False):
    # all simulations will be done in 3 replicas
    replicas = 3

    for t in targets.target_list:
        if t['name'] == target:
            targetDir = t['dir']
            break
    else:
        print('Target not found')
        return 0

    # set paths relative to base
    inputpath = './input'
    mdppath = './PLBenchmarks/workflow/mdp'

    # path where to find coordinates, topologies, hybrid topologies, simulations...
    ligPath = f'PLBenchmarks/data/{targetDir}/03_docked/'
    topPath = f'PLBenchmarks/data/{targetDir}/04_topo/{forcefield}/'
    hybPath = f'PLBenchmarks/data/{targetDir}/05_hybrid/{forcefield}'

    runpath =  {'em' : f'PLBenchmarks/data/{targetDir}/06_em/{forcefield}',
                'nvt': f'PLBenchmarks/data/{targetDir}/07_nvt/{forcefield}',
                'eq' : f'PLBenchmarks/data/{targetDir}/08_eq/{forcefield}',
                'morphes': f'PLBenchmarks/data/{targetDir}/09_morphes/{forcefield}'
            }
    
    def getRunCoord(runtype, run=1, target='xxx', edge='yyy_zzz', wp='water', state='stateA'):
        if runtype == 'em':
            return f'{hybPath}/{wp}/{edge}/ions{run}.pdb'
        elif runtype == 'nvt':
            return f'{runpath["em"]}/{wp}/{edge}/{state}/em{run}/em{run}.gro'
        elif runtype == 'eq': 
            return f'{runpath["nvt"]}/{wp}/{edge}/{state}/nvt{run}/nvt{run}.gro'
        elif runtype == 'morphes': 
            return f'{runpath["eq"]}/{wp}/{edge}/{state}/eq{run}/eq{run}.gro'
        else:
            print('runtype not known')
            return ''


    # workpath/[water|protein]/edge* - every edge has its own folder
    waterProtein = ['water','protein']
    # workpath/[water|protein]/edge*/state[A|B] - two states will be considered for every edge
    states = ['stateA','stateB']

    ################################################
    # energy minimization #
    ################################################
    
    # for testing set this variable to True
    # read pre-defined edges

   
    edgesToUse = read_edges( target )

    #----- TEST ------
    bTest = False
    if bTest==True:
        edgesToUse = testedges
    #----- TEST ------    

    bDeleteRunFiles = True
    
    for edge in edgesToUse.keys():
        for wp in waterProtein:
            for state in states:
                # 3 replicas
                for run in range(1,replicas+1):
                    printInfo(runtype=runtype, run=run, target=target, edge=edge, wp=wp, state=state)
                
                    # create folder
                    create_folder( f'{runpath[runtype]}/{wp}/{edge}/{state}/', f'{runtype}{run}/' )
                
                    if bDeleteRunFiles:
                        # everything is deleted!
                        toclean = glob.glob(f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/*.*')
                        for clean in toclean:
                            os.remove(clean)
                    else:
                        if os.path.isfile(f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/{runtype}{run}.log'):
                            # run log exists, so nothing is deleted
                            print('Energy minimization has been run already, nothing is deleted. Continue with next simulation.')
                            continue
                        else:
                            toclean = glob.glob(f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/*.*')
                            for clean in toclean:
                                os.remove(clean)

                    # specify input files
                    mdp = f'{mdppath}/{runtype}_{state}.mdp'
                    topology = f'{hybPath}/{wp}/{edge}/topol{run}.top'
                    coord = getRunCoord(runtype, run=run, target=target, edge=edge, wp=wp, state=state)
                    
                    # specify output files
                    tprfile = f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/{runtype}{run}.tpr' # temporary tpr file  
                    mdout = f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/mdout.mdp'
                
                    # call to gmx grompp
                    process = subprocess.Popen(['gmx','grompp',
                                                '-p',topology,
                                                '-c',coord,
                                                '-o',tprfile,
                                                '-f',mdp,
                                                '-po',mdout,
                                                '-maxwarn',str(3)],
                                               stdout=subprocess.PIPE, 
                                               stderr=subprocess.PIPE)
                    process.wait()    
            
                    if verbose:
                        out = process.communicate()
                        print('STDOUT{} '.format(out[0].decode("utf-8")))
                        print('STDERR{} '.format(out[1].decode("utf-8")))
            
                    # set variables
                    jobscriptFile = f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/{runtype}{run}.sh'
                    simpath = f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/'
                    jobname = f'{runtype}{run}_{target}_{wp}_{edge}_{state}'

                    # write submission file
                    simtime,simcpu = decide_on_resources( wp, runtype )

                    if queueType == 'sge':
                        create_SGE_jobscript( jobscriptFile, simpath, jobname, runtype, run, simtime=simtime, simcpu=simcpu )
                    elif queueType == 'slurm':
                        create_SLURM_jobscript( jobscriptFile, simpath, jobname, runtype, run, simtime=simtime, simcpu=simcpu )


def checkSimulations(targetID, target, forcefield, runtype, verbose=False):
    # for testing set this variable to True
    bTest = False
    
    #----- TEST ------
    edgesToUse = edges
    if bTest==True:
        edgesToUse = testedges
        #----- TEST ------    
    
    # set paths relative to base
    inputpath = './input'
    mdppath = './mdppath'

    # path where to find coordinates, topologies, hybrid topologies, simulations...
    ligPath = f'PLBenchmarks/data/{targetDir}/03_docked/'
    topPath = f'PLBenchmarks/data/{targetDir}/04_topo/{forcefield}/'
    hybPath = f'PLBenchmarks/data/{targetDir}/05_hybrid/{forcefield}'

    runpath =  {'em' : f'PLBenchmarks/data/{targetDir}/06_em/{forcefield}',
                'nvt': f'PLBenchmarks/data/{targetDir}/07_nvt/{forcefield}',
                'eq' : f'PLBenchmarks/data/{targetDir}/08_eq/{forcefield}',
                'morphes': f'PLBenchmarks/data/{targetDir}/09_morphes/{forcefield}'
            }
    


    # workpath/[water|protein]/edge* - every edge has its own folder
    waterProtein = ['water','protein']
    # workpath/[water|protein]/edge*/state[A|B] - two states will be considered for every edge
    states = ['stateA','stateB']

    for edge in edgesToUse.keys():
        for wp in waterProtein:
            for state in states:
                # 3 replicas
                for run in range(1,replicas+1):
                    printInfo(runtype=runtype, run=run, target=target, edge=edge, wp=wp, state=state)
                
                    if os.path.isfile(f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/{runtype}{run}.log'):
                        with open(f'{runpath[runtype]}/{wp}/{edge}/{state}/{runtype}{run}/{runtype}{run}.log') as f:
                            for line in f.readlines():
                                if re.search('Steepest Descents converged', line):
                                    print(line.strip())
                                    break
                            else:
                                print('\33[31mSimulation not (yet) converged.\33[0m')
                    else:
                        print('\33[31mNo log file available, simulation probably not finished\33[0m')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', 
                        '--target', 
                        metavar = 'TARGET',
                        type=str,
                        default='01_jnk1',
                        help='The target protein.')
    parser.add_argument('-s', 
                        '--simulationtype', 
                        metavar = 'SIMULATION_TYPE',
                        type=str,
                        default='em',
                        choices = ['em', 'eq', 'nvt', 'morphes'],
                        help='The simulation type.')
    parser.add_argument('-f', 
                        '--forcefield', 
                        metavar = 'FORCEFIELD',
                        type=str,
                        default='smirnoff99Frosst-1.1.0.offxml',
                        choices = ['smirnoff99Frosst-1.1.0.offxml', 'openff-1.0.0.offxml', 'gaff2'],
                        help='The force field used.')
    parser.add_argument('-q', 
                        '--queuetype', 
                        metavar = 'QUEUE',
                        type=str,
                        default='slurm',
                        choices = ['slurm', 'sge'],
                        help='The queue type of the job scripts.')
    parser.add_argument('-v', 
                        '--verbose', 
                        type=bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Turn on verbose output.')
    args = parser.parse_args()

    
    # forcefield
    forcefield = args.forcefield

    runtype = args.simulationtype

    prepareSimulations(args.target, forcefield, runtype, queueType=args.queuetype, verbose=args.verbose)

                

    



