#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=MoMSeg
#################  
#a file for job output, you can check job progress
#SBATCH --output=MoM.out
#################
# a file for errors from the job
#SBATCH --error=MoM.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=20:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=2
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=4000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1 --cpus-per-tasks=16" instead
#SBATCH --ntasks-per-node=16
#################

#now run normal batch commands
module load matlab
matlab -nodisplay -r "mainloop('segTest_a3'), quit"

