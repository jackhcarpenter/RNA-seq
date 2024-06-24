#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=1      # for multi-threaded jobs
#SBATCH --mem-per-cpu=32G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=carpenterj3@cardiff.ac.uk      # email
#SBATCH --mail-type=BEGIN,END,FAIL      # email on job start, end, and/or failure

#################################################################################
# Print Slurm Parameters to Console
#################################################################################

echo "Usable Environment Variables:"
echo "============================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

#################################################################################
# Modulels to Load and Setup
#################################################################################

module load fastqc/v0.11.9

export workingdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/fastq/merged

echo "working dir =" $workingdir

export exportdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/fastqc/

echo "export dir =" $exportdir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMDs
#################################################################################

# Loop variables

for file in $workingdir/*
do
        if [[ $file == *.fastq ]]
        then

                echo "============================="
                echo "RUNNING fastqc"
                echo ${file} "= running"

                fastqc \
                        ${file} \
                        -o $exportdir

                echo ${file} "= complete"


            else
                    echo "============================="
                    echo "ignoring" ${file} "-  not a fastq file"

            fi
    done

echo "fastqc complete"
echo "============================="

#################################################################################
# End
#################################################################################