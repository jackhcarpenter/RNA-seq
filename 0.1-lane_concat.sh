#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=8      # for multi-threaded jobs
#SBATCH --mem-per-cpu=8G      # in megabytes, unless unit explicitly stated
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

export workingdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/fastq

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMDs
#################################################################################

## Creating file name list

list=("Col-0_11d_1_S1" \
        "Col-0_11d_2_S2" \
        "Col-0_11d_3_S3" \
        "Col-0_12h_1_S19" \
        "Col-0_12h_2_S20" \
        "Col-0_12h_4_S21" \
        "Col-0_24h_1_S28" \
        "Col-0_24h_3_S29" \
        "Col-0_24h_4_S30" \
        "Col-0_3h_1_S10" \
        "Col-0_3h_2_S11" \
        "Col-0_3h_3_S12" \
        "STM_11d_1_S4" \
        "STM_11d_2_S5" \
        "STM_11d_3_S6" \
        "STM_12h_1_S22" \
        "STM_12h_2_S23" \
        "STM_12h_3_S24" \
        "STM_24h_2_S31" \
        "STM_24h_3_S32" \
        "STM_24h_5_S33" \
        "STM_3h_1_S13" \
        "STM_3h_2_S14" \
        "STM_3h_3_S15" \
        "TCP_11d_1_S7" \
        "TCP_11d_2_S8" \
        "TCP_12h_1_S25" \
        "TCP_12h_2_S26" \
        "TCP_12h_4_S27" \
        "TCP_24h_1_S34" \
        "TCP_24h_3_S35" \
        "TCP_24h_4_S36" \
        "TCP_3h_1_S16" \
        "TCP_3h_2_S17" \
        "TCP_3h_3_S18" \
        "TCP4_11d_5_S9")

## Merging command

for i in ${list[@]}
do
        # combining lane 1 and lane to for forward reads
        echo "==========================================="
        echo "merging" ${i}"_L001_R1 and "${i}"_L002_R1"
        cat $workingdir/${i}_L001_R1* $workingdir/${i}_L002_R1* >> $workingdir/merged/${i}_R1.fastq
        # combining lane 1 and lane to for reverse reads
        echo "merging" ${i}"_L001_R2 and "${i}"_L002_R2"
        cat $workingdir/${i}_L001_R2* $workingdir/${i}_L002_R2* >> $workingdir/merged/${i}_R2.fastq

done
echo "==========================================="
echo "Complete"