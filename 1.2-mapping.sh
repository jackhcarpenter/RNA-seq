#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=8      #
#SBATCH --mem-per-cpu=4G     # in megabytes, unless unit explicitly stated
#SBATCH --time=48:00:00
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=carpenterj3@cardiff.ac.uk  # email
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

module load STAR/2.7.6a

# point to the directory containing the reference genome

export refdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/reference_genome

# define the working directory

export workingdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/fastp

# define the export directory

export exportdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/STAR

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMDs
#################################################################################

# Loop variables

declare -a files

for file in $workingdir/*
do
        if [[ $file == *R1.fastp ]]
        then

                files+=("$(basename ${file::-9})")

        fi

done

echo $files

# Map forward and reverse reads to indexed genome

echo "RUNNING ALIGNMENT"

for file in ${files[@]}
do

        echo ${file} " = running"

        STAR \
                --outMultimapperOrder Random \
                --outSAMmultNmax 1 \
                --runThreadN ${SLURM_CPUS_PER_TASK} \
                --runMode alignReads \
                --outSAMtype BAM Unsorted \
                --quantMode GeneCounts \
                --outFileNamePrefix $exportdir/${file}_unsort. \
                --genomeDir $refdir/ \
                --readFilesIn $workingdir/${file}_R1.fastp \
                $workingdir/${file}_R2.fastp

        echo ${file} " = complete"
done

echo "ALIGNMENT COMPLETE"
echo "============================="

#################################################################################
# End
#################################################################################
