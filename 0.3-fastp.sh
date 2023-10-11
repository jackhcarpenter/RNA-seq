#!/bin/bash

#SBATCH --partition= defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=4G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user= carpenterj3@cardiff.ac.uk      # email
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

module load fastp/v0.20

export workingdir=/mnt/scratch/c1831460/RNA-seq

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

#################################################################################
# Main CMDs
#################################################################################

# Loop variables

lane=("L001" \
	"L002")

list=("STM_C3_S7" \
	"STM_C4_S8" \
	"STM_C5_S9" \
	"STM_CD3_S10" \
	"STM_CD4_S11" \
	"STM_CD5_S12" \
	"STM_D3_S4" \
	"STM_D4_S5" \
	"STM_D5_S6" \
	"STM_M3_S1" \
	"STM_M4_S2" \
	"STM_M5_S3" \
	"TCP_C2_S19" \
	"TCP_C4_S20" \
	"TCP_C5_S21" \
	"TCP_CO2_S22" \
	"TCP_CO4_S23" \
	"TCP_CO5_S24" \
	"TCP_D2_S16" \
	"TCP_D4_S17" \
	"TCP_D5_S18" \
	"TCP_M2_S13" \
	"TCP_M4_S14" \
	"TCP_M5_S15" \
	"Undetermined_S0")

# Trim low quality reads, remove adapters, and poly Gs

echo "RUNNING fastp"

for x in ${lane[@]}
do
    for i in ${list[@]}
    do
        echo ${i} "= running"

        fastp \
            -i $workingdir/${i}_${x}.R1.fastq.gz \
			-I $workingdir/${i}_${x}.R2.fastq.gz \
            --detect_adapted_for_pe \
            --trim_poly_g \
            --correction \
            -o $workingdir/trimmed_reads/${i}_${x}.fp1.gz \
            -O $workingdir/trimmed_reads/${i}_${x}.fp2.gz

    echo ${i} "= complete"

    done
done

echo "fastp COMPLETE" 
echo "============================="