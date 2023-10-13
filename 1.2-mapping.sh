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

export refdir=/mnt/scratch/c1831460/RNA-seq/reference_genome

# define the working directory

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

# Map forward and reverse reads to indexed genome

echo "RUNNING ALIGNMENT"

for x in ${lane[@]}
do
	for i in ${list[@]}
	do
		echo ${i} " = running"

		STAR \
			--outMultimapperOrder Random \
			--outSAMmultNmax 1 \
			--runThreadN ${SLURM_CPUS_PER_TASK} \
			--runMode alignReads \
			--outSAMtype BAM Unsorted \
			--quantMode GeneCounts \
			--outFileNamePrefix ${i}_${x}_unsort. \
			--genomeDir $refdir/ \
			--readFilesIn /mnt/scratch/c1831460/RNA-seq/fastp/${i}_${x}.fastp1.gz \
			$workingdir/fastp/${i}_${x}.fastp2.gz \
			--readFilesCommand zcat

		echo ${i} " = complete"
	done
done

echo "ALIGNMENT COMPLETE"
echo "============================="

#################################################################################
# End
#################################################################################