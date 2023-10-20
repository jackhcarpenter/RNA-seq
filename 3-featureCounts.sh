#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=2      #  
#SBATCH --mem-per-cpu=4G       # in megabytes, unless unit explicitly stated
#SBATCH --time=6:0:0
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=carpenterj3@cardiff.ac.uk   # email address used for event notification
#SBATCH --mail-type=BEGIN,END,FAIL # email on job start, end, and/or failure

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

module load subread-2.0.0-gcc-8.3.1-l7x34bp

# define the reference directroy

export refdir=/mnt/scratch/xxxxxx/RNA-seq/reference_genome

# define the working directory

export workingdir=/mnt/scratch/xxxxxx/RNA-seq

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

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
	"TCP_M5_S15")

for i in ${list[@]}
do
    echo "============================="
    echo ${i} "markdup FC = running"

# Count how many genomic features are present in your sequencing reads (example 
# is gene reads, and only exons because it's RNA)

    featureCounts \
                -T ${SLURM_CPUS_PER_TASK} \
                -p \
                -F GTF \
                -t exon \
                -g gene_id \
                -a $refdir/Arabidopsis_thaliana.TAIR10.57.gtf \
                -o $workingdir/featureCounts/${i}.markdup.featurecount \
                $workingdir/markdup/${i}.markdup.bam

    echo ${i} "markdup FC = complete"

    echo ${i} "remove duplicate FC = running"

    featureCounts \
                -T ${SLURM_CPUS_PER_TASK} \
                -p \
                -F GTF \
                -t exon \
                -g gene_id \
                -a $refdir/YOUR_REFERENCE_ANNOTATION \
                -o $workingdir/featureCounts/${i}.rmdup.featurecount \
                $workingdir/markdup/${i}.rmdup.bam

    echo ${i} "remove duplicate FC = complete" 
    echo "============================="

done

#################################################################################
# End
#################################################################################