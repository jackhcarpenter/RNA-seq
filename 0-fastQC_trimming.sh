#!/bin/bash

#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     # for parallel distributed jobs
#SBATCH --cpus-per-task=4      # for multi-threaded jobs
#SBATCH --mem-per-cpu=4G      # in megabytes, unless unit explicitly stated
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=USERNAME@INSTITUTIONAL_ADDRESS      # email
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
module load fastp/v0.20     ## NOTE: if you are unable to disable paralisation
module load multiqc/v/1.9   ## then you must run fastp and multiqc seperatly as
                            ## they use conflicting versions of python

export workingdir=/your/working/dir

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

mkdir fastqc
mkdir fastqc/trimmed_reads

#################################################################################
# Main CMDs
#################################################################################

# List of sequences
list=("sample1" "sample2" "samplen")

# fastqc of paired end data

echo "============================="
echo "RUNNING fastqc and multiqc"

for i in ${list[@]}
do
    echo ${i} "= running"

    fastqc \
        $workingdir/${i}.R1.fastq.gz \
        -o $workingdir/fastqc/
    fastqc $workingdir/${i}.R2.fastq.gz \
        -o $workingdir/fastqc/

    echo ${i} "= complete"
done

echo "fastqc COMPLETE"

# Summarise the WC data for all the sequences

multiqc \
    -i "PROJECT_NAME_RAW_SEQUENCES" \
    $workingdir/fastqc

echo "multiqc COMPLETE" 
echo "============================="

# Trim low quality reads, remove adapters, and poly Gs

echo "RUNNING fastp"

for i in ${list[@]}
do
    echo ${i} "= running"

    fastp \
        -i $workingdir/${i}.R1.fastq.gz \
        -I $workingdir/${i}.R2.fastq.gz \
        --detect_adapted_for_pe \
        --trim_poly_g \
        --correction \
        -o $workingdir/trimmed_reads${i}.fp1.gz \
        -O $workingdir/trimmed_reads${i}.fp2.gz

   echo ${i} "= complete"

done

echo "fastp COMPLETE" 
echo "============================="

## Do another round of QC

echo "RUNNING fastqc and multiqc"

# fastqc of paired end data

for i in ${list[@]}
do
    echo ${i} "= running"

    fastqc $workingdir/trimmed_reads${i}.fp1.gz
    fastqc $workingdir/trimmed_reads${i}.fp2.gz

    echo ${i} "= complete"
done

echo "fastqc COMPLETE"

# Summarise the WC data for all the sequences

multiqc \
    -i "PROJECT_NAME_RAW_SEQUENCES_trimmed"

echo "multiqc COMPLETE" 
echo "============================="
#################################################################################
# End
#################################################################################