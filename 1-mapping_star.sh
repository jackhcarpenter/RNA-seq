#!/bin/bash

#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=8      #   
#SBATCH --mem-per-cpu=4G     # in megabytes, unless unit explicitly stated
#SBATCH --time=48:00:00
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=USERNAME@INSTITUTIONAL_ADDRESS  # email
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

module load STAR/2.7.3a

# point to the directory containing the reference genome

export $refdir=/your/reference/directory

# define the working directory

export $workingdir=/your/working/directory

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

mkdir star

#################################################################################
# Main CMDs
#################################################################################

echo "============================="
echo "RUNNING INDEXING"

STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --limitGenomeGenerateRAM \
    --rnMode genomeGenerage \
    --genomeDir $refdir/ \
    --genomeFastaFiles $redir/REFERENCE_GENOME.dna.toplevel.fa \
    --sjdbGTFfile $refdir/YOUR_REFERENCE_ANNOTATION.gtf \
    --sjdbOverhang 75

# Note: Change --sjdbOverhang to length of your sequence data/2 minus

echo "INDEXING COMPLETE"
echo "============================="
# List of sequences to align

list=("sample1" "sample2" "samplen")

# Map forward and reverse reads to indexed genome

echo "RUNNING ALIGNMENT"

for i in ${list[@]}
do
    echo ${i} " = running"

    STAR \
        --outMultimapperOrder Random \
        --outSAMmultNmax 1 \
        --runThreadsN ${SLURM_CPUS_PER_TASK} \
        --runMode alignReads \
        --outSAMtype BAM Unsorted \
        --quantMode GeneCounts \
        --outFileNamePrefix $workingdir/star/${i}_unsort. \
        --genomeDir $refdir \
        --readFilesIn $workingdir/trimmed/${i}_fp1.fastq.gz \
                    $workingdir/trimmed/${i}_fp2.fastq.gz \
        --readFilesCommand zcat

    echo ${i} " = complete"
done

echo "ALIGNMENT COMPLETE\n" "============================="

#################################################################################
# End
#################################################################################