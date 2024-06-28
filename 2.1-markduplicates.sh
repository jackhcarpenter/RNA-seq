#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=4      #
#SBATCH --mem-per-cpu=16G       # in megabytes, unless unit explicitly stated
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

module load picard/2.22.2
module load bamtools/v2.5.1
module load samtools/1.10

# point to the directory containing the reference genome

export refdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/reference_genome
echo "Reference directory =" $refdir

# define the working directory

export workingdir=/mnt/scratch/xxxxxx/RNA-seq_TCP4_STM/STAR
echo "Working directory =" $workingdir

# define the export directory

export "Export directory =" exportdir=/mnt/scratch/xxxxxxx/RNA-seq_TCP4_STM/markdup
echo $exportdir

##REMEMBER: set up any directories that the software needs in this script in case
##it is unable to do so itself

#################################################################################
# Main CMD
#################################################################################

# Loop variables

declare -a files

for file in $workingdir/*_unsort.Log.out
do
        echo ${file}
        files+=("$(basename ${file::-15})")

done

echo ${files}

# Loops

for file in ${files[@]}
do

        echo "============================="
        echo ${file} "samtoolsort = running"

        # Sort sequences so that they are organized by genomic coordinates
        samtools sort \
        -@ ${SLURM_CPUS_PER_TASK} \
        -o $workingdir/${file}_sorted.bam \
        $workingdir/${file}_unsort.Aligned.out.bam

        samtools index \
        $workingdir/${file}_sorted.bam

        echo ${file} "samtoolsort = complete"
        echo "============================="

        echo ${file} "markdup = running"

        ## Mark the duplicated reads
        java -jar $PICARD MarkDuplicates \
                I=$workingdir/${file}_sorted.bam \
                O=$exportdir/${file}_markdup.bam \
                M=$exportdir/${file}_metrics.markdup.txt \
                REMOVE_DUPLICATES=false \
                VALIDATION_STRINGENCY=SILENT

        bamtools stats \
                -in $exportdir/${file}_markdup.bam \
                > $exportdir/${file}_markdup_dupstats.txt

        echo ${file} "markdup = complete"

        echo ${file} "remove duplicate = running"

        ## Remove the duplicated reads
        java -jar $PICARD MarkDuplicates \
                I=$workingdir/${file}_sorted.bam \
                O=$exportdir/${file}_rmdup.bam \
                M=$exportdir/${file}_metrics.rmdup.txt \
                REMOVE_DUPLICATES=true \
                VALIDATION_STRINGENCY=SILENT

        bamtools stats \
                -in $exportdir/${file}_rmdup.bam \
                > $exportdir/${file}_rmdup.dupstats.txt

        echo ${file} "remove duplicate = complete"
        echo "============================="

## Now look at the files to see if it is better to keep or remove duplicates

done

#################################################################################
# End
#################################################################################
