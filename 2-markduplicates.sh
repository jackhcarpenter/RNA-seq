#!/bin/bash

#SBATCH --partition=QUEUE_NAME       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=4      #  
#SBATCH --mem-per-cpu=8G       # in megabytes, unless unit explicitly stated
#SBATCH --time=6:0:0
#SBATCH --error=logs/%J.err         # redirect stderr to this file
#SBATCH --output=logs/%J.out        # redirect stdout to this file
#SBATCH --mail-user=USERNAME@INSTITUTIONAL_ADDRESS # email address used for event notification
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

module load picard/2.22.2
module load bamtools/v2.5.1
module load samtools/1.10

# point to the directory containing the reference genome

export $refdir=/your/reference/directory

# define the working directory

export $workingdir=/your/working/directory

##REMEMBER: set up any directories that the software needs in this script in case 
##it is unable to do so itself

mkdir markdup

#################################################################################
# Main CMD
#################################################################################

list=(
        "sample_1" \
        "sample_2"\
        "sample_n")

for i in ${list[@]}

do
    echo "============================="
    echo ${i} "markdup = running"

## Mark the duplicated reads
        java -jar $PICARD MarkDuplicates \
        I=$workingdir/bowtie/${i}.sorted.bam \
        O=$workingdir/markdup/${i}.markdup.bam \
        M=$workingdir/markdup/${i}.metrics.markdup.txt \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=SILENT

        bamtools stats \
        -in $workingdir/markdup/${i}.markdup.bam \
        > $workingdir/markdup/${i}.markdup.dupstats.txt

    echo ${i} "markdup = complete"

    echo ${i} "remove duplicate = running"

## Remove the duplicated reads
        java -jar $PICARD MarkDuplicates \
        I=$workingdir/bowtie/${i}.sorted.bam \
        O=$workingdir/markdup/${i}.rmdup.bam \
        M=$workingdir/markdup/${i}.metrics.rmdup.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT

        bamtools stats \
        -in $workingdir/markdup/${i}.rmdup.bam \
        > $workingdir/markdup/${i}.rmdup.dupstats.txt

    echo ${i} "remove duplicate = complete"
    echo "============================="

## Now look at the files to see if it is better to keep or remove duplicates
done

#################################################################################
# End
#################################################################################