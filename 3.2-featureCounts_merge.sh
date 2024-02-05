#!/bin/bash

#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=1      #   
#SBATCH --mem-per-cpu=1G     # in megabytes, unless unit explicitly stated
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
# Setup
#################################################################################

# defining the working directory

export workingdir=/mnt/scratch/xxxxxx/RNA-seq/featureCounts

dup=("markdup" \
	"rmdup")

# the featureCount output format needs to be converted into tables containing 
# just gene IDs and corresponding counts for each

# extract the first and last column from the featurecounts.csv files (skipping
# the header)

for x in ${dup[@]}

do 
	echo "============================="
	echo "${x}"
	
	for i in {1..24}

	do
		echo "S"${i}" = running"
		awk 'NR>1' $workingdir/*S${i}.${x}.featurecounts.txt | \
		       	awk '{print $1, $7}' FS='\t' \
			>> joined/S${i}.${x}.txt
		echo ${i} " = complete"
	done

done
echo ${x} " first and last column extraction  = complete"
echo "============================="S13.markdup.bam


# using join to merge STM and TCP4 featureCount files into gene-specific 
# dataframes/files that only contains the counts field which was just extracted 

for x in ${dup[@]}
do
	# create the initial merge file
	
	echo "${x}"
	join \
		$workingdir/joined/S1.${x}.txt \
		$workingdir/joined/S2.${x}.txt \
		>> $workingdir/joined/temp.S2.${x}.txt

	# now add all second columns from the other files to the dataframe

	sum=1
	
        for i in {3..12}

        do
            echo "joining S"${i} "to STM_"${x}

		    previous=$(($i-$sum))
		
		# this method creates a lot of files for join to work properly
		# so the files created are marked as temp. to be removed later
                
		join \
			joined/temp.S${previous}.${x}.txt \
			joined/S${i}.${x}.txt \
			>> joined/temp.S${i}.${x}.txt
	
	done

	echo "Converting to tab delimited format"

	sed \
		's/ /\t/g' \
		$workingdir/joined/temp.S12.${x}.txt \
		>> $workingdir/joined/STM.${x}.csv

	echo "============================="
	
done

# Now do the same for TCP4

for x in ${dup[@]}
do
    # create the initial merge file
        
	echo "${x}"
	join \
		$workingdir/joined/S13.${x}.txt \
		$workingdir/joined/S14.${x}.txt \
		>> $workingdir/joined/temp.${i}.${x}.txt

	# now add data from the other files to it
        
	sum=13

	for i in {15..24}

        do
            echo "joining S"${i} "to TCP4_"${x}

            previous=$(($i-$sum))

            join  \
		    joined/temp.S${previous}.${x}.txt \
		    joined/S${i}.${x}.txt \
		    >> joined/temp.S${i}.${x}.txt
	done

	echo "Converting to tab delimited format"

	sed \
		's/ /\t/g' \
		$workingdir/joined/temp.S24.${x}.txt \
		>> $workingdir/joined/TCP4.${x}.csv

	echo "============================="
	
done

# removing unwanted files

echo "cleaning up"

#rm $workingdir/joined/temp.*.*.csv

