#!/bin/bash

# Based on the SeqCap pipline from:
# Targeted capture of complete coding regions across divergent species
# Ryan K Schott, Bhawandeep Panesar, Daren C Card, Matthew Preston, Todd A Castoe, Belinda SW Chang

clean_up() {
	#Deletes broken files if they were created and ends program
	printf "Creating $1" >&2
	rm -f $1
	shift
	for file in $@
	do
		printf ", $file" >&2
		rm -f $file
	done
	printf " failed\n" >&2
	exit 1
}

USAGE="Usage: bash BWA.sh [-s scriptdir] ref_genome.fas threads [fwd_read.fastq rvrse_read.fastq]..."
scriptdir=$(pwd)
seqs=1

while getopts ':hs:' flag
do
	case $flag in
		h) # Help
			echo "$USAGE"
			exit 1
			;;
		s) # Change script directory
			scriptdir=$OPTARG
			;;
		\?)# Unexpected flag
			echo "Unexpected option -$OPTARG" >&2
			exit 1
			;;
		:) # Missing argument
			echo "Option -$OPTARG requires an argument" >&2
			exit 1
			;;
		*) # Unexpected error
		    echo "Unexpected error in getopts" >&2
			exit 1
			;;
	esac
done
shift $( expr $OPTIND - 1 )

if [ $# -le 1 ]
then
	echo $USAGE
	exit 1
fi

date +"Time started: %F %T"
echo "Using Reference Genome: $1"
ref=$1
shift

echo "Using $1 threads"
threads=$1
shift

if [ -f $ref.bwt ]
then
   echo "Reference index exists, skipping reference index build"
else
   echo "Making bwa index"
   bwa index $ref || clean_up $ref.bwt
fi

if [ $# -gt 0 ]
then
	tmp=$( dirname $ref | cut -d / -f 1 )
	dirs=("." "..")
	dirs+=( $(ls -F | grep ".*/" | sed 's|/||') )
	for d in ${dirs[@]}
	do
		if [ $tmp = $d ]
		then
			ref="../../$ref"
			break
		fi
	done

	tmp=$(basename $ref)
	refname=${tmp%.*}
	echo "Making ${refname}_BWA directory"
	if [ ! -d ${refname}_BWA ]
	then
		mkdir ${refname}_BWA
	fi
	cd ${refname}_BWA || { echo "Changing Directory failed"; exit 1; }
fi

while test ${#} -gt 0
do
	date +"Pair start: %T"
	echo "Making Directories"

	tmp=$(basename $1)
	seqname=${tmp%_R*}
	echo "Using Sequence Name: $seqname"
	if [ ! -d "$seqname" ]
	then
		mkdir $seqname
	fi

	forward=$1
	echo "Forward read: $forward"
	shift
	reverse=$1
	echo "Reverse read: $reverse"
	shift
	echo "Changing to directory $seqname"
	cd $seqname || { echo "Changing Directory failed"; exit 1; }

	if [ -f aligned_out_$seqname.sam ]
	then
	   echo "BWA output exists, skipping alignment"
	else
	   echo "Running BWA: bwa mem -B 2 -M -t $threads $ref $forward $reverse"
	   bwa mem -B 2 -M -t $threads $ref $forward $reverse > aligned_out_$seqname.sam || clean_up aligned_out_$seqname.sam
	fi

	if [ -f aligned_out_$seqname.bam ]
	then
	   echo "BAM exists, skipping conversion."
	else
	   echo "Converting SAM to BAM"
	   samtools view -Sb aligned_out_$seqname.sam > aligned_out_$seqname.bam || clean_up aligned_out_$seqname.bam
	fi

	if [ -f stats_$seqname.txt ]
	then
	   echo "Stats output exists, skipping"
	else
	   echo "Generating Stats"
	   samtools flagstat aligned_out_$seqname.bam > stats_$seqname.txt || clean_up stats_$seqname.txt
	fi

	if [ -f aligned_out_mapped_$seqname.bam ]
	then
	    echo "Mapped output exists, skipping"
	else
	    echo "Filtering out only the mapped seqs"
	    samtools view -b -F 4 -q 3 aligned_out_$seqname.bam > aligned_out_uniqmapped_$seqname.bam || clean_up aligned_out_uniqmapped_$seqname.bam
	fi

	if [ -f aligned_out_sorted_$seqname.bam ]
	then
	    echo "Sorted output exists, skipping"
	else
	    echo "Sorting bam file"
	    samtools sort aligned_out_uniqmapped_$seqname.bam aligned_out_uniqmapped_sorted_$seqname || aligned_out_uniqmapped_sorted_$seqname.bam
	fi

	if [ -f aligned_out_sorted_$seqname.bam.bai ]
	then
	    echo "Index exists, skipping"
	else
	    echo "Making samtools index for bam file"
	    samtools index aligned_out_uniqmapped_sorted_$seqname.bam || aligned_out_uniqmapped_sorted_$seqname.bam.bai
	fi

	if [ -f depth_of_coverage_$seqname.csv ]
	then
		echo "Depth of Coverage exists, skipping"
	else
		echo "Calculating depth of coverage"
		samtools idxstats aligned_out_uniqmapped_sorted_$seqname.bam > readcount_$seqname.csv
		python $scriptdir/depth_of_coverage_impl.py -b aligned_out_uniqmapped_sorted_$seqname.bam -o depth_of_coverage_$seqname.csv -t || clean_up depth_of_coverage_$seqname.csv
	fi

	date +"Pair end: %T"
	echo "_________________________________________________"
	echo "FINISHED. Moving on to next pair"
	echo "_________________________________________________"
	cd ..  || { echo "Changing Directory failed"; exit 1; }
	seqs=$((seqs+1))
done

echo "All Done here!"
date +"Time ended: %F %T"
