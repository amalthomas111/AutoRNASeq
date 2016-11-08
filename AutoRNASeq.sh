function AutoRNASeq ()  {

	# Loads config variables
	echo "----------------------------------"
	source config.sh $1 $2 $3 $4
	echo "----------------------------------"
	proj=$1
	genome=$2
	r1=$3
	r2=$4
	readdir=`dirname $r1`
	outdir="${readdir}/output-${proj}"
	tempdir="${readdir}/_temp-${proj}"
	ref_genome=/home/cmb-panasas2/desenabr/ref_genomes/${genome}/STAR_index
	gtf=/home/cmb-panasas2/desenabr/ref_genomes/${genome}/bowtie_index/${genome}_.gtf

	# Check if number of arguments is correct and things exist
	if [ "$#" -lt 3 ]; then
		echo "[ERROR] insufficient arguments"
		echo "Usage:	./AutoRNASeq.sh [project-name] [ref-genome] [fastq1] [optional-fastq2]"
		return 1
	fi

	# Reference genome
	if [ ! -d /home/cmb-panasas2/desenabr/ref_genomes/${organism} ]; then
		echo "[ERROR] reference genome directory not found : ${organism}"
		return 1
	fi

	# Fastq1 
	if [ ! -f $r1 ]; then
		echo "[ERROR] fastq read 1 file not found: $r1"
		return 1
	fi

	# Fastq2 (if given)
	if [ "$#" -eq 4 ] && [ ! -f $r2 ]; then
		echo "[ERROR] fastq read 2 file not found: $2"
		return 1
	fi

	echo "----------------------------------"

	echo "Project name: 		$proj"
	echo "Reference Genome: 	$genome"
	if [ "$#" -eq 4 ]
	then
		echo "Read type: 		paired-end"
		echo "Read 1: 		$r1"
		echo "Read 2: 		$r2"
	else
		echo "Read type: 		single-end"
		echo "Read: 			$r1"
	fi
	echo "----------------------------------"

	# creates directories
	if [ -d $outdir ]
	then
		rm -rf $outdir
	fi

	if [ -d $tempdir ]
	then
		rm -rf $tempdir
	fi

	mkdir -p $tempdir
	mkdir -p $outdir

	# fastqc
	fastqcjob=$(./jobs/1-fastqc.sh $PWD $proj $readdir $tempdir $outdir $r1 $r2)
	echo "fastqc job : 		$fastqcjob"

	# trim_galore
	trimjob=$(./jobs/2-trim.sh $fastqcjob $PWD $proj $tempdir $outdir $r1 $r2)
	echo "trim galore id : 	$trimjob"	

	# STAR mapping
	if [ "$#" -eq 3 ]; then
		mapjob=$(./jobs/3-map-single.sh $trimjob $PWD $proj $tempdir $ref_genome $r1)
		echo "STAR id (single end) : 	$mapjob"	
	else
		mapjob=$(./jobs/3-map-paired.sh $trimjob $PWD $proj $tempdir $ref_genome $r1 $r2)
		echo "STAR id (paired end): 	$mapjob"	
	fi

	# htseq-count
	countjob=$(./jobs/4-count.sh $mapjob $PWD $proj $tempdir $outdir $gtf)
	echo "count job: 		$countjob"	

	# Postprocessing job
	postjob=$(./jobs/5-postprocessing.sh $countjob $PWD $proj $tempdir $outdir)
	echo "postprocessing id: 	$postjob"
}

AutoRNASeq $1 $2 $3 $4
