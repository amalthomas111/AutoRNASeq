AutoRNASeqAll () {
	if [ "$#" -ne 2 ]; then
		echo "[ERROR] AutoRNASeqAll requires two argumens. Usage: AutoRNASeqAll [organism] [directory]"
		return 1 
	fi

	echo "-----------------------------------"
	source config.sh
	echo "Directory: 	"$2
	echo "-----------------------------------"

	for f in `ls -d $2/*`; do 
		numfastqfiles=`ls -1 ${f}/*.fastq | wc -l | cut -f1`
		if [  $numfastqfiles -eq 1 ]
		then
			# Single end read
			rnaseqallr1=`ls ${f}/*.fastq`
			AutoRNASeq `basename $f` $1 $rnaseqallr1; 
		elif [ $numfastqfiles -eq 2 ]
		then
			rnaseqallr1=`ls ${f}/*R1*.fastq`
			rnaseqallr2=`ls ${f}/*R2*.fastq`
			if [ -z $rnaseqallr1 ] || [ -z $rnaseqallr2 ]
			then
				echo "[ERROR] Paired-end reads detected in directory `basename $f` but reads do not end with R1 and R2"
				return 1
			else
				AutoRNASeq `basename $f` $1 $rnaseqall1 $rnaseqall2
			fi
		else
			echo "[ERROR] Directory `basename $f` doesn't have one or two fastq files"
			return 1
		fi
	done

	return 0
}

