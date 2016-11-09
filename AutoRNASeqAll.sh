AutoRNASeqAll () {
	if [ "$#" -ne 2 ]; then
		echo "[ERROR] AutoRNASeqAll requires two argumens. Usage: AutoRNASeqAll [organism] [directory]"
		return 1 
	fi

	echo "-----------------------------------"
	source config.sh
	echo "Directory: 	"$2
	echo "-----------------------------------"

	for f in `ls -d $2/*`; do AutoRNASeq `basename $f` $1 $f/`basename $f`.fastq; done

	return 0
}

