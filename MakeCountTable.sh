MakeCountTable() {
	numdirectories=`ls -1 | wc -l | cut -f1`
	numvectors=`ls -1 */output*/*.vector | wc -l | cut -f1`
	
	# Verbose
	echo "Num directories : $numdirectories"
	echo "Num vectors : $numvectors"

	if [ $numdirectories  -ne $numvectors ]
	then
		echo "[ERROR] not all directories have count vectors"
		return 1
	fi

	# Make row names
	printf "\t`ls -d1 */ | cut -f1 -d '/' | tr '\n' '\t'`\n" >_rownames.txt

	# Make count table with no title
	paste $PANASAS/ref_genomes/$1/${1}.genenames */output*/*.vector >count-table-nonames.txt

	# Merges both
	cat _rownames.txt count-table-nonames.txt >count-table-`basename $PWD`.txt

	# Remove temporary files
	rm _rownames.txt count-table-nonames.txt
	
	return 0

}
