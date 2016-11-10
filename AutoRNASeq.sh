function AutoRNASeqFastQC () { 
    proj=$2;
    readdir=$3;
    tempdir=$4;
    outdir=$5;
    r1=$6;
    r2=$7;
    error=$tempdir"/1-fastqc.error";
    output=$tempdir"/1-fastqc.output";
    cat  <<EOS | qsub -
    #!/bin/sh
    #PBS -N fqc-${proj}
    #PBS -S /bin/sh
    #PBS -e $error
    #PBS -o $output
    #PBS -l mem=30000mb
    #PBS -l vmem=30000mb
    #PBS -l pmem=30000mb
    #PBS -l nodes=1:ppn=1
    #PBS -l walltime=10:00:00
    #PBS -q cmb

    cd ${1}
    $run_fastqc $r1 $r2
    mv ${readdir}/*.html ${readdir}/*.zip $outdir

EOS
}

function AutoRNASeqTrim () { 
    proj=$3;
    tempdir=$4;
    outdir=$5;
    r1=$6;
    r2=$7;
    fqbase=${2}/$(cut -d "." -f 1 <<< `basename $r1`);
    error=$tempdir"/2-trimgalore.error";
    output=$tempdir"/2-trimgalore.output";
    cat  <<EOS | qsub -
    #!/bin/sh
    #PBS -N trim-${proj}
    #PBS -S /bin/sh
    #PBS -e $error
    #PBS -o $output
    #PBS -W depend=afterok:${1}
    #PBS -l mem=50000mb
    #PBS -l vmem=50000mb
    #PBS -l pmem=50000mb
    #PBS -l nodes=1:ppn=1
    #PBS -l walltime=20:00:00
    #PBS -q cmb

    cd ${2}
    if [ -z $r2 ]
    then
        $run_trim_galore --length 10 \
                         --path_to_cutadapt $run_cutadapt \
                          --output_dir $tempdir $r1
 
    else
        $run_trim_galore --length 10 \
                         --path_to_cutadapt $run_cutadapt \
                         --output_dir $tempdir \
                         --paired $r1 $r2 
    fi

EOS
}

function AutoRNASeqMapSingle () { 
    proj=$3;
    tempdir=$4;
    ref_genome=$5;
    r1=$6;
    fq=$tempdir/$(cut -d "." -f 1 <<< `basename $r1`)_trimmed.fq;
    error=$tempdir"/3-map.error";
    output=$tempdir"/3-map.output";
    suffix="${tempdir}/${proj}.";
    cat  <<EOS | qsub -
    #!/bin/sh
    #PBS -N map-${proj}
    #PBS -S /bin/sh
    #PBS -e $error
    #PBS -o $output
    #PBS -W depend=afterok:${1}
    #PBS -l mem=120000mb
    #PBS -l vmem=120000mb
    #PBS -l pmem=120000mb
    #PBS -l nodes=1:ppn=1
    #PBS -l walltime=40:00:00
    #PBS -q cmb

    cd ${2}
    $run_star --outFileNamePrefix $suffix \
              --quantMode GeneCounts \
              --outSAMtype BAM SortedByCoordinate \
              --bamRemoveDuplicatesType UniqueIdentical \
              --runThreadN 1 \
              --outBAMsortingThreadN 1 \
              --genomeDir $ref_genome \
              --readFilesIn $fq
EOS
}

function AutoRNASeqMapPaired () { 
    proj=$3;
    tempdir=$4;
    ref_genome=$5;
    r1=$6;
    r2=$7;
    error=$tempdir"/3-map.error";
    output=$tempdir"/3-map.output";
    suffix="${tempdir}/${proj}.";
    cat  <<EOS | qsub -
    #!/bin/sh
    #PBS -N map-${proj}
    #PBS -S /bin/sh
    #PBS -e $error
    #PBS -o $output
    #PBS -W depend=afterok:${1}
    #PBS -l mem=120000mb
    #PBS -l vmem=120000mb
    #PBS -l pmem=120000mb
    #PBS -l nodes=1:ppn=1
    #PBS -l walltime=40:00:00
    #PBS -q cmb

    cd ${2}
    $run_star --outFileNamePrefix $suffix \
              --quantMode GeneCounts \
              --outSAMtype BAM SortedByCoordinate \
              --bamRemoveDuplicatesType UniqueIdentical \
              --runThreadN 1 \
              --outBAMsortingThreadN 1 \
              --genomeDir $ref_genome \
              --readFilesIn $tempdir/$(cut -d "." -f 1 <<< `basename $r1`)_val_1.fq \
                            $tempdir/$(cut -d "." -f 1 <<< `basename $r2`)_val_2.fq
EOS
}

function AutoRNASeqCount () { 
    proj=$3;
    tempdir=$4;
    outdir=$5;
    gtf=$6;
    error=$tempdir"/4-count.error";
    output=$tempdir"/4-count.output";
    cat  <<EOS | qsub -
    #!/bin/sh
    #PBS -N count-${proj}
    #PBS -S /bin/sh
    #PBS -e $error
    #PBS -o $output
    #PBS -W depend=afterok:${1}
    #PBS -l mem=30000mb
    #PBS -l vmem=30000mb
    #PBS -l pmem=30000mb
    #PBS -l nodes=1:ppn=1
    #PBS -l walltime=10:00:00
    #PBS -q cmb

    cd $2
    samtools view $tempdir/${proj}.Aligned.sortedByCoord.out.bam | \
                  $run_htseq --stranded=no - $gtf >$tempdir/${proj}.counts
    
    cut -f2 $tempdir/${proj}.counts | head -n -5 > $outdir/${proj}.vector
EOS
}


function AutoRNASeqPostprocessing () { 
    proj=$3;
    tempdir=$4;
    outdir=$5;
    error=$tempdir"/5-postprocessing.error";
    output=$tempdir"/5-postprocessing.output";
    cat  <<EOS | qsub -
    #!/bin/sh
    #PBS -N postprocessing-${proj}
    #PBS -S /bin/sh
    #PBS -e $error
    #PBS -o $output
    #PBS -W depend=afterok:${1}
    #PBS -l mem=5000mb
    #PBS -l vmem=5000mb
    #PBS -l pmem=5000mb
    #PBS -l nodes=1:ppn=1
    #PBS -l walltime=10:00:00
    #PBS -q cmb

    cd $2
    cat $tempdir/*.error >$outdir/readme.log
    cp $tempdir/*.final.out $outdir/${proj}.STAR-Mapping-Log
    # rm -rf $tempdir
EOS
}

AutoRNASeq () { 
    if [ "$#" -lt 3 ]; then
        echo "[ERROR] insufficient arguments";
        echo "Usage:    ./AutoRNASeq.sh [project-name] [ref-genome] \
                                        [fastq1] [optional-fastq2]";
        return 1;
    fi;

    # Loads variables
    proj=$1;
    genome=$2;
    r1=`realpath $3`;

    # 4th argument means it's paired end read
    if [ "$#" -eq 4 ]; then
        r2=`realpath $4`;
    fi;

    # Outpit/temp directory paths
    readdir=`dirname $r1`;
    outdir="${readdir}/output-${proj}";
    tempdir="${readdir}/_temp-${proj}";

    # Reference genome index and annotation
    ref_genome=/home/cmb-panasas2/desenabr/ref_genomes/${genome}/STAR_index;
    gtf=/home/cmb-panasas2/desenabr/ref_genomes/${genome}/bowtie_index/${genome}_.gtf;

    # Check if files exist
    if [ ! -d /home/cmb-panasas2/desenabr/ref_genomes/${genome} ]; then
        echo "[ERROR] reference genome directory not found : ${genome}";
        return 1;
    fi;
    if [ ! -f $r1 ]; then
        echo "[ERROR] fastq read 1 file not found: $r1";
        return 1;
    fi;
    if [ "$#" -eq 4 ] && [ ! -f $r2 ]; then
        echo "[ERROR] fastq read 2 file not found: $2";
        return 1;
    fi;

    # Verbose
    echo "----------------------------------";
    echo "Project name:         $proj";
    echo "Reference Genome:     $genome";
    if [ "$#" -eq 4 ]; then
        echo "Read type:         paired-end";
        echo "Read 1:         $r1";
        echo "Read 2:         $r2";
    else
        echo "Read type:         single-end";
        echo "Read:             $r1";
    fi;

    # Creates empty directories
    if [ -d $outdir ]; then
        rm -rf $outdir;
    fi;
    if [ -d $tempdir ]; then
        rm -rf $tempdir;
    fi;
    mkdir -p $tempdir;
    mkdir -p $outdir;
    echo "";

    fastqcjob=$(AutoRNASeqFastQC $PWD $proj $readdir $tempdir $outdir $r1 $r2);
    echo "fastqc job :         $fastqcjob";

    trimjob=$(AutoRNASeqTrim $fastqcjob $PWD $proj $tempdir $outdir $r1 $r2);
    echo "trim galore id :     $trimjob";

    if [ "$#" -eq 3 ]; then
        mapjob=$(AutoRNASeqMapSingle $trimjob $PWD $proj $tempdir $ref_genome $r1);
        echo "STAR id (single end) :     $mapjob";
    else
        mapjob=$(AutoRNASeqMapPaired $trimjob $PWD $proj $tempdir $ref_genome $r1 $r2);
        echo "STAR id (paired end):     $mapjob";
    fi;

    countjob=$(AutoRNASeqCount $mapjob $PWD $proj $tempdir $outdir $gtf);
    echo "count job:         $countjob";

    postjob=$(AutoRNASeqPostprocessing $countjob $PWD $proj $tempdir $outdir);

    echo "postprocessing id:     $postjob";
    echo "----------------------------------";
    return 0
}

AutoRNASeqAll () 
{ 
    if [ "$#" -ne 2 ]; then
        echo "[ERROR] AutoRNASeqAll requires two argument. \
              Usage: AutoRNASeqAll [organism] [directory]";
        return 1;
    fi;

    projectpath=`realpath $2`;

    echo "======================================";
    echo "AutoRNASeqAll Directory:     "$projectpath;
    echo "======================================";

    for f in `ls -d ${projectpath}/*`;
    do
        numfastqfiles=`ls -1 ${f}/*.fastq | wc -l | cut -f1`;
        if [ $numfastqfiles -eq 1 ]; then
            rnaseqallr1=`ls ${f}/*.fastq`;
            AutoRNASeq `basename $f` $1 $rnaseqallr1;
        else if [ $numfastqfiles -eq 2 ]; then
            rnaseqallr1=`ls ${f}/*R1*.fastq`;
            rnaseqallr2=`ls ${f}/*R2*.fastq`;
            if [ -z $rnaseqallr1 ] || [ -z $rnaseqallr2 ]; then
                echo "[ERROR] Paired-end reads detected in directory \
                     `basename $f` but reads do not end with R1 and R2";
                return 1;
            else
                 AutoRNASeq `basename $f` $1 $rnaseqallr1 $rnaseqallr2;
            fi;
            else
                echo "[ERROR] Directory `basename $f` \
                      doesn't have one or two fastq files";
                return 1;
            fi;
        fi;
    done;
    return 0
}

MakeCountTable ()  { 
    numdirectories=`ls -d1 */ | wc -l | cut -f1`;
    numvectors=`ls -1 */output*/*.vector | wc -l | cut -f1`;

    echo "Num directories : $numdirectories";
    echo "Num vectors : $numvectors";

    # Checks if all output subdirectories have a *.vector file
    if [ $numdirectories -ne $numvectors ]; then
        echo "[ERROR] not all directories have count vectors";
        return 1;
    fi;

    # Creates row names based on directory titles
    printf "\t`ls -d1 */ | cut -f1 -d '/' | tr '\n' '\t'`\n" > _rownames.txt;

    # Merges vectors
    paste $PANASAS/ref_genomes/$1/${1}.genenames */output*/*.vector > count-table-nonames.txt;

    # Concatenates row names and vectors
    cat _rownames.txt count-table-nonames.txt > count-table-`basename $PWD`.txt;

    # Remove temporary files
    rm _rownames.txt count-table-nonames.txt;
    return 0
}

