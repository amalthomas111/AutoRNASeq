#!/bin/sh
proj=$3
tempdir=$4
ref_genome=$5
r1=$6
r2=$7
error=$tempdir"/3-map.error"
output=$tempdir"/3-map.output"
suffix="${tempdir}/${proj}."

cat <<EOS | qsub -
#!/bin/sh
#PBS -N map-${proj}
#PBS -S /bin/sh
#PBS -e $error
#PBS -o $output
#PBS -l mem=120gb
#PBS -W depend=afterok:${1}
#PBS -l vmem=120gb
#PBS -l pmem=120gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=40:00:00
#PBS -q cmb

cd ${2}
$run_star --outFileNamePrefix $suffix \
      --quantMode GeneCounts \
      --readFilesCommand cat \
      --outSAMtype BAM SortedByCoordinate \
      --bamRemoveDuplicatesType UniqueIdentical \
      --runThreadN 1 \
      --outBAMsortingThreadN 1 \
      --genomeDir $ref_genome \
      --readFilesIn $tempdir/$(cut -d "." -f 1 <<< `basename $r1`)_val_1.fq $tempdir/$(cut -d "." -f 1 <<< `basename $r2`)_val_2.fq
EOS

