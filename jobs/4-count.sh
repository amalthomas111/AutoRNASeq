#!/bin/sh
proj=$3
tempdir=$4
outdir=$5
gtf=$6
error=$tempdir"/4-count.error"
output=$tempdir"/4-count.output"

cat <<EOS | qsub -
#!/bin/sh
#PBS -N count-${proj}
#PBS -S /bin/sh
#PBS -e $error
#PBS -o $output
#PBS -l mem=120gb
#PBS -W depend=afterok:${1}
#PBS -l vmem=30gb
#PBS -l pmem=30gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -q cmb

cd $2
samtools view $tempdir/${proj}.Aligned.sortedByCoord.out.bam | $run_htseq --stranded=no - $gtf >$tempdir/${proj}.counts
cut -f2 $tempdir/${proj}.counts | head -n -5 > $outdir/${proj}.vector

EOS
