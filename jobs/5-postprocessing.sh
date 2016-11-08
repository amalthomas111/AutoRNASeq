#!/bin/sh
proj=$3
tempdir=$4
outdir=$5
error=$tempdir"/5-postprocessing.error"
output=$tempdir"/5-postprocessing.output"

cat <<EOS | qsub -
#!/bin/sh
#PBS -N postprocessing-${proj}
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
cat $tempdir/*.error >$outdir/readme.log
cp $tempdir/*.final.out $outdir/${proj}.STAR-Mapping-Log
rm -rf $tempdir
EOS

