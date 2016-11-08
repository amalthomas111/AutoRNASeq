#!/bin/sh
proj=$2
readdir=$3
tempdir=$4
outdir=$5
r1=$6
r2=$7
error=$tempdir"/1-fastqc.error"
output=$tempdir"/1-fastqc.output"

cat <<EOS | qsub -
#!/bin/sh
#PBS -N fqc-${proj}
#PBS -S /bin/sh
#PBS -e $error
#PBS -o $output
#PBS -l mem=50gb
#PBS -l vmem=50gb
#PBS -l pmem=50gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -q cmb

cd ${1}
$run_fastqc $r1 $r2
mv ${readdir}/*.html ${readdir}/*.zip $outdir

EOS

