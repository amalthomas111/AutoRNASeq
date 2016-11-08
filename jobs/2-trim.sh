#!/bin/sh
proj=$3
tempdir=$4
outdir=$5
r1=$6
r2=$7
fqbase=${2}/$(cut -d "." -f 1 <<< `basename $r1`)
error=$tempdir"/2-trimgalore.error"
output=$tempdir"/2-trimgalore.output"

cat <<EOS | qsub -
#!/bin/sh
#PBS -N trim-${proj}
#PBS -S /bin/sh
#PBS -e $error
#PBS -o $output
#PBS -W depend=afterok:${1}
#PBS -l mem=50gb
#PBS -l vmem=50gb
#PBS -l pmem=50gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -q cmb

cd ${2}
if [ -z $r2 ]
then
	$run_trim_galore --path_to_cutadapt $run_cutadapt --output_dir $tempdir $r1
else
	$run_trim_galore --path_to_cutadapt $run_cutadapt --output_dir $tempdir --paired $r1 $r2
fi

EOS
