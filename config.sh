#Change the paths to things here
export run_fastqc=`which fastqc`
export run_cutadapt=`which cutadapt`
export run_trim_galore=`which trim_galore`
export run_star=`which STAR`
export run_htseq=`which htseq-count`

# Verbose
echo "[CONFIG] Fastqc path: 	"$run_fastqc
echo "[CONFIG] Cutadapt path: "$run_cutadapt
echo "[CONFIG] Trim galore: 	"$run_trim_galore
echo "[CONFIG] STAR path: 	"$run_star
echo "[CONFIG] HTSeq count: 	"$run_htseq
