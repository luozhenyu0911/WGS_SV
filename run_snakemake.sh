#!/usr/bin/env bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta
pwd=/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/tmp/01.SV/output/00114031182M22BFF2
smk=/BIGDATA2/gzfezx_shhli_2/software/script/WGS_SV/run.all.smk
echo "start = $(date)"
start=$(date)
echo "$(hostname)"
/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/snakemake -d $pwd -j 24 -pk --configfile ${pwd}/config.yaml -s ${smk} 2> ${pwd}/snakemake.err.txt
echo "start = $start"
echo "end = $(date)"
# echo "last status $? $(hostname)" |mail -s "job done: $(pwd)" xx.email
