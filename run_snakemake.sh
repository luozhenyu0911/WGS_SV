#!/usr/bin/env bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta

smk=/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV/run.all.smk
echo ""
echo "start = $(date)"
echo "$(hostname)"
/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/snakemake -j 24 -pk -s ${smk} 2> snakemake.err.txt
echo "end = $(date)"
# echo "last status $? $(hostname)" |mail -s "job done: $(pwd)" xx.email
