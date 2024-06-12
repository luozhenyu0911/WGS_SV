#!/usr/bin/env bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
conda activate /BIGDATA2/gzfezx_shhli_2/miniconda3/envs/manta
pwd=/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20240507_SV/mul_test2
smk=/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/WGS_SV_mul/run.all.smk
echo ""
echo "start = $(date)"
echo "$(hostname)"
/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/snakemake -j 24 -pk --configfile ${pwd}/config.yaml -s ${smk} 2> ${pwd}/snakemake.err.txt
echo "end = $(date)"
# echo "last status $? $(hostname)" |mail -s "job done: $(pwd)" xx.email
