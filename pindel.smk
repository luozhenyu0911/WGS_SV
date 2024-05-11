pindel -i pindel.conf -f $REF -o $sample
# pindel2vcf -r <参考基因组文件> -R <参考基因组名称> -d 参考基因组日期 -p <pindel输出文件> -e <最小的reads数>
pindel2vcf -r $REF -R hg38 -d 20140111 -P $sample -v $sample.vcf

rule 

rule pindel_step1:
    input:
        bam = "data/{id}.bam",
        ref = config["REF"]
    output:
        "pindel/{id}.vcf"

    params:
        config["pindel"]

    shell:
        "{params}/"