rule lumpy_step1_get_discordants_bam:
    input:
        "data/{id}.bam"
    output:
        "{id}/lumpy/{id}.discordants.unsorted.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools view -@ {params.threads} -b -F 1294 {input} > {output}"

rule lumpy_step1_sort_discordants_bam:
    input:
        "{id}/lumpy/{id}.discordants.unsorted.bam"
    output:
        "{id}/lumpy/{id}.discordants.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools sort -@ {params.threads} {input} -o {output}"

rule lumpy_step2_get_splitters_bam:
    input:
        "data/{id}.bam"
    output:
        "{id}/lumpy/{id}.splitters.unsorted.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools view -h {input} | {params.lumpy_path}/extractSplitReads_BwaMem -i stdin | {params.samtools_path}/samtools view -Sb - > {output}"

rule lumpy_step3_sort_splitters_bam:
    input:
        "{id}/lumpy/{id}.splitters.unsorted.bam"
    output:
        "{id}/lumpy/{id}.splitters.bam"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"],
        samtools_path = config["params"]["samtools"]
    shell:
        "{params.samtools_path}/samtools sort -@ {params.threads} {input} -o {output}"

rule lumpy_step4_run_lumpy:
    input:
        discordants = "{id}/lumpy/{id}.discordants.bam",
        splitters = "{id}/lumpy/{id}.splitters.bam",
        raw_bam = "data/{id}.bam"
    output:
        "{id}/lumpy/{id}.lumpy.vcf"
    params:
        lumpy_path = config["params"]["lumpy"],
        threads = config["threads"]
    shell:
        """
        {params.lumpy_path}/lumpyexpress -B {input.raw_bam} -S {input.splitters} -D {input.discordants} -o {output}
        """

rule lumpy_genotyped_vcf:
    input:
        vcf = "{id}/lumpy/{id}.lumpy.vcf",
        bam = "data/{id}.bam"
    output:
        "{id}/lumpy/{id}.genotyped.vcf"
    params:
        env = config["params"]["breakdancer"],
        threads = config["threads"],
        # sample = config["samples"]["id"]
    shell:
        """
        vcftools --vcf {input.vcf} \
        --indv {wildcards.id} --recode --recode-INFO-all  \
        --out {wildcards.id}/lumpy/{wildcards.id} && \
        {params.env}/svtyper \
        -i {wildcards.id}/lumpy/{wildcards.id}.recode.vcf \
        -B {input.bam} \
        -o {output}
        """