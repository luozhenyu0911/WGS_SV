
rule manta_step1_configmanta:
    input:
        bam = "data/{}.bam".format(config['samples']['id']),
        ref = REF 
    output:
        "manta/{id}.manta.vcf.gz"
    params:
        config["params"]['manta']
    threads:
        24
    shell:
        """
        {params}/python {params}/configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir manta/ && \
        {params}/python manta/runWorkflow.py -j 24 && cp manta/results/variants/diploidSV.vcf.gz {output}
        """
rule manta_step2_filter:
    input:
        "manta/{id}.manta.vcf.gz"
    output:
        "manta/{id}.manta.filtered.vcf"
    params:
        python3 = config["params"]['python3'],
        src = config["params"]['smk_path']
    shell:
        """
        {params.python3} {params.src}/src/manta_filter.py -i {input} -o {output}
        """
