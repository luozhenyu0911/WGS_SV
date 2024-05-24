
rule manta_step1_configmanta:
    input:
        bam = "data/{}.bam".format(config['samples']['id']),
        ref = REF 
    output:
        "manta/{id}.manta.vcf.log"
    params:
        config["params"]['manta']
    threads:
        24
    shell:
        """
        {params}/python {params}/configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir manta/ && \
        {params}/python manta/runWorkflow.py -j 24 &> {output} 
        """
