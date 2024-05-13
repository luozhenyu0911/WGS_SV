
rule manta_step1_configmanta:
    input:
        bam = "data/{}.bam".format(config['samples']['id']),
        ref = REF 
    output:
        "manta/manta_step1_configmanta.log"
    params:
        config["params"]['manta']
    threads:
        24
    shell:
        """
        {params}/python {params}/configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir manta/ && \
        {params}/python manta/runWorkflow.py -j 24 &> {output} 
        """

rule manta_step2_runmanta:
    input:
        "manta/manta_step1_configmanta.log"
    output:
        "manta/{id}.manta.vcf.log"
    params:
        config["params"]['manta']
    threads:
        24
    shell:
        """
        cp {input} {output}  && {params}/python manta/runWorkflow.py -j 24 
        """
