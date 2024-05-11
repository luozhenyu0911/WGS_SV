
rule breakdancer_step1:
    input:
        "data/{id}.bam"
    output:
        "breakdancer/{id}.cfg"
    params:
        config['params']['breakdancer']
    shell:
        "{params}/perl {params}/bam2cfg.pl -g -h {input} > {output}"

rule breakdancer_step2:
    input:
        "breakdancer/{id}.cfg"
    output:
        "breakdancer/{id}.cfg.SV.output"
    params:
        config['params']['breakdancer']
    shell:
        "{params}/breakdancer-max {input} > {output}"