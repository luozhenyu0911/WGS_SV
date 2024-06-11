
rule breakdancer_step1:
    input:
        "data/{id}.bam"
    output:
        "{id}/breakdancer/{id}.cfg"
    params:
        config['params']['breakdancer']
    shell:
        "{params}/perl {params}/bam2cfg.pl -g -h {input} > {output}"

rule breakdancer_step2:
    input:
        "{id}/breakdancer/{id}.cfg"
    output:
        "{id}/breakdancer/{id}.cfg.SV.output"
    params:
        config['params']['breakdancer']
    shell:
        "{params}/breakdancer-max {input} > {output} && mv {wildcards.id}*.insertsize_histogram* {wildcards.id}/breakdancer/"

