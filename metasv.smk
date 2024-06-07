from pathlib import Path

def get_cram(wildcards):
    fq_path = Path(config['samples']['cram_path'])
    name = config['samples']['id']
    cram_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches
        if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".cram"):
            print(fp)
            cram_files.append(str(fp))
    return cram_files 

def get_cram_index(wildcards):
    fq_path = Path(config['samples']['cram_path'])
    name = config['samples']['id']
    cram_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches 
        if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".crai"):
            print(fp)
            cram_files.append(str(fp))
    return cram_files

# link the cram and cram index files to the data directory
rule link_cram:
    input:
        get_cram
    output:
        "data/{id}.cram"
    run:
        shell("ln -s {input} {output}")

rule link_cram_index:
    input:
        get_cram_index
    output:
        "data/{id}.cram.crai"
    run:
        shell("ln -s {input} {output}")

# convert the cram file to bam file
rule bam_from_cram:
    input:
        "data/{id}.cram"
    output:
        "data/{id}.bam"
    threads:
        config['threads']
    params:
        samtools_path = config['params']['samtools']
    shell:
        "{params.samtools_path}/samtools view -@ {threads} -bh {input} -o {output} && {params.samtools_path}/samtools index -@ {threads} {output}"


rule metaSV_merge_vcf:
    input:
        manta_vcf = "manta/{id}.manta.vcf.gz",
        lumpyinput = "lumpy/{id}.genotyped.vcf",
        breakdancerinput= "breakdancer/{id}.cfg.SV.output",
        cnvnatorinput= "cnvnator/{id}.cnvnator.vcf",
        ref = REF
    output:
        metasv_vcf = "metasv/{id}.SV.vcf.gz"
        # breakdancer = "metasv/breakdancer.vcf.gz",
        # manta = "metasv/manta.vcf.gz",
        # cnvnator = "metasv/cnvnator.vcf.gz",
        # lumpy = "metasv/lumpy.vcf.gz"
    params:
        env = config['params']['metasv'],
        threads = config['threads'],
        sample = config['samples']['id'],
        readlength = config['params']['read_length']
    run:
        shell(
        "{params.env}/python {params.env}/run_metasv.py "
            "--reference {input.ref} --sample {params.sample} --disable_assembly --num_threads {params.threads} "
            "--enable_per_tool_output --keep_standard_contigs --mean_read_length {params.readlength} "
            "--outdir metasv --workdir metasv/tmp_work "
            "--overlap_ratio 0.5 --minsvlen 50 --maxsvlen 10000000 "
            "--breakdancer_native {input.breakdancerinput} "
            "--manta_vcf {input.manta_vcf} "
            "--cnvnator_vcf {input.cnvnatorinput} "
            "--lumpy_vcf {input.lumpyinput} && "
        "mv metasv/variants.vcf.gz {output.metasv_vcf} && mv metasv/variants.vcf.gz.tbi {output.metasv_vcf}.tbi")

rule add_genotype_to_metasv_lumpy:
    input:
        lumpy_vcf = "lumpy/{id}.genotyped.vcf",
        # modify_vcf = "metasv/lumpy.vcf.gz",
        metasv_vcf = "metasv/{id}.SV.vcf.gz"
    output:
        "metasv/{id}.lumpy.gt.vcf.gz"
    params:
        python = config['params']['python3'],
        src = config['params']['smk_path']
    shell:
        "{params.python} {params.src}/src/modify_genotype.py -r {input.lumpy_vcf} -m metasv/lumpy.vcf.gz -o {output}"

rule add_genotype_to_metasv_manta:
    input:
        manta_vcf = "manta/{id}.manta.vcf.gz",
        # modify_vcf = "metasv/manta.vcf.gz",
        metasv_vcf = "metasv/{id}.SV.vcf.gz"
    output:
        "metasv/{id}.manta.gt.vcf.gz"
    params:
        python = config['params']['python3'],
        src = config['params']['smk_path']
    shell:
        "{params.python} {params.src}/src/modify_genotype.py -r {input.manta_vcf} -m metasv/manta.vcf.gz -o {output}"

rule done:
    input:
        "metasv/{id}.SV.vcf.gz",
        "metasv/{id}.lumpy.gt.vcf.gz",
        "metasv/{id}.manta.gt.vcf.gz"
    output:
        "{id}_have_done.txt"
    shell:
        "touch {output}"
# rule split_vcf_by_svtype:
#     input:
#         "metasv/{id}.metasv.genotype.vcf".format(id = config['samples']['id'])
#     output:
#         expand("metasv/{id}_{svtype}.vcf", id = config['samples']['id'], svtype=["DEL", "DUP", "INV", "INS"])
#     params:
#         sample = config['samples']['id'],
#         python = config['params']['python3'],
#         src = config['params']['smk_path']
#     run:
#         shell(
#             "{params.python} {params.src}/src/split_vcf_bysvtype.py {input} metasv/{params.sample}")
