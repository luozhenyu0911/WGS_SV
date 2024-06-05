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


# import os
# if os.path.exists("manta/manta.vcf.log")
#     Manta_vcf = 'manta/results/variants/diploidSV.vcf.gz'

rule metaSV_merge_vcf:
    input:
        manta_vcf = "manta/{id}.manta.sv.vcf",
        lumpyinput = "lumpy/{id}.lumpy.vcf",
        # breakdancerinput= "breakdancer/{id}.cfg.SV.output",
        # cnvnatorinput= "cnvnator/{id}.cnvnator.vcf",
        # pindelinput= "pindel/{id}.pindel.vcf",
        ref = REF
    output:
        "metasv/{id}.metasv.vcf.gz"
    params:
        env = config['params']['metasv'],
        threads = config['threads'],
        sample = config['samples']['id'],
        readlength = config['params']['read_length']
    run:
        if config['params']['depth'] == 7:
            shell(
                "{params.env}/python {params.env}/run_metasv.py "
                    "--reference {input.ref} --sample {params.sample} --disable_assembly --num_threads {params.threads} "
                    "--enable_per_tool_output --keep_standard_contigs --mean_read_length {params.readlength} "
                    "--minsvlen 50 --maxsvlen 10000000 "
                    "--outdir metasv --workdir metasv/tmp_work "
                    "--lumpy_vcf {input.lumpyinput} && "
                "mv metasv/variants.vcf.gz {output}")
        elif config['params']['depth'] == 30:
            shell(
                "{params.env}/python {params.env}/run_metasv.py "
                    "--reference {input.ref} --sample {params.sample} --disable_assembly --num_threads {params.threads} "
                    "--enable_per_tool_output --keep_standard_contigs --mean_read_length {params.readlength} "
                    "--minsvlen 50 --maxsvlen 10000000 "
                    "--outdir metasv --workdir metasv/tmp_work "
                    "--manta_vcf {input.manta_vcf}  && "
                "mv metasv/variants.vcf.gz {output} && "
                "mv metasv/variants.vcf.gz.tbi {output}.tbi")

rule modify_genotype:
    input:
        vcf = "metasv/{id}.metasv.vcf.gz",
        lumpy_vcf = "lumpy/{id}.genotyped.vcf",
        manta_vcf = "manta/{id}.manta.sv.vcf"
    output:
        "metasv/{id}.metasv.genotype.vcf"
    params:
        python = config['params']['python3'],
        src = config['params']['smk_path']
    run:
        if config['params']['depth'] == 7:
            shell("{params.python} {params.src}/src/modify_genotype.py -r {input.lumpy_vcf} -m {input.vcf} -o {output}")
        elif config['params']['depth'] == 30:
            shell("{params.python} {params.src}/src/modify_genotype.py -r {input.manta_vcf} -m {input.vcf} -o {output}")

rule split_vcf_by_svtype:
    input:
        "metasv/{id}.metasv.genotype.vcf".format(id = config['samples']['id'])
    output:
        expand("metasv/{id}_{svtype}.vcf", id = config['samples']['id'], svtype=["DEL", "DUP", "INV", "INS"])
    params:
        sample = config['samples']['id'],
        python = config['params']['python3'],
        src = config['params']['smk_path']
    run:
        shell(
            "{params.python} {params.src}/src/split_vcf_bysvtype.py {input} metasv/{params.sample}")
