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
        manta_vcf = "manta/{id}.manta.vcf.log",
        lumpyinput = "lumpy/{id}.lumpy.vcf",
        breakdancerinput= "breakdancer/{id}.cfg.SV.output",
        cnvnatorinput= "cnvnator/{id}.cnvnator.vcf",
        pindelinput= "pindel/{id}.pindel.vcf",
        ref = REF
    output:
        "metasv/{id}.metasv.log"
    params:
        env = config['params']['metasv'],
        threads = config['threads'],
        sample = config['samples']['id'],
        readlength = config['params']['read_length']
    shell:
        """
        {params.env}/python {params.env}/run_metasv.py \
            --reference {input.ref} --sample {params.sample} --disable_assembly --num_threads {params.threads} \
            --enable_per_tool_output --keep_standard_contigs --mean_read_length {params.readlength} \
            --outdir metasv --workdir metasv/tmp_work \
            --breakdancer_native {input.breakdancerinput} \
            --manta_vcf manta/results/variants/diploidSV.vcf.gz \
            --cnvnator_vcf {input.cnvnatorinput} \
            --lumpy_vcf {input.lumpyinput} \
            --pindel_vcf {input.pindelinput} &> {output}
        """

rule metaSV_merge_vcf_no_pindel:
    input:
        manta_vcf = "manta/{id}.manta.vcf.log",
        lumpyinput = "lumpy/{id}.lumpy.vcf",
        breakdancerinput= "breakdancer/{id}.cfg.SV.output",
        cnvnatorinput= "cnvnator/{id}.cnvnator.vcf",
        # pindelinput= "pindel/{id}.pindel.vcf",
        ref = REF
    output:
        "metasv_no_pindel/{id}.SV.vcf.gz"
    params:
        env = config['params']['metasv'],
        threads = config['threads'],
        sample = config['samples']['id'],
        readlength = config['params']['read_length']
    shell:
        # --pindel_vcf {input.pindelinput} \
        """
        {params.env}/python {params.env}/run_metasv.py \
            --reference {input.ref} --sample {params.sample} --disable_assembly --num_threads {params.threads} \
            --enable_per_tool_output --keep_standard_contigs --mean_read_length {params.readlength} \
            --outdir metasv_no_pindel --workdir metasv_no_pindel/tmp_work \
            --overlap_ratio 0.2 --minsvlen 50 --maxsvlen 10000000 \
            --breakdancer_native {input.breakdancerinput} \
            --manta_vcf manta/results/variants/diploidSV.vcf.gz \
            --cnvnator_vcf {input.cnvnatorinput} \
            --lumpy_vcf {input.lumpyinput} && \
        mv metasv_no_pindel/variants.vcf.gz {output} && mv metasv_no_pindel/variants.vcf.gz.tbi {output}.tbi
        """   

rule SV_SeparateFilter:
    input:
        "metasv_no_pindel/{id}.SV.vcf.gz"
    output:
        "metasv_no_pindel/{id}.SV.pass.vcf"
    params:
        src = config['params']['smk_path'],
        sample = config['samples']['id']
    shell:
        """
        python {params.src}/src/VCF_SeparateFilter.py -i {input} -p {params.sample}_pass -o {output}
        """