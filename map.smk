
# Rule for mapping
rule map_reads:
    input:
        ref = REF,
        fq_1 = expand("{sample}", sample=config['samples']['fq_1']),
        fq_2 = expand("{sample}", sample=config['samples']['fq_2'])
    output:
        sam = "Align/{}.aln_mem.sam".format(config['samples']['id'])
    threads:
        config['threads']['hisat2']
    shell:
        # map
        "hisat2 -p {threads} -x {input.ref} -1 {input.fq_1} -2 {input.fq_2} -S {output.sam} 2>Align/aln.err"

# generate raw bam 
rule sam2bam:
    input:
        "Align/{id}.aln_mem.sam"
    output:
        "Align/{id}.raw.bam"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools view -@ {threads} -bhS {input} -o {output}"

# filter reads with mapping quality < INT and unmapped
rule get_uniq_bam:
    input:
        "Align/{id}.raw.bam"
    output:
        "Align/{id}.uniq.bam"
    threads:
        config['threads']['hisat2']
    params:
        reads_map_q = config['params']['reads_map_q']
    shell:
        "samtools view -@ {threads} -bhS -q {params.reads_map_q} -F 0x400 {input} -o {output}"

# sort uniq bam
rule sort_uniq_bam:
    input:
        "Align/{id}.uniq.bam"
    output:
        "Align/{id}.uniq.sort.bam"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools sort -@ {threads} {input} -O bam -o {output} &&"
        "samtools index -@ {threads} {output}"

# samtools flagstat for uniq.bam
rule flagstat_uniq_bam:
    input:
        "Align/{id}.uniq.sort.bam"
    output:
        "Align/{id}.uniq.flagstat.txt"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools flagstat -@ {threads} {input} >> {output}"
        
rule flagstat_raw_bam:
    input:
        "Align/{id}.raw.bam"
    output:
        "Align/{id}.raw.flagstat.txt"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools flagstat -@ {threads} {input} >> {output}"


# get reads count of gene_id or transcript_id 
rule get_reads_count:
    input:
        "Align/{}.uniq.sort.bam".format(config['samples']['id'])
    output:
        gene_count = "Align/{}_gene_count.txt".format(config['samples']['id']),
        transcript_count = "Align/{}_transcript_count.txt".format(config['samples']['id'])
    threads:
        config['threads']['featureCounts']
    params:
        gtf = config['params']['gtf'],
        featureCounts = config['params']['featureCounts']
        
    shell:
        "{params.featureCounts} -T {threads} -p -a {params.gtf} -g gene_name -o {output.gene_count} {input} |"
        "{params.featureCounts} -T {threads} -p -a {params.gtf} -g transcript_id -o {output.transcript_count} {input}"

# This looks at the coverage across the genome, as well as percent coverage at particular depths (4X, 10X, 30X)
rule coverage_depth:
    input:
        "Align/{id}.uniq.sort.bam"
    output:
        "Align/{id}.uniq.bam_coverage_depth.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        ref = config['params']['ref_fa']
    shell:
        "perl {params.toolsdir}/tools/depthV2.0.pl -l $({params.toolsdir}/tools/fasta_non_gapped_bases.py {params.ref}) {input} Align > {output}"
        



















