# include the config file
configfile: "config.yaml"

# Define the ref based on the config file
# Sort of acts like a global variable so you don't need to always type the whole thing
REF = config['params']['ref_fa']
SAMPLE = config['samples']['id']

# define a function to return target files based on config settings
def run_all_input(wildcards):

    run_all_files = []
    
    run_all_files.extend(expand("metasv/{id}_{svtype}.vcf", id = config['samples']['id'],svtype=["DEL", "DUP"]))
    run_all_files.append("metasv/{}.metasv.genotype.vcf".format(config['samples']['id']))
    # run_all_files.append("cnvnator/{}.cnvnator.vcf".format(config['samples']['id']))
    # run_all_files.append("manta/{}.manta.vcf.log".format(config['samples']['id']))

    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input

smk_path = config['params']['smk_path']
include: smk_path+"/manta.smk"
include: smk_path+"/metasv2.smk"
include: smk_path+"/lumpy.smk"
include: smk_path+"/annotation.smk"
