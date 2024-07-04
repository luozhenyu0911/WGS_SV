#coding=utf-8
from __future__ import print_function 
from __future__ import division
from collections import defaultdict
import gzip
import argparse
import re
example_text = '''example:
   python *
   '''
parser = argparse.ArgumentParser(description="The script is to filter SVs with PASS and separete SV types from VCF file.",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--input_vcf','-i',type=str,help="",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="",required= True,metavar='')
parser.add_argument('--output_vcf', '-o',type= str,help="",required= True,metavar='')
args = parser.parse_args()

input_vcf = args.input_vcf
prefix = args.prefix
output_vcf = args.output_vcf

# prefix = output.split(".")[0]
SV_len =re.compile(r';SVLEN=(-?\d+?);')
annotation = []
SVtypes_dict = defaultdict(list)
SVtypes_stat = defaultdict(dict)
sv_len_dict = defaultdict(list)

annotation = []
with gzip.open(input_vcf, 'rt') as f:
    for line in f:
        fields = line.strip().split('\t')
        if line.startswith('#'):
            annotation.append(line)
annotation.insert(2, '##INFO=<ID=natorP1,Number=.,Type=Integer,Description="maybe its related to the CNVnator">\n')
annotation.insert(2, '##INFO=<ID=natorP2,Number=.,Type=Integer,Description="maybe its related to the CNVnator">\n')
annotation.insert(2, '##INFO=<ID=natorP3,Number=.,Type=Integer,Description="maybe its related to the CNVnator">\n')
annotation.insert(2, '##INFO=<ID=natorP4,Number=.,Type=Integer,Description="maybe its related to the CNVnator">\n')
annotation.insert(2, '##INFO=<ID=natorQ0,Number=.,Type=Integer,Description="maybe its related to the CNVnator">\n')
annotation.insert(2, '##INFO=<ID=natorRD,Number=.,Type=Integer,Description="maybe its related to the CNVnator">\n')

with open(output_vcf, 'w') as outf:
    outf.writelines(annotation)
    with gzip.open(input_vcf, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if line.startswith('#'):
                pass
            elif "PASS" in fields[6]:
                sv_len = SV_len.findall(line.strip())[0]
                outf.write(line)
                SVtypes_dict[fields[4]].append(line)
                SVtypes_stat[fields[4]] = SVtypes_stat.get(fields[4], {'count': 0, 'len': 0})
                SVtypes_stat[fields[4]]['count'] += 1
                SVtypes_stat[fields[4]]['len'] += abs(int(sv_len))
                sv_len_dict[fields[4]].append(abs(int(sv_len)))
                

for SVtype in SVtypes_dict:
    with open(prefix + "_" + SVtype.strip("<").strip(">") + ".vcf", 'w') as outf:
        outf.writelines(annotation)
        outf.writelines(SVtypes_dict[SVtype])

with open(prefix + '_svtypr_stats.txt', 'w') as f:
    print("SVtype", "Count", "Length", sep="\t", file=f)
    for svtype in SVtypes_stat:
        print(svtype.strip(">").strip("<"), SVtypes_stat[svtype]['count'], SVtypes_stat[svtype]['len'], file=f, sep="\t")
        
with open(prefix + '_svtype_len.txt', 'w') as f:
    for i in sv_len_dict:
        for j in sv_len_dict[i]:
            print(i.strip(">").strip("<"), j, prefix, file=f, sep="\t")
            