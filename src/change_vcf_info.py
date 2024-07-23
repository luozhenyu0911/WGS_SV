
import sys
import gzip
import re

def OpenVcf(vcf_file):
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')

def change_info(inpf, outf):
    ENDpos =re.compile(r'\bEND=(-?\d+)\D')
    SVLEN =re.compile(r'\bSVLEN=(-?\d+)\D')
    SVTYPE =re.compile(r'\bSVTYPE=(-?\w+)\b')
    for line in inpf:
        if line.startswith('#'):
            outf.write(line)
        else:
            chrom, pos, id, ref, alt, qual, filter, info, *_ = line.strip().split('\t')
            try:
                end = ENDpos.findall(info)[0]
            except:
                continue
                # end = "None"
            try:
                svlen = SVLEN.findall(info)[0]
                end = ENDpos.findall(info)[0]
            except:
                svlen = int(end) - int(pos)
            svtype = SVTYPE.findall(info)[0]
            info = "END=" + end + ";SVLEN=" + str(svlen) + ";SVTYPE=" + svtype
            print(chrom, pos, id, ref, alt, qual, filter, info,*_, sep="\t", file=outf)
            
if __name__ == '__main__':
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    inpf = OpenVcf(input_vcf)
    with open(output_vcf, 'w') as outf:
        change_info(inpf, outf)
    inpf.close()
