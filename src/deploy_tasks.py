import os
import argparse
example_text = '''example: python * '''
parser = argparse.ArgumentParser(description="The script is .",formatter_class= argparse.RawTextHelpFormatter, usage = '%(prog)s [-h]', epilog=example_text)

parser.add_argument('--sample_info','-i',type=str,help=" the sample id in the first cloumn", required= True,metavar='')
parser.add_argument('--config_file', '-c',type= str,help="the template of config.yaml",required= True,metavar='')
parser.add_argument('--snakemakefile', '-s',type= str,help="the template of snakemakefile",required= True,metavar='')
parser.add_argument('--depth', '-d',type= int,help="the depth of sample bam file",required= True,metavar='')
parser.add_argument('--read_length', '-l',type= int,help="the read_length of sample bam file",required= True,metavar='')
parser.add_argument('--output_path', '-p',type= str,help="the output path of the pipeline (like, ./output)",required= True,metavar='')


def file2list(file_path, count=6):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        sample_list = [line.strip().split('\t')[0] for line in lines]
    return [sample_list[i:i + 6] for i in range(0, len(sample_list), 6)]


def modify_config(config_file, runfile ,sample_list, depth, read_length, path):
    os.makedirs('shell', exist_ok=True)
    for i, sample in enumerate(sample_list):
        config_name = 'shell/w'+str(i+1)+'.'+str(depth) +'X.'+ config_file.split('/')[-1]
        with open(config_name, 'w') as outf:
            with open(config_file, 'r') as configf:
                for line in configf:
                    if line.strip(" ").startswith("id"):
                        print('    id_list: '+ str(sample),file = outf)
                    elif line.strip(" ").startswith("depth"):
                        print('    depth: '+str(depth),file = outf)
                    elif line.strip(" ").startswith("read_length"):
                        print('    read_length: '+str(read_length),file = outf)
                    elif line.strip(" ").startswith("PWD"):
                        print('    PWD: '+'"'+str(path)+'"',file = outf)
                    else:
                        outf.write(line)
                        
        runfile_name = 'shell/w'+str(i+1)+'.'+str(depth) +'X.'+ runfile.split('/')[-1]
        
        with open(runfile_name, 'w') as outf:
            with open(runfile, 'r') as runf:
                for line in runf:
                    if line.strip(" ").startswith("tasknum"):
                        id = 'w'+str(i+1)+'.'+str(depth) +'X'
                        print('tasknum='+'"'+id+'"',file = outf)
                    else:
                        outf.write(line)
                        
if __name__ == '__main__':
    args = parser.parse_args()
    sample_list = file2list(args.sample_info)
    modify_config(args.config_file, args.snakemakefile, sample_list, args.depth, args.read_length, args.output_path)