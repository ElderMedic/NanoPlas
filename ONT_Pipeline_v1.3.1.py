#!/usr/bin/env python
# author: Changlin Ke, NGS BI, Macrogen Europe
# date: 2023-06-02

import os, re, sys, time
import subprocess
import pandas as pd
import numpy as np
import argparse
import json
import glob
from datetime import datetime, timedelta
from tqdm import tqdm
import pickle
from Bio import SeqIO

from zip_results import zip_results

'''
main function of pipeline

work to do:
1. add dorado basecalling/barcoding support
2. fix fastq file name matching in mergefq
3. basecalling qc notebook imported

'''
# Work list


# find corresponding barcode id for each input barcode name
def get_selected_barcode(sample_config)->"pd.DataFrame":
    # take the last two digits at the end of Barcode assigned col in sample config table
    sample_config['barcode_guppy'] = 'barcode' + sample_config['Barcode assigned'].str.extract('(\d+)', expand=False).str.zfill(2)
    return sample_config

def parse_args()->"args":
    """
        Parses inputs from commandline and returns them as a Namespace object
        :return: inputs as a Namespace object
    """
    parser = argparse.ArgumentParser(description='''

███╗   ██╗ █████╗ ███╗   ██╗ ██████╗ ██████╗ ██╗      █████╗ ███████╗    
████╗  ██║██╔══██╗████╗  ██║██╔═══██╗██╔══██╗██║     ██╔══██╗██╔════╝    
██╔██╗ ██║███████║██╔██╗ ██║██║   ██║██████╔╝██║     ███████║███████╗    
██║╚██╗██║██╔══██║██║╚██╗██║██║   ██║██╔═══╝ ██║     ██╔══██║╚════██║    
██║ ╚████║██║  ██║██║ ╚████║╚██████╔╝██║     ███████╗██║  ██║███████║    
╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝╚══════╝    

#########################################################################
## MacrogenEU Nanapore Plasmid Sequencing Pipeline                     
## NGS-BI                                
## ver 1.3.1                                                          
## Workflow:  
## 1. Guppy basecalling and demultiplexing
## 2. QC of basecalled fastq files
## 3. Alignment to consensus assembly, plannotate
## 4. QC of alignment bam files
## 5. Variant calling
## 6. QC of variant calling vcf files
## Planed: 
## (add min-overlap manual setting/calculate by fraction to flye
## 7. Annotation of vcf files
## 8. QC of annotation files
## 9. Report generation
##
## Usage: at base environment (to be updated!), run with python 3.10
## Output: at working directory, a folder named with the sample ID will be created in /Results
#########################################################################
                                                                         
                                                                         ''', formatter_class=argparse.RawTextHelpFormatter)
    # parser.add_argument('-i', '--input', help='Input folder', required=False)
    parser.add_argument('-wd', '--working-dir', help='Output folder and working directory to find config files', required=True)
    parser.add_argument('-c', '--config', help='Config file name', required=False, default='pipeline_config.json')
    parser.add_argument('-s', '--sample-config', help='Sample config file name', required=False, default='sample_config.tsv')
    parser.add_argument('--no-guppy', help='skip guppy basecalling and barcoding', required=False, action='store_true', default=False)
    # don't show any cmd running info
    parser.add_argument('--quiet', help='No cmd running info in debug', required=False, action='store_true', default=False)

    # skip a certain step, feed a list of steps to skip
    # keep this choice list updated
    parser.add_argument('--skip', help='skip a certain step, select in a list', required=False, choices=['mergefq', 'assembly', 'plannotate', 'minimap2', 'samtools', 'bcftools', 'variant_call','copy'], nargs='+', default=[])

    # parser.add_argument('-n', '--name', help='Name', required=True)
    # parser.add_argument('-e', '--email', help='Email', required=True)

    args = parser.parse_args()
    return args

# subprocess run mode: define a shell cmd runnning function in python environment
def run_command(command, shell_use=False)->"ListOutput":
    start_time = time.time()
    start_datetime = datetime.now()
    print(f"Start date and time: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stdout)

    # running the cmd
    if not quiet:
        if shell_use:
            print("Running command: ",command, file=sys.stdout)
        else:
            print("Running command: "," ".join(command), file=sys.stdout)
    process = subprocess.run(command, check=False, capture_output=True, shell=shell_use)

    end_time = time.time()
    end_datetime = datetime.now()
    print(f"End date and time: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stdout)

    elapsed_time = end_time - start_time
    elapsed_hms = str(timedelta(seconds=elapsed_time))
    print(f"Elapsed time: {elapsed_hms}\n", file=sys.stdout)
    return process.stdout, process.stderr, process.returncode

# prepare pipeline_config and sample_config file
def get_config(base_dir, **kwargs)->"dict":
    path_sample_config = os.path.join(base_dir, sample_config_file)
    sample_config = pd.read_csv(path_sample_config, sep="\t")
    print('Reading Sample config file: ', path_sample_config, file=sys.stdout)    
    sample_config = get_selected_barcode(sample_config) # only exceute once

    # a pipeline_config.json file is required at base_dir for the pipeline
    with open(os.path.join(base_dir, config_file)) as f:
        pipeline_config = json.load(f)
    print('Reading pipeline config file: ', os.path.join(base_dir, config_file), file=sys.stdout)
    pipeline_config["sample_config"] = sample_config

    return pipeline_config

# show output at error
def show_output_atError(*process):
    for p in process:
        if p is not None and p[2] != 0:
            print("===== Debugging info: =====\n")
            print("StdErr: \n", p[0].decode("utf-8"),"\n")
            print("StdOut: \n", p[1].decode("utf-8"),"\n")


'''
Running each step/tool of the pipeline
'''
# define an expected output object/format of functions "ListOutput"
# ListOutput = [stdout, stderr, returncode]

# Guppy basecaller and demultiplexing
def run_guppy_basecaller(config)->"ListOutput":
    guppy_basecalling_command = [
    'guppy_basecaller',
    '--disable_pings',
    '-i', config['input'],
    '-s', config['output'],
    '-c', config['config'],
    '-x', 'auto', # auto detect gpu
    '-r',
    '--num_callers','8',
    '--gpu_runners_per_device','10',
    '--detect_adapter',
    '--detect_mid_strand_adapter',
    '--do_read_splitting',
    '--trim_adapters',
    '--min_qscore', config['min_qscore'],
    '--compress_fastq'
]   
    # print(" ".join(guppy_basecalling_command))
    process_basecalling = list(run_command(guppy_basecalling_command))
    return process_basecalling

def run_guppy_barcoding(config)->"ListOutput":
    guppy_barcoding_command = [
    'guppy_barcoder',
    '--disable_pings',
    '-i', config['input'],
    '-s', config['output'],
    '-c', config['config'],
    '-x', 'auto',
    '-r',
    '--detect_barcodes',
    '--detect_mid_strand_barcodes',
    '--barcode_kits', config['barcode_kits'],
    '--enable_trim_barcodes',
    '--fastq_out','--compress_fastq'
]
    # print(" ".join(guppy_barcoding_command))
    process_barcoding = list(run_command(guppy_barcoding_command))
    return process_barcoding


def basecall_qc(sub_dir:str, fastq_file:str, pipeline_config)->"ListOutput":
    # Concatenate fastq.gz files for all runs in the same barcode subdirectory
    files_to_merge_path = os.path.join(sub_dir, pipeline_config['mergefq']['files_to_merge'])
    command_mergefq = f'cat {files_to_merge_path} > {fastq_file}'

    command_nanoplot = [
        'NanoPlot',
        '--fastq',f'{fastq_file}', 
        '-o',f'{sub_dir}',
        '--tsv_stats', 
        #'--drop_outliers',
    ]

    process_mergefq = list(run_command(command_mergefq, shell_use=True)) #single string command + shell=True, since no wildcards allowed in list
    process_nanoplot = list(run_command(command_nanoplot))
                    
    return process_mergefq, process_nanoplot
    # process_fastqc


# Assembly and polish, annotate, copy to results folder
def assembly(sub_dir, fastq_file, dir_name, pipeline_config:dict)->"ListOutput":
    sample_config = pipeline_config['sample_config'] # sample config table
    plasmid_size = sample_config[sample_config['barcode_guppy']==dir_name].loc[:,'Plasmid size (bp)'].values[0]
    # selected_barcode = sample_config['Barcode assigned'].unique().tolist() # existing barcodes in table 
    command_flye = [
        'flye',
        pipeline_config['flye']['data_type'], f'{fastq_file}',
        '--genome-size', str(plasmid_size), # take only the first one
        '--iterations', pipeline_config['flye']['iteration'],
        '-o',f'{sub_dir}',
        '--threads', pipeline_config['flye']['threads'],
    ]

    # add optional parameters
    if pipeline_config['flye']['asm_coverage']:
        command_flye.extend(['--asm-coverage', str(sample_config[sample_config['barcode_guppy']==dir_name].loc[:,'asm_coverage'].values[0])]) # take only the first one
    if pipeline_config['flye']['meta']:
        command_flye.append('--meta')
    if pipeline_config['flye']['no-alt']:
        command_flye.append('--no-alt-contigs')

    # to add: gpu with medaka
    command_medaka = [
        'medaka_consensus',
        '-f',
        '-x',
        '-i',f'{fastq_file}', 
        '-o',f'{sub_dir}',
        '-d',f'{os.path.join(sub_dir, "assembly.fasta")}',
        '-m', pipeline_config['medaka']['model'],
        '-t', pipeline_config['medaka']['threads'],
    ]
    

    # print('Executed cmd: ', " ".join(command_flye))
    # print('Executed cmd: ', " ".join(command_medaka))
    process_flye = list(run_command(command_flye))
    process_medaka = list(run_command(command_medaka))
    return process_flye, process_medaka
 
# detect if there's only one read per fasta file otherwise split the fasta file by read, name it as consensus_[read_name].fasta
# run plannotate on each fasta file
# merge all plannotate results into one html file
def run_plannotate(dir_name, sub_dir, pipeline_config):
    sample_config = pipeline_config['sample_config'] # sample config table
    plasmid_size = sample_config[sample_config['barcode_guppy']==dir_name].loc[:,'Plasmid size (bp)'].values[0]

    def detect_single_read(fasta_file):
        # Parse the FASTA file and count the sequences
        records = list(SeqIO.parse(fasta_file, 'fasta'))
        num_sequences = len(records)

        if num_sequences == 1:
            return True
        else:
            return False
        
    def split_fasta_by_read(fasta_file, sub_dir, plasmid_size):
        # Parse the FASTA file and create separate files for each read
        records = SeqIO.parse(fasta_file, 'fasta')
        closest_read = None
        closest_length_diff = float('inf')
        list_split = []

        for record in records:
            read_name = record.id
            read_length = len(record.seq)
            filename = os.path.join(sub_dir, f"consensus_{read_name}.fasta")
            list_split.append(filename)
            SeqIO.write(record, filename, 'fasta')
            print(f"Read: {read_name} | Length: {read_length}")

            # select the read with length closest to the plasmid size
            length_diff = abs(read_length - plasmid_size)
            if length_diff < closest_length_diff:
                closest_length_diff = length_diff
                closest_read = read_name
        print(f"Closest to target size read: {closest_read} | Length difference: {closest_length_diff}")

        return list_split, closest_read

    def run_plannotate_on_fasta(consensus_file, sub_dir, pipline_config):
        command_plannotate = [
            "plannotate", "batch",
            '-i', consensus_file, # sub_dir+consensus.fasta
            '-o', sub_dir, 
            '--html',
            ]
        custom_db = pipeline_config['plannotate']['customdb']
        if custom_db is not None:
            command_plannotate.extend(['-y', custom_db])
        process_plannotate = list(run_command(command_plannotate))
        return process_plannotate
    
    # detect if there's consensus.fasta
    if not os.path.exists(os.path.join(sub_dir, "consensus.fasta")):
        print('No consensus.fasta found in ', sub_dir)
        return None
    else:    
        consensus_file = os.path.join(sub_dir, 'consensus.fasta')
        if detect_single_read(consensus_file):
            print("The consensus file contains one single read.")
            process_plannotate = run_plannotate_on_fasta(consensus_file, sub_dir, pipeline_config)
        else:
            print("The consensus file contains multiple reads.\n", file=sys.stderr)
            list_split_consensus_file, closest_read = split_fasta_by_read(consensus_file, sub_dir, plasmid_size) # possible to use closest read only
            for split_consensus_file in list_split_consensus_file:
                print(f"Split the consensus FASTA file into {split_consensus_file}.\n", file=sys.stderr)
                # if split_consensus_file.endswith(f"consensus_{closest_read}.fasta"):
                process_plannotate = run_plannotate_on_fasta(split_consensus_file, sub_dir, pipeline_config)
                # else:
                    # print(f"Skip {split_consensus_file}.\n", file=sys.stderr)
    
        return process_plannotate



# run minimap2 for all barcodes and skip if there's no consensus.fasta
# example: minimap2 -ax map-ont "/media/macrogen-ont/ONT_DATA/ONT_example_data_23052023/Barcoding/barcode25/consensus.fasta" "/media/macrogen-ont/ONT_DATA/ONT_example_data_23052023/Barcoding/barcode25/barcode25.fastq.gz" > /media/macrogen-ont/ONT_DATA/ONT_example_data_23052023/Results/SP307/barcode25.sam
def run_minimap2(sub_dir, fastq_file, dir_name, pipeline_config)->"ListOutput":
    # detect if there's consensus.fasta
    if not os.path.exists(os.path.join(sub_dir, "consensus.fasta")):
        print('No consensus.fasta found in ', sub_dir)
        return None
    else:    
        command_minimap2 = [
            'minimap2',
            '-ax', 'map-ont',
            f'{os.path.join(sub_dir, "consensus.fasta")}',
            f'{fastq_file}',
            '-o', f'{os.path.join(sub_dir, dir_name+".sam")}',
        ]
        process_minimap2 = list(run_command(command_minimap2))
        return process_minimap2
    

# To convert SAM from Minimap2 resutls to BAM file:
# samtools view -bS alignment.sam > alignment.bam
# samtools sort <input.bam> -o <output.bam>
# samtools index <input.bam>
# samtools rmdup <input.bam> <output.bam>
def ConcertBAM(sub_dir, dir_name)->"ListOutput":
    # detect if there's SAM file
    if not os.path.exists(os.path.join(sub_dir, dir_name+".sam")):
        print('No sam file found in ', sub_dir)
        return None, None, None, None
    else:
        command_samtools_view = [
            'samtools',
            'view',
            '-bS',
            f'{os.path.join(sub_dir, dir_name+".sam")}',
            '-o', f'{os.path.join(sub_dir, dir_name+".bam")}',
        ]

        command_samtools_sort = [
            'samtools',
            'sort',
            f'{os.path.join(sub_dir, dir_name+".bam")}',
            '-o', f'{os.path.join(sub_dir, dir_name+".sorted.bam")}',
        ]

        command_samtools_rmdup = [
            'samtools',
            'rmdup',
            f'{os.path.join(sub_dir, dir_name+".sorted.bam")}',
            f'{os.path.join(sub_dir, dir_name+".sorted.rmdup.bam")}',
        ]

        command_stamtools_index = [
            'samtools',
            'index',
            f'{os.path.join(sub_dir, dir_name+".sorted.rmdup.bam")}',
        ]

        process_samtools_view = list(run_command(command_samtools_view))
        process_samtools_sort = list(run_command(command_samtools_sort))
        process_samtools_rmdup = list(run_command(command_samtools_rmdup))
        process_stamtools_index = list(run_command(command_stamtools_index))

        return process_samtools_view, process_samtools_sort, process_samtools_rmdup, process_stamtools_index
    
    
# bcftools mpileup -Ob -o <study.bcf> -f <ref.fa> <sample1.bam> <sample2.bam> <sample3.bam>
# bcftools call -vmO z -o <study.vcf.gz> <study.bcf>
# bcftools view -e 'ALT="."' -Oz -o test_barcode02.vcf
# tabix -p vcf <study.vcf.gz>
def run_bcftools(sub_dir, dir_name, pipeline_config)->"ListOutput":
    # detect if there's BAM file
    if not os.path.exists(os.path.join(sub_dir, dir_name+".sorted.rmdup.bam")):
        print('No BAM file found in ', sub_dir)
        return None, None, None, None
    else:
        command_bcftools_mpileup = [
            'bcftools',
            'mpileup',
            '-Ou', '-B',
            '-d', '10000', # max per-file depth
            '-q', '0', # min mapping Q
            '-Q', pipeline_config['guppy_basecaller']['min_qscore'], # min per-sample BQ
            '-o', f'{os.path.join(sub_dir, dir_name+"_bcftools.bcf")}',
            '-f', f'{os.path.join(sub_dir, "consensus.fasta")}',
            f'{os.path.join(sub_dir, dir_name+".sorted.rmdup.bam")}',
            '--threads', pipeline_config['bcftools']['threads'],
        ]

        command_bcftools_call = [
            'bcftools',
            'call',
            '-Ov',
            '-c', #consensus-caller
            '-A', 
            '--ploidy', '1',
            '-o', f'{os.path.join(sub_dir, dir_name+"_bcftools.vcf")}',
            f'{os.path.join(sub_dir, dir_name+"_bcftools.bcf")}',
            '--threads', pipeline_config['bcftools']['threads'],
        ]

        command_view = [
            'bcftools',
            'view',
            '-e','ALT="."', # exclude sites with no ALTs
            '-Oz',
            '-o', f'{os.path.join(sub_dir, dir_name+"_bcftools.vcf.gz")}',
            f'{os.path.join(sub_dir, dir_name+"_bcftools.vcf")}',
        ]

        command_tabix = [
            'tabix',
            '-p', 'vcf',
            '-f',
            f'{os.path.join(sub_dir, dir_name+"_bcftools.vcf.gz")}',
        ]

        process_bcftools_mpileup = list(run_command(command_bcftools_mpileup))
        process_bcftools_call = list(run_command(command_bcftools_call))
        process_view = list(run_command(command_view))
        process_tabix = list(run_command(command_tabix))

        return process_bcftools_mpileup, process_bcftools_call, process_view, process_tabix
    

# bacterial ploidy-1 variant calling with Medaka
def medaka_variant_call(sub_dir, fastq_file, pipeline_config)->"ListOutput":
    # detect if there's consensus.fasta
    if not os.path.exists(os.path.join(sub_dir, "consensus.fasta")):
        print('No consensus.fasta found in ', sub_dir)
        return None
    else:
        command_medaka_variant = [
            'medaka_haploid_variant',
            '-r', f'{os.path.join(sub_dir, "consensus.fasta")}',
            '-i', f'{fastq_file}',
            '-o', f'{sub_dir}',
            '-m', pipeline_config['medaka']['variant_model'],
            '-t', pipeline_config['medaka']['threads'],
        ]

        process_medaka_variant = list(run_command(command_medaka_variant))

        return process_medaka_variant
    

# SV calling with sniffles
# sniffles --input mapped_input.bam --vcf output.vcf --non-germline
def sniffles_variant_call(sub_dir, dir_name, pipeline_config)->"ListOutput":
    # detect if there's BAM file
    if not os.path.exists(os.path.join(sub_dir, dir_name+".sorted.rmdup.bam")):
        print('No BAM file found in ', sub_dir)
        return None
    else:
        command_sniffles_variant = [
            'sniffles',
            '--input', f'{os.path.join(sub_dir, dir_name+".sorted.rmdup.bam")}',
            '--vcf', f'{os.path.join(sub_dir, dir_name+"_sniffles.vcf")}',
            '--non-germline',
            '--allow-overwrite',
        ]

        process_sniffles_variant = list(run_command(command_sniffles_variant))

        return process_sniffles_variant

# convert all vcf files in sub_dir to csv files using parse_vcf_to_csv.sh
# sh /media/macrogen-ont/ONT_DATA/parse_vcf_to_csv.sh [input.vcf] [output.csv]   
def parse_vcf_csv(sub_dir, dir_name, pipeline_config)->"ListOutput":
    # detect if there's any vcf file in sub_dir regardless of the name
    if not glob.glob(os.path.join(sub_dir, "*.vcf")):
        print('No VCF file found in ', sub_dir)
        return None        
    else:
        # store all found vcf files in a list
        vcf_files = []
        for file in os.listdir(sub_dir):
            if file.endswith(".vcf"):
                vcf_files.append(file)
        
        # run comdand of parse_vcf_to_csv.sh for each vcf file
        for vcf_file in vcf_files:
            command_parse_vcf_csv = [
                'sh',
                f'{pipeline_config["parse_vcf_csv"]["script"]}',
                f'{os.path.join(sub_dir, vcf_file)}',
                f'{os.path.join(sub_dir, vcf_file.replace(".vcf", "_vcf.csv"))}', # replace suffix .vcf with _vcf.csv
            ]
            
            process_parse_vcf_csv = list(run_command(command_parse_vcf_csv))

        return process_parse_vcf_csv


# copy files list in pipeline_config['files_to_copy'] to f'{base_dir}/Results' folder in one go
def copy_results(base_dir, fastq_file, sub_dir, dir_name, pipeline_config)->"ListOutput":
    files_to_copy = pipeline_config['files_to_copy']
    sampleid = pipeline_config['sample_config'][pipeline_config['sample_config']['barcode_guppy']==dir_name].loc[:,'Sample name'].values[0]
    print('Copying files to Sample ID: ', sampleid)

    # remove spaces in sampleid
    sampleid = sampleid.replace(' ','_')
    dest_dir = os.path.join(base_dir, 'Results', sampleid) # base_dir/Results/RBK83/
    os.makedirs(dest_dir, exist_ok=True)

    # detect if there's any file in sub_dir matching file names in files_to_copy
    list_files_to_copy = []
    for i in files_to_copy:
        list_files_to_copy.extend(glob.glob(os.path.join(sub_dir, i)))
    if not list_files_to_copy:
        print('No files found in barcoding folder: ', sub_dir)
    
    # copy files to dest_dir
    copy_files = 'cp'+ ' '+ " ".join(list_files_to_copy) + ' ' + f'{dest_dir}' # this command will fail if no file found
    
    # detect if all files in files_to_copy are in dest_dir, output a list of files not found
    files_not_found = [i for i in files_to_copy if not glob.glob(os.path.join(dest_dir, i))]
    if files_not_found:
        print('Files not found: ', files_not_found, file=sys.stdout)
    else:
        print('All files copied to ', dest_dir, file=sys.stdout)

    # print('Executed cmd: ', copy_files)
    process_copy = list(run_command(copy_files, shell_use=True))
    return process_copy, files_not_found, sampleid

    


'''
Main function
workflow control

'''
def main(base_dir, **kwargs): 
    # skip_variant_call=False, skip_copy=False
    # get config file
    pipeline_config = get_config(base_dir, sample_config_file = sample_config_file, config_file = config_file)

    # create subfolders
    for major_dir in [pipeline_config['guppy_basecaller']['output'],
                      pipeline_config['guppy_barcoder']['output'], 
                      os.path.join(base_dir, "Results"),
                      ]:
        os.makedirs(major_dir, exist_ok=True)
    
    # run guppy for all barcodes
    # about 30min for 27 barcodes
    if not no_guppy:
        process_basecalling = run_guppy_basecaller(pipeline_config['guppy_basecaller'])
        show_output_atError(process_basecalling)
        process_barcoding = run_guppy_barcoding(pipeline_config['guppy_barcoder'])
        show_output_atError(process_barcoding)
    else:
        print("Skip guppy basecalling and barcoding! \n", file=sys.stdout)
    
    # if skip_list is not empty list, print out
    if skip_list:
        print("Skip the following steps: ", skip_list, '\n', file=sys.stdout)
    
    # run pipeline for each barcode
    pattern = re.compile(r'barcode\d{2}') # pattern to match barcodeXX, two digits at the end
    dict_process = {}
    dict_results_files = {}
    selected_barcode = pipeline_config['sample_config']['barcode_guppy'].unique().tolist() # existing barcodes in table 

    # for each barcode, run pipeline:
    for root, dirs, _ in os.walk(f'{base_dir}/Barcoding'):
        filtered_dirs = [i for i in dirs if pattern.match(i) and i in selected_barcode] # if 'barcodeXX' selected in sample config in subfolder name
        for dir_name in tqdm(filtered_dirs):
    # for dir_name in dirs:
        # if pattern.match(dir_name): # if 'barcodeXX' in subfolder name
        #     if dir_name in selected_barcode: # if barcode is selected
            sub_dir = os.path.join(root, dir_name)
            fastq_file = os.path.join(sub_dir, f'{dir_name}.fastq.gz')
            
            # save all processes in a dictionary if the process is not skipped
            dict_process[dir_name] = {}
            # main workflow
            
            if 'mergefq' not in skip_list:
                process_mergefq, process_nanoplot = basecall_qc(sub_dir, fastq_file, pipeline_config)
                show_output_atError(process_mergefq, process_nanoplot)
                dict_process[dir_name]['mergefq'] = process_mergefq
                dict_process[dir_name]['nanoplot'] = process_nanoplot

            if 'assembly' not in skip_list:
                process_flye, process_medaka = assembly(sub_dir, fastq_file, dir_name, pipeline_config)
                show_output_atError(process_flye, process_medaka)
                dict_process[dir_name]['flye'] = process_flye
                dict_process[dir_name]['medaka'] = process_medaka
            
            if 'plannotate' not in skip_list:
                process_plannotate = run_plannotate(dir_name, sub_dir, pipeline_config)
                show_output_atError(process_plannotate) # no output for plannotate
                dict_process[dir_name]['plannotate'] = process_plannotate


            if 'minimap2' not in skip_list:
                process_minimap2 = run_minimap2(sub_dir, fastq_file, dir_name, pipeline_config)
                show_output_atError(process_minimap2)
                dict_process[dir_name]['minimap2'] = process_minimap2

            if 'samtools' not in skip_list:
                process_samtools_view, process_samtools_sort, process_stamtools_index, process_samtools_rmdup = ConcertBAM(sub_dir, dir_name)
                show_output_atError(process_samtools_view, process_samtools_sort, process_stamtools_index, process_samtools_rmdup)
                dict_process[dir_name]['samtools_view'] = process_samtools_view
                dict_process[dir_name]['samtools_sort'] = process_samtools_sort
                dict_process[dir_name]['stamtools_index'] = process_stamtools_index
                dict_process[dir_name]['samtools_rmdup'] = process_samtools_rmdup

            if 'bcftools' not in skip_list:
                process_bcftools_mpileup, process_bcftools_call, process_view, process_tabix = run_bcftools(sub_dir, dir_name, pipeline_config)
                show_output_atError(process_bcftools_mpileup, process_bcftools_call,process_view, process_tabix)
                dict_process[dir_name]['bcftools_mpileup'] = process_bcftools_mpileup
                dict_process[dir_name]['bcftools_call'] = process_bcftools_call
                dict_process[dir_name]['bcftools_view'] = process_view
                dict_process[dir_name]['bcftools_tabix'] = process_tabix
            
            if 'variant_call' not in skip_list:
                process_medaka_variant = medaka_variant_call(sub_dir, fastq_file, pipeline_config)
                show_output_atError(process_medaka_variant)
                dict_process[dir_name]['medaka_variant'] = process_medaka_variant

                process_sniffles_variant = sniffles_variant_call(sub_dir, dir_name, pipeline_config)
                show_output_atError(process_sniffles_variant)
                dict_process[dir_name]['sniffles_variant'] = process_sniffles_variant

                process_parse_vcf_csv = parse_vcf_csv(sub_dir, dir_name, pipeline_config)
                show_output_atError(process_parse_vcf_csv)
                dict_process[dir_name]['parse_vcf_csv'] = process_parse_vcf_csv

            if 'copy' not in skip_list:
                process_copy, files_not_found, sampleid = copy_results(base_dir, fastq_file, sub_dir, dir_name, pipeline_config)
                show_output_atError(process_copy)
                dict_process[dir_name]['copy'] = process_copy    
                dict_results_files[sampleid] = files_not_found # add files not found to the list                                     

            # time.sleep(2) # wait for 2 seconds before next barcode, noneed if run in parallel
            # else:
            #     print(f'{dir_name} is not selected in sample config table, skipped.', file=sys.stdout)
            #     continue

    
    # show dict_results_files per barcode 
    if dict_results_files:
        print('=============================================================\n', file=sys.stdout)
        num_inc = 0
        for key, value in dict_results_files.items():
            if value:
                print(f'{key} Results files not generated: {value}', file=sys.stdout)
                num_inc += 1
            else:
                print(f'{key} All results files generated.', file=sys.stdout)
        print('Number of barcodes with incomplete results files: ', num_inc, '\n', file=sys.stdout)
    
    # zip all folders in base_dir/Results/ according to sampleid and Order No. in sample_config.csv, named as [Order No.].zip. 
    # put samples in the same order as in sample_config.csv into the same zip file
    zip_results(base_dir, pipeline_config['sample_config'])

    # save all processes in a dictionary if the process is not skipped
    # dump the dictionary to a pickle file
    print('Dumping processes stdout/err/returncode to pickle file... at ', f'{base_dir}/Results/processes.pkl \n', file=sys.stdout)
    with open(f'{base_dir}/Results/processes.pkl', 'wb') as f:
        pickle.dump(dict_process, f)

    # print current time
    print('All processes finished. Good luck! \n', file=sys.stdout)
    print('Current time: ', datetime.now(), file=sys.stdout)
    print('=============================================================', file=sys.stdout)


if __name__ == "__main__":
    args = parse_args()
    base_dir = args.working_dir
    config_file = args.config
    sample_config_file = args.sample_config
    no_guppy = args.no_guppy
    skip_list = args.skip
    
    global quiet 
    quiet = args.quiet # global variable for quiet mode
    main(base_dir, config_file = config_file, sample_config_file = sample_config_file, no_guppy = no_guppy, skip_list = skip_list, quiet = quiet)
