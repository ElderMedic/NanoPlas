#!/usr/bin/env python
# author: Changlin Ke, NGS BI, Macrogen Europe
# date: 2023-10-31, v1.4
# history: 2023-10-30, v1.3 preview

import os, re, sys, time
import subprocess
import pandas as pd

'''
main function of pipeline

work to do:
1. add order number to sample config. copy and zip files based on order number
2. basecalling qc notebook imported
3. improve readability of std.out/error

'''

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
## ver 1.3 alpha preview.                                                         
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
## Usage: at base environment (to be updated!), run with python3.10
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

