#!/usr/bin/env python
# author: Changlin Ke, NGS BI, Macrogen Europe

import pandas as pd
import time, sys, os, re
import subprocess
from datetime import datetime, timedelta



# find corresponding barcode id for each input barcode name
def get_selected_barcode(sample_config)->pd.DataFrame:
    # take the last two digits at the end of Barcode assigned col in sample config table
    sample_config['barcode_guppy'] = 'barcode' + sample_config['Barcode assigned'].str.extract('(\d+)', expand=False).str.zfill(2)
    return sample_config


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


# running a step of pipeline, reusable 
class StepRunner:
    def __init__(self, command:str, shell_use=False):
        self.command = command
        self.shell_use = shell_use
        self.stdout = None
        self.stderr = None
        self.returncode = None

    # run a shell tool with subprocess.run and return stdout, stderr, returncode
    def run_command(self)->list:
        start_time = time.time()
        start_datetime = datetime.now()
        print(f"Start date and time: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stdout)

        # running the cmd
        if not quiet: # a global var defined in main.py
            if self.shell_use:
                print("Running command: ",self.command, file=sys.stdout)
            else:
                print("Running command: "," ".join(self.command), file=sys.stdout)
        process = subprocess.run(self.command, check=False, capture_output=True, shell=self.shell_use)

        end_time = time.time()
        end_datetime = datetime.now()
        print(f"End date and time: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stdout)

        elapsed_time = end_time - start_time
        elapsed_hms = str(timedelta(seconds=elapsed_time))
        print(f"Elapsed time: {elapsed_hms}\n", file=sys.stdout)
        
        self.stdout, self.stderr, self.returncode = process.stdout, process.stderr, process.returncode
        # Check for errors and display error output if any
        self.show_output_atError()

    # show output at error
    def show_output_atError(self):
        if self.returncode != 0:
            print(f"Error: Command failed with return code {self.returncode}", file=sys.stderr)
            print(f"Stderr: {self.stderr.decode('utf-8')}", file=sys.stderr)

    # providing a way for calling code to access the output of the command after it has been run.
    def get_output(self):
        return self.stdout, self.stderr, self.returncode






