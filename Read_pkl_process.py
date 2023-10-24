"""
Read in pickle file and show log for troubleshooting.
Author: Changlin Ke, NGS BI, MGEU

1. input: pickle file, list of barcode, list of tool, if no barcode and tool, use all
2. Read in pickle file:
each pickle contains a dictionary with format: {barcodes: {tool: [stderr, stdout, returncode]}}
3. print out the dictionary by barcodes and by tools for examination and troubleshooting

"""

import pickle
import os
import argparse

# parse arguments
def arg_parse():
    parser = argparse.ArgumentParser(description='Read in pickle file and show log for troubleshooting.')
    parser.add_argument('-p', '--pickle', help='pickle file path', required=True)
    parser.add_argument('-b', '--barcode', help='list of barcode', nargs='+', required=False)
    parser.add_argument('-t', '--tool', help='list of tool', nargs='+', required=False)
    args = parser.parse_args()
    return args

def print_process(process):
    print('stderr: ', process[0])
    print('stdout: ', process[1])
    print('returncode: ', process[2])
    print('\n')

if __name__ == '__main__':
    args = arg_parse()
# read in pickle file
    with open(args.pickle, 'rb') as f:
        data = pickle.load(f)

# print out the dictionary by barcodes and by tools for examination and troubleshooting
    if args.barcode:
        barcode_selected = args.barcode
    else:
        barcode_selected = list(data.keys())
    
    if args.tool:
        tool_selected = args.tool
    else:
        # all tools must present and be the same for all barcodes!
        # TypeError: 'dict_keys' object is not subscriptable
        tool_selected = list(data[barcode_selected[0]].keys())

    for barcode in barcode_selected:
        print('barcode: ', barcode)
        for tool in tool_selected:
            print('tool: ', tool)
            if data[barcode][tool] == None:
                print('No process recorded for this barcode and tool.')
            else:
                print_process(data[barcode][tool])
            print('=============\n')





