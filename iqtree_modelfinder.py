#!/usr/bin/env python3
'''
iqtree_modelfinder.py -- find the best scheme of evolution models for given alignments using modelfinder implemented in iqtree.

Author: yanpengch@qq.com
Date  : 2023-06-27
'''

import os
import sys
import argparse

parser = argparse.ArgumentParser(
                    prog='iqtree_modelfinder.py',
                    description=__doc__)

parser.add_argument('-i', '--input',
                    metavar='alignment.fna',
                    type=str,
                    nargs='+',
                    required=True,
                    help='multiple sequence alignment files that must be in FASTA format')
parser.add_argument('-o', '--outdir',
                    metavar='outdir',
                    type=str,
                    required=True,
                    help='specify output directory')
parser.add_argument('-m', '--model_restriction',
                    metavar='mrbayes|iqtree',
                    choices=['mrbayes', 'iqtree'],
                    default='mrbayes',
                    help='restrict search to models supported by other programs: mrbayes|iqtree')

args = parser.parse_args()

def get_partition(args_input_tuple):
    '''get the intervals of input alignments
    '''
    partition_len_dict = {}
    for alignment in args_input_tuple:
        partition_label = os.path.basename(alignment).split('.')[0]
        num_seq = 0
        len_alignment = 0
        with open(alignment, 'rt') as infh:
            for line in infh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    num_seq += 1
                    continue
                else:
                    if num_seq == 1:
                        len_alignment += len(line)
                    else:
                        break
        partition_len_dict[partition_label] = len_alignment
    return partition_len_dict

def run_moelfinder(args_input_tuple, args_model_restriction, args_outdir):
    '''to find the best evolution model for each input
    '''
    partition_model_dict = {}
    for alignment in args_input_tuple:
        partition_label = os.path.basename(alignment).split('.')[0]
        if args_model_restriction == 'mrbayes':
            iqtree_command = f'iqtree2 -s {alignment} -T 4 --prefix {args_outdir}/{partition_label}_modelfinder -m TESTONLY --mset mrbayes --msub nuclear'
        else:
            iqtree_command = f'iqtree2 -s {alignment} -T 4 --prefix {args_outdir}/{partition_label}_modelfinder -m TESTONLY --msub nuclear --redo'
        try:
            os.system(iqtree_command)
        except:
            print(f'failed run iqtree', file=sys.stderri, flush=True)
        with open(f'{args_outdir}/{partition_label}_modelfinder.log', 'rt') as infh:
            for line in infh:
                if line.startswith('Bayesian Information Criterion:'):
                    model = line.rstrip('\n').split()[-1]
                    break
        partition_model_dict[partition_label] = model
    return partition_model_dict

def get_best_scheme(partition_len_dict, partition_model_dict, args_outdir):
    '''
    output best partition scheme
    '''
    outfile_name = f'{args_outdir}/best_scheme.txt'
    outfh = open(outfile_name, 'wt')
    outfh.write('#nexus\n')
    outfh.write('begin sets;\n')

    end = 0
    for partition_label in partition_len_dict:
        start = end + 1
        end = end + partition_len_dict[partition_label]
        outfh.write(f'charset {partition_label} = {start}-{end};\n')

    charpartition = ''
    for partition_label in partition_model_dict:
        model = partition_model_dict[partition_label]
        charpartition += f'{model}:{partition_label}, '
    charpartition = 'charpartition ModelFinder = ' + charpartition.rstrip(', ') + ';'
    outfh.write(f'{charpartition}\n')
    outfh.write('end;\n')
    outfh.close()
    return print(f'Best scheme:{outfile_name}', file=sys.stdout, flush=True)

if __name__ == '__main__':
    partition_len_dict = get_partition(args.input)
    partition_model_dict= run_moelfinder(args.input, args.model_restriction, args.outdir)
    get_best_scheme(partition_len_dict, partition_model_dict, args.outdir)
    sys.exit(0)