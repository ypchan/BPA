#!/usr/bin/env python3
'''
iqtree_modelfinder.py -- find the best scheme of evolution models for given alignments using modelfinder implemented in iqtree.

Author: yanpengch@qq.com
Date  : 2023-06-27
'''

import os
import sys
import argparse

# convert iqtree models to mrbayes definitions
model_map = {'GTR':'nst=6',
             'GTR+I':'nst=6 rates=propinv',
             'GTR+G':'nst=6 rates=gamma',
             'GTR+I+G':'nst=6 rates=invgamma',
             'SYM':['nst=6', 'statefreqpr=fixed(equal)'],
             'SYM+I':['nst=6 rates=propinv', 'statefreqpr=fixed(equal)'],
             'SYM+G':['nst=6 rates=gamma','statefreqpr=fixed(equal)'],
             'SYM+I+G':['nst=6 rates=invgamma','statefreqpr=fixed(equal)'],
             'HKY':'nst=2',
             'HKY+I':'nst=2 rates=propinv',
             'HKY+G':'nst=2 rates=gamma',
             'HKY+I+G':'nst=2 rates=invgamma',
             'K2P':['nst=2', 'statefreqpr=fixed(equal)'],
             'K2P+I':['nst=2 rates=propinv', 'statefreqpr=fixed(equal)'],
             'K2P+G':['nst=2 rates=gamma','statefreqpr=fixed(equal)'],
             'K2P+I+G':['nst=2 rates=invgamma','statefreqpr=fixed(equal)'],
             'F81':'nst=1',
             'F81+I':'nst=1 rates=propinv',
             'F81+G':'nst=1 rates=gamma',
             'F81+I+G':'nst=1 rates=invgamma',
             'JC':['nst=1', 'statefreqpr=fixed(equal)'],
             'JC+I':['nst=1 rates=propinv', 'statefreqpr=fixed(equal)'],
             'JC+G':['nst=1 rates=gamma','statefreqpr=fixed(equal)'],
             'JC+I+G':['nst=1 rates=invgamma','statefreqpr=fixed(equal)']}

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

parser.add_argument('--mrbayes_nexus',
                    action='store_true',
                    help='generate a nexus file for running mrbayes')

parser.add_argument('--outgroup',
                    metavar='Outgroup',
                    help='specify the outgroup')

parser.add_argument('-m', '--model_restriction',
                    metavar='mrbayes|iqtree',
                    choices=['mrbayes', 'iqtree'],
                    default='mrbayes',
                    help='restrict search to models supported by other programs: mrbayes|iqtree')

args = parser.parse_args()

def concatenate_msa(msafile_tuple, args_outdir):
    '''input multiple MSA files, return a concatenated a python3 dictionary
    '''
    #1 get taxa list, and store sequence into dictionary
    all_msa_dict = {}
    taxa_lst = []
    for msa in msafile_tuple:
        all_msa_dict[msa] = {}
        with open(msa, 'rt') as fafh:
            for line in fafh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    taxa_id = line.lstrip('>')
                    taxa_lst.append(taxa_id)
                    all_msa_dict[msa][taxa_id] = []
                else:
                    all_msa_dict[msa][taxa_id].append(line)
    taxa_lst = list(set(taxa_lst))

    # multiple line fasta to single-line fasta
    for msa, fa_dict in all_msa_dict.items():
        all_msa_dict[msa] = {k: ''.join(v) for k, v in fa_dict.items()}

    #2 get the length of each msa files
    msa_length_dict = {}
    for msa, fa_dict in all_msa_dict.items():
        for taxa_id, seq in fa_dict.items():
            fa_length = len(seq)
            msa_length_dict[msa] = fa_length
            break

    #3 add missing data into msa
    for msa, fadict in all_msa_dict.items():
        msa_length = msa_length_dict[msa]
        for taxa_id in taxa_lst:
            if taxa_id not in fadict:
                all_msa_dict[msa][taxa_id] = '?' * msa_length

    #4 concatenate all MSA files
    concatenated_dict = {}
    for taxa_id in taxa_lst:
        concatenated_dict[taxa_id] = []
        for msa, fa_dict in all_msa_dict.items():
            seq = fa_dict[taxa_id]
            concatenated_dict[taxa_id].append(seq)

    #5 output concatenated FASTA file
    with open(args_outdir + '/concatenated.fna', 'wt') as outfh:
        for taxa_id, seq in concatenated_dict.items():
            outfh.write(f'>{taxa_id}\n{"".join(seq)}\n')
    return concatenated_dict

def mrbayes_template(partition_model_dict, partition_len_dict, concatenated_dict, args_outdir, args_outgroup='outgroup_label'):
    '''convert modelfinder models to mrbayes style
    '''
    mrbayes_file = f'{args_outdir}/run_mrbayes.nexus'
    ntaxa = len(concatenated_dict)
    nchar = sum(list(partition_len_dict.values()))
    concatenated_dict = {k: ''.join(seq_lst) for k, seq_lst in concatenated_dict.items()}

    with open(mrbayes_file, 'wt') as ofh:
        ofh.write('#NEXUS\n')
        ofh.write(f'  DIMENSIONS NTAX={ntaxa} NCHAR={nchar};\n')
        ofh.write(f'  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;\n')

        ofh.write(f'  MATRIX\n')
        longest_taxa_label_length = max([len(taxa_label) for taxa_label in concatenated_dict.keys()])

        start = 0
        end = 80
        while end < nchar:
            for taxa_label, sequence in concatenated_dict.items():
                seq_part = sequence[start:end]
                ofh.write(f'  {taxa_label:{longest_taxa_label_length}} {seq_part}\n')
            start = end
            end += 80
            ofh.write('\n')

        for taxa_label, sequence in concatenated_dict.items():
            seq_part = sequence[start:nchar]
            ofh.write(f'  {taxa_label:<{longest_taxa_label_length}} {seq_part}\n')

        ofh.write(';\n')
        ofh.write('END;\n')
        ofh.write('\n')
        ofh.write('BEGIN MRBAYES;\n')
        ofh.write('  log start filename=run_mrbayes.log;\n')
        ofh.write('  [! non-interactive,no prompts, no warnings]\n')
        ofh.write('  set autoclose=yes nowarnings=yes;\n')
        ofh.write('  [! add outgroup]\n')
        ofh.write(f'  outgroup {args_outgroup};\n')
        ofh.write('\n')

        start = 1
        end = 0
        num_subset = 0
        for _, partition_length in partition_len_dict.items():
            end += partition_length
            num_subset += 1
            ofh.write(f'  charset Subset{num_subset} = {start}-{end};\n')
            start += partition_length

        partition_string = ', '.join(['Subset' + str(i) for i in range(1, num_subset+1)])
        ofh.write(f'  partition ModelFinder = {num_subset}:{partition_string};\n')
        ofh.write('  set partition=ModelFinder;\n')
        ofh.write('\n')

        num_subset = 0
        for model in partition_model_dict.values():
            num_subset += 1
            model = model.rstrip('4').replace('+F', '')
            try:
                if len(model_map[model]) == 2:
                    ofh.write(f'  lset applyto=({num_subset}) {model_map[model][0]};\n')
                    ofh.write(f'  prset applyto=({num_subset}) {model_map[model][1]};\n')
                else:
                    ofh.write(f'  lset applyto=({num_subset}) {model_map[model]};\n')
            except:
                sys.exit(f'Error: unknown model {model}')

        ofh.write('\n')
        ofh.write('  prset applyto=(all) ratepr=variable;\n')
        ofh.write('  unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);\n')
        ofh.write('\n')

        ofh.write('  mcmc ngen=20000000 Stoprule=yes Stopval=0.01;\n')
        ofh.write('  sump;\n')
        ofh.write('  sumt;\n')
        ofh.write('  log stop;\n')
        ofh.write('END;\n')

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
    concatenated_dict = concatenate_msa(args.input, args.outdir)
    if args.mrbayes_nexus:
        mrbayes_template(partition_model_dict, partition_len_dict, concatenated_dict, args.outdir, args.outgroup)
    sys.exit(0)
