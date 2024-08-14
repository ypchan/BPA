#!/usr/bin/env python3

# update: 2024-02-01
# bugs  : yanpengch@qq.com

import os
import sys
import argparse

# check the requirements first.
if not os.system('cd-hit -h &> /dev/null'):
    print('Error: cd-hit is required. Please install it.', file=sys.stderr, flush=True)
    sys.exit(1)
if os.system('blastn -h &> /dev/null'):
    print('Error: blastn is required. Please install it.', file=sys.stderr, flush=True)
    sys.exit(1)

# add command-line argument parser
parser = argparse.ArgumentParser(description="extract barcode sequences from genome assemblies.")
parser.add_argument('--sequence', '-s',
    required=True,
    action="append",
    help='barcode sequences in multi-fasta format.')

parser.add_argument('--barcode_label', '-b',
    action="append",
    required=True,
    help='barcode labels. for instances: ITS TUB TEF RPB2')

parser.add_argument('--flank_length', '-f',
    type=int,
    default=100,
    help='length of flanking sequences surrounded the blastn hit region.')

parser.add_argument('--genome', '-g',
    action="append",
    help='genume assembles that you want to extracted barcode seqeucnes from')

parser.add_argument('--genome', '-g',
    required=True,
    help='genume assembles that you want to extracted barcode seqeucnes from')


parser.add_argument('--prefix', '-p',
    type=str,
    help='prefix the output: prefix_ITS.fasta')

args = parser.parse_args()

# confirm that each barcode sequence have the corresponding barcode label
def check_input(args_sequence, args_barcode_label):
    if len(args_barcode_label) != len(args_sequence):
        sys.exit('Error: Each barcode sequence must be associated with a corresponding barcode label.' )

def run_cdhit(args_sequence, args_barcode_label):
    os.system('rm -rf temp_barcode_fisher')
    os.system('mkdir temp_barcode_fisher')
    so.system('mkdir temp_barcode_fisher/temp_cdhit')
    print(args_sequence, args_barcode_label)
    for fasta, label in zip(args_sequence, args_barcode_label):
        print(fasta, label)
        os.system(f"cd-hit -i {fasta} -aL 0.8 -o temp_barcode_fisher/temp_cdhit/{label}.cdhit.I0.9L0.8 &>/dev/null")

    for label in args_barcode_label:
        with open(f"temp_barcode_fisher/temp_cdhit/{label}.cdhit.I0.9L0.8", 'rt') as infh, open(f"temp_barcode_fisher/temp_cdhit/{label}.cdhit.I0.9L0.8.fasta", 'wt') as outfh:
            count = 0
            for line in infh:
                if line.startswith('>'):
                    count += 1
                    outfh.write(f">{label}_{count}\n")
                    seq_lst = []
                else:
                    outfh.write(f"{line}")

def run_blastn(args_barcode_label, args_genome):
    os.system('mkdir temp_barcode_fisher/temp_blastn')
    for label in args_barcode_label:
        for genome in args_genome:
            genome_abbre = os.path.basename(genome)
            os.system(f"blastn -query temp_barcode_fisher/temp_cdhit/{label}.cdhit.I0.9L0.8.fasta -subject {genome} -outfmt 6 -out temp_barcode_fisher/temp_blastn/{label}_{genome_abbre}.blastn")

def read_genome_2dict(genome):
    genome_dict = {}
    with open(genome, 'rt') as infh:
        for line in infh:
            if line.startswith('>'):
                contig_id = line.rstrip('\n').split()[0].lstrip('>')
                genome_dict[contig_id] = []
            else:
                genome_dict[contig_id].append(line.rstrip('\n'))
    genome_dict = {k:''.join(v) for k,v in genome_dict.items()}
    return genome_dict
    
def check_args_genome(args_genome):
    genome_lst = []
    is_genome_list = True
    with open(args_genome, 'rt') infh:
        for line in infh:
            if line.startswith('>'):
                is_genome_list = False
                break
    if is_genome_list:
        with open(args_genome, 'rt') infh:
            for line in infh:
                genome_lst.append(line.rstrip('\n'))
    else:
        genome_lst = [args_genome]
     return genome_lst

def extract_barcode_seq(args_barcode_label, genome_lst, args_flank_length):
    
    barcode_dict = {}
    for barcode in args_barcode_label:
        barcode_dict[barcode] = {}
        for genome in genome_lst:
            genome_abbre = os.path.basename(genome)
            if os.stat(f"temp_barcode_fisher/temp_blastn/{barcode}_{genome_abbre}.blastn").st_size == 0:
                barcode_dict[barcode][genome_abbre] = 0
            else:
                genome_dict = read_genome_2dict(genome)
                with open(f"temp_barcode_fisher/temp_blastn/{barcode}_{genome_abbre}.blastn", 'rt') as infh:
                    line_lst = infh.readlines()[0].split()

                    contig_id = line_lst[1]
                    sstart = int(line_lst[8])
                    send = int(line_lst[9])

                if sstart > send:
                    sstart, send = send, sstart
                if (sstart - args_flank_length) <0:
                    sstart = 0
                else:
                    sstart -= args_flank_length
                if (send - args_flank_length) > len(genome_dict[contig_id]):
                    send = len(genome_dict[contig_id])
                else:
                    send += args_flank_length

                barcode_dict[barcode][genome_abbre] = genome_dict[contig_id][sstart-1:send]
    return barcode_dict

def output(barcode_dict, args_prefix):
    for barcode in barcode_dict:
        if args.prefix:
            outfile = f"{args_prefix}_{barcode}.fasta"
        else:
            outfile = f"{barcode}.fasta"
        with open(outfile, 'wt') as ofh:
            for genome_abbrev, seq in barcode_dict[barcode].items():
                genome_abbrev = genome_abbrev.strip('.fasta').strip('.fna')
                if seq != 0:
                    ofh.write(f">{genome_abbrev}\n{seq.upper()}\n")

if __name__ == '__main__':
    check_input(args.sequence, args.barcode_label)
    run_cdhit(args.sequence, args.barcode_label)
    genome_lst = check_args_genome(args_genome)
    run_blastn(args.barcode_label, genome_lst)
    barcode_dict = extract_barcode_seq(args.barcode_label, args.genome, args.flank_length)
    output(barcode_dict, args.prefix)
    sys.exit(0)
