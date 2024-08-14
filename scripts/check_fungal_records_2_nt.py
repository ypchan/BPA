#!/usr/bin/env python3
import sys
import argparse
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Check fungal records in GenBank nt database.')
parser.add_argument('--species',
                    nargs='+',
                    type=str,
                    required=True,
                    help='Species name(s) to query. For example: Endocalyx_cinctus. No space allowed.')
parser.add_argument('--country',
                    type=str,
                    help='if give the country name: only the record have this inofrmation can be present.')
parser.add_argument('--host',
                    type=str,
                    help='if give the host name: only the record from this host can be present.')
args = parser.parse_args()


def query_genbank(args_species):
    # query GenBank using esearch and efetch
    # define searching expressions
    species_name = args_species.replace("_", " ")
    search_expression = f"\"{species_name}[Organism]) AND (fungi[filter] AND biomol_genomic[PROP] AND is_nuccore[filter] AND (100[SLEN] : 2000[SLEN])\""
    query_command = f"esearch -db nucleotide -query {search_expression} | efetch -format gb > tmp.gbk"
    process = subprocess.Popen(query_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()

    if process.returncode != 0:
        print(f"Error querying GenBank: {err.decode()}", file=sys.stderr)
        return None

    return

def parse_genbank_data(species_name):
    records = []
    for record in SeqIO.parse('tmp.gbk', "genbank"):
        # pase every record one by one
        submit_date = record.annotations.get('date', '')
        organism = record.annotations.get('organism', '')
        country = ';'.join(record.features[0].qualifiers.get('country', ''))
        culture_collection = ';'.join(record.features[0].qualifiers.get('culture_collection', ''))
        note = ';'.join(record.features[0].qualifiers.get('note', ''))
        specimen_voucher = ';'.join(record.features[0].qualifiers.get('specimen_voucher', ''))
        host = ';'.join(record.features[0].qualifiers.get('host', ''))

        records.append([species_name, submit_date, organism, country, host, note, culture_collection, specimen_voucher])
    return records

if __name__ == "__main__":
    out_header = 1
    for species in args.species:
        query_genbank(species)
        try:
            records = parse_genbank_data(species)
        except:
            print('No records.', file=sys.stdout, flush=True)
            sys.exit(0)
        head_lst = ['#query_name', 'submit_date', 'organism', 'country','host', 'note', 'culture_collection','specimen_voucher']
        if out_header == 1:
            print('\t'.join(head_lst), file=sys.stdout, flush=True)
            out_header = 0
        for rec in records:
            print('\t'.join(rec), file=sys.stdout, flush=True)
