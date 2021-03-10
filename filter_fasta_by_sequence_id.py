import sys
from Bio import SeqIO
import argparse

################

#Extract or remove sequences from a fasta file by sequence IDs

#Example usage: python filter_fasta_by_sequence_id.py -h

################

arg_parse = argparse.ArgumentParser(description='Extract or remove sequences from a fasta file by sequence IDs.')

arg_parse.add_argument('Input_fasta', type=str, nargs=1, help='Input fasta file (required)')
arg_parse.add_argument('Sequence_ID_list', type=str, nargs=1, help='Input file containing sequence IDs (required)')
arg_parse.add_argument('Output_fasta', type=str, nargs=1, help='Name of output fasta file (required)')
arg_parse.add_argument('--remove', action='store_const', default=1, const=0, help='Remove listed sequences (Default: extract listed sequences)')

args = arg_parse.parse_args()

fasta_file = args.Input_fasta[0]
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

id_file = args.Sequence_ID_list[0]
ids = open(id_file).read().rstrip("\n").split("\n")

outfile_name = args.Output_fasta[0]


id_list = []
output = []

print(len(seq_dict), "sequences read from input fasta file")

for id in ids:
    if id not in id_list:
         id_list.append(id)

print(len(id_list), "unique sequence IDs read from input ID list")

if args.remove == 1:
    for record in id_list:
        id = record.split()[0]

        if id not in seq_dict.keys():
            sys.stderr.write("Warning: ")
            sys.stderr.write(record)
            sys.stderr.write(" not found in fasta file\n")

        else:
            output.append(seq_dict[id])

elif args.remove == 0:
    for record in id_list:
        id = record.split()[0]

        if id not in seq_dict.keys():
            sys.stderr.write("Warning: ")
            sys.stderr.write(id)
            sys.stderr.write(" not found in fasta file\n")

        else:
            del seq_dict[id]

    for seq in seq_dict.values():
        output.append(seq)

SeqIO.write(output, outfile_name, "fasta")

print(len(output), "sequences written in output fasta")
