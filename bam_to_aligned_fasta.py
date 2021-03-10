#!/usr/bin/env python

import os, sys, re, argparse, collections
from subprocess import Popen, PIPE, DEVNULL
from Bio import SeqIO


#################

#Generate multiple sequence alignments in fasta format from bam file

#Usage example: python bam_to_aligned_fasta.py sample.bam ref.fa

#################


parser = argparse.ArgumentParser(description='Generate multiple sequence alignments in fasta format from bam file.')

parser.add_argument('input_bam', type=str, nargs=1, help='bam file produced by pbalign (required)')
parser.add_argument('reference_seq', type=str, nargs=1, help='fasta file containing the reference sequence used for pbalign (required)')
parser.add_argument('--min_insert_freq', type=float, default=0.2, metavar='<float>', help='minimum frequency required to call real insertions (default: 0.2)')
parser.add_argument('--noInsertionFilter', action='store_const', default=1, const=0, help='keep all insertions (default is to remove putatively artefactual insertions based on read frequency)')

args = parser.parse_args()

bam_file = args.input_bam[0]
ref_seq_fasta_file = args.reference_seq[0]
min_insert_freq = args.min_insert_freq


try:
    Popen(['samtools'], stdout=DEVNULL, stderr=DEVNULL).communicate()
except OSError:
    print('Cannot find samtools. Please add samtools to $PATH', file=sys.stderr)


if not os.path.exists(bam_file):
    print('Error:', bam_file, 'does not exist', file=sys.stderr)
    sys.exit(1)

if not os.path.exists(ref_seq_fasta_file):
    print('Error:', ref_seq_fasta_file, 'does not exist', file=sys.stderr)
    sys.exit(1)


ref_seq = list(SeqIO.parse(ref_seq_fasta_file, "fasta"))

if len(ref_seq) > 1:
    print('Error: Multiple reference sequences provided; only 1 allowed', file=sys.stderr)
    sys.exit(1)

ref_seq_len = len(ref_seq[0].seq)
ref_sequence = ref_seq[0].seq


samtools_view = ['samtools', 'view', bam_file]

proc = Popen(samtools_view, stdout=PIPE, stderr=PIPE, text=True).communicate()
sam_data = proc[0].rstrip("\n").split("\n")


if sam_data == ['']:
    print('Error: Input bam file is empty', file=sys.stderr)
    sys.exit(1)


out_prefix = os.path.splitext(bam_file)[0]


def read_cigar(cg):
    cg = re.sub('=', 'M', cg)
    cg_op = re.sub('[0-9]+', ' ', cg).split()
    cg_num = re.sub('[A-Z]+', ' ', cg).split()
    return(cg_op,cg_num)


read_dict = {}
insertion_dict = {}

for line in sam_data:

    line_split = line.split("\t")

    name = line_split[0]
    seq_id = ">" + name

    align_start = line_split[3]
    seq = line_split[9]
    cigar = line_split[5]

    [cg_op,cg_num] = read_cigar(cigar)

    if cg_op[0] == "S":
        seq = seq[int(cg_num[0]):]
        del cg_op[0]
        del cg_num[0]

    if cg_op[-1] == "S":
        seq = seq[:-int(cg_num[-1])]
        del cg_op[-1]
        del cg_num[-1]

    if align_start != 1:
        seq = "N"*(int(align_start) - 1) + seq

    pos = int(align_start) - 1

    for i in range(0,len(cg_num)):

        cig_char = cg_op[i]      

        if cig_char == "M":
            pos = pos + int(cg_num[i])
        elif cig_char == "X":
            pos = pos + int(cg_num[i])
        elif cig_char == "D":
            gap = "-"*int(cg_num[i])
            seq = seq[:pos] + gap + seq[pos:]
            pos = pos + int(cg_num[i])
        elif cig_char == "I":
            insertion = "_".join([str(pos), seq[pos:(pos+int(cg_num[i]))]])
            seq = seq[:pos] + seq[(pos+int(cg_num[i])):]

            if insertion in insertion_dict:
                insertion_dict[insertion].append(seq_id)
            else:
                insertion_dict[insertion] = [seq_id]

        else:
            print('Error: Found unexpected CIGAR character {}'.format(cig_char), file=sys.stderr)
            sys.exit(1)

    if len(seq) < ref_seq_len:
        seq = seq + "N"*(ref_seq_len - len(seq))

    if seq_id not in read_dict:
        read_dict[seq_id] = seq   
    else:
        print('Error: Found non-unique read ID: {}'.format(seq_id), file=sys.stderr)
        sys.exit(1)

#Filter insertions

real_insertions = {}
ins_size = {}

for ins in insertion_dict:

    read_count = len(insertion_dict[ins])

    ins_pos = ins.split("_")[0]

    if args.noInsertionFilter == 0:
        real_insertions[ins] = insertion_dict[ins]
    elif read_count > len(read_dict) * min_insert_freq:
        real_insertions[ins] = insertion_dict[ins]

for ins in real_insertions:

    ins_pos = ins.split("_")[0]
    ins_len = len(ins.split("_")[1])

    if ins_pos not in ins_size:
        ins_size[ins_pos] = ins_len
    elif ins_size[ins_pos] < ins_len:
        ins_size[ins_pos] = ins_len

for pos in collections.OrderedDict(sorted(ins_size.items(), key=lambda i: int(i[0]), reverse = True)):

    size = ins_size[pos]
    gap = "-"*size

    ref_sequence = ref_sequence[:int(pos)] + gap + ref_sequence[int(pos):]

    no_ins_reads = list(read_dict.keys())

    for ins, reads in real_insertions.items():

        ins_pos = ins.split("_")[0]

        if ins_pos == pos:

            no_ins_reads = [ r for r in no_ins_reads if r not in reads ]

            ins_seq = ins.split("_")[1]
            ins_gap = "-"*(size - len(ins_seq))

            for r in reads:
                read_seq = read_dict[r]
                read_seq_n = read_seq[:int(ins_pos)] + ins_seq + ins_gap + read_seq[int(ins_pos):]
                read_dict[r] = read_seq_n

    for r in no_ins_reads:
        read_seq = read_dict[r]
        read_seq_n = read_seq[:int(pos)] + gap + read_seq[int(pos):]
        read_dict[r] = read_seq_n


out_fasta_file = out_prefix + ".aligned.fa"
out_fasta = open(out_fasta_file, "w")

for seq_id, seq in read_dict.items():
    out_fasta.write(seq_id)
    out_fasta.write("\n")
    out_fasta.write(seq)
    out_fasta.write("\n")

out_fasta.close()
