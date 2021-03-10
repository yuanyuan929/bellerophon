#!/usr/bin/env python

import os, sys, re, operator, argparse, collections
from subprocess import Popen, PIPE, DEVNULL
import pandas as pd
from Bio import SeqIO


#################

#Call alleles for individual sample

#Usage example: python bellerophon.py sample.bam ref.fasta --db db.fa --blacklist blacklist.txt

#################


parser = argparse.ArgumentParser(description='Call alleles from gene amplicon data generated on PacBio system.')

parser.add_argument('input_bam', type=str, nargs=1, help='bam file produced by pbalign (required)')
parser.add_argument('reference_seq', type=str, nargs=1, help='fasta file containing the reference sequence used for pbalign (required)')
parser.add_argument('--db', type=str, nargs=1, default=1, metavar='known_allele_db', help='fasta file containing known alleles (optional)')
parser.add_argument('--blacklist', type=str, nargs=1, default=1, metavar='blacklist', help='region(s) in the sequences to exclude from variable site calling (optional)')
parser.add_argument('--min_var_freq', type=float, default=0.05, metavar='<float>', help='minimum frequency required to call nucleotide variant (default: 0.05)')
parser.add_argument('--min_read_perc', type=float, default=0.01, metavar='<float>', help='minimum percentage of reads required to consider a sequence variant as a candidate (default: 0.01)')
parser.add_argument('--evidence', action='store_const', default=1, const=0, help='generates an additional output file listing read IDs of each allele')


args = parser.parse_args()

bam_file = args.input_bam[0]
ref_seq_fasta_file = args.reference_seq[0]
min_var_freq = args.min_var_freq
min_read_perc = args.min_read_perc



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

if args.db != 1:
    if not os.path.exists(args.db[0]):
        print('Error:', args.db[0], 'does not exist', file=sys.stderr)
        sys.exit(1)

if args.blacklist != 1:
    if not os.path.exists(args.blacklist[0]):
        print('Error:', args.blacklist[0], 'does not exist', file=sys.stderr)
        sys.exit(1)


ref_seq = list(SeqIO.parse(ref_seq_fasta_file, "fasta"))

if len(ref_seq) == 0:
    print('Error: No reference sequence found', file=sys.stderr)
    sys.exit(1)

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

if args.db != 1:
    allele_db_file = args.db[0]
    allele_db = list(SeqIO.parse(allele_db_file, "fasta"))

if args.blacklist != 1:
    blacklist = open(args.blacklist[0]).read().rstrip("\n").split("\n")


out_prefix = os.path.splitext(bam_file)[0]


def read_cigar(cg):
    cg = re.sub('=', 'M', cg)

    cg_op = re.sub('[0-9]+', ' ', cg).split()

    cg_num = re.sub('[A-Z]+', ' ', cg).split()

    return(cg_op,cg_num)


def trim_seq(seq):
    trimmed_seq = []

    for pos in var_list:
        trimmed_seq.append(seq[pos])

    return(trimmed_seq)


def restore_seq(variable_sites):
    seq = rep_nt_list

    for i in range(0,len(var_list)):
        seq[var_list[i]] = variable_sites[i]

    seq = "".join(seq)

    return(seq)


def check_allele(seq):
    allele_name = []

    for record in allele_db:
        if seq in record.seq:
            allele_name.append(record.name)

    allele = ";".join(allele_name)

    return(allele)


def degap(seq):
     seq = seq.replace("-", "")
     seq = seq.replace("N", "")
     return(seq)


print("Generating multiple sequence alignments")

#seq_list = []
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

    #Insert gaps at the end of sequence for soft-clipped reads
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

    #########
    #Regions containing long streches of a single type of nucleotides may need to be blocked due to high rates of indel errors caused by PCR & sequencing which can interfere with allele calling
    #
    #if int(ins_pos) < 15 or 375 < int(ins_pos) < 384:
    #    continue
    #
    if args.blacklist != 1:

        bl_flag = 0

        for n in range(0,len(blacklist)):
            bl_start = int(blacklist[n].split("\t")[0]) - 1
            bl_end = int(blacklist[n].split("\t")[1])

            if bl_start <= int(ins_pos) < bl_end:
                bl_flag += 1

        if bl_flag > 0:
            continue
        
    #########

    if read_count > len(read_dict) * min_var_freq * 4:
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


print("Generating list of variable sites")

seq_list = list(read_dict.values())

variable_sites = []
rep_nt_list = []

for i in range(0,len(seq_list)):
    seq_list[i] = " ".join(seq_list[i]).split(" ")

seq_list_t = pd.DataFrame(seq_list).transpose().values.tolist()

for i in range(0,len(seq_list_t)):
    A=seq_list_t[i].count("A")
    C=seq_list_t[i].count("C")
    G=seq_list_t[i].count("G")
    T=seq_list_t[i].count("T")
    gap=seq_list_t[i].count("-")

    if A+C+G+T+gap != 0:

        fA = float(A)/float(A+C+G+T+gap)
        fC = float(C)/float(A+C+G+T+gap)
        fG = float(G)/float(A+C+G+T+gap)
        fT = float(T)/float(A+C+G+T+gap)
        fgap = float(gap)/float(A+C+G+T+gap)

        freq = [("A",fA), ("C",fC), ("G",fG), ("T",fT), ("-",fgap)]
        freq.sort(key=operator.itemgetter(1), reverse = True)

        if freq[1][1] >= min_var_freq:
            variable_site = "\t".join([str(i+1), "%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (fA, fC, fG, fT, fgap)])
            variable_sites.append(variable_site)

        rep_nt = freq[0][0]

    else:
        rep_nt = "N"

    rep_nt_list.append(rep_nt)


out_alleles_file = out_prefix + ".alleles.fa"
out_alleles = open(out_alleles_file, "w")

if not variable_sites:

    print("Warning: no variable sites detected")
    print("Writing allele fasta file")

    #seq_degap = [ n for n in rep_nt_list if n != "-" and n != "N" ]

    seq = degap("".join(rep_nt_list))

    allele_name = []

    if args.db != 1:
        allele_name = check_allele(seq)

    if not allele_name:
        seq_title = "_".join([">allele_1", str(len(seq_list)), "reads"])
    else:
        seq_title = "_".join([">allele_1", str(len(seq_list)), "reads", allele_name])

    out_alleles.write(seq_title)
    out_alleles.write("\n")
    out_alleles.write(seq)
    out_alleles.write("\n")
    out_alleles.close()

    print("Run completed")

    sys.exit(0)


out_variableSites_file = out_prefix + ".variableSites.txt"
out_variableSites = open(out_variableSites_file, "w")

out_variableSites.write("\t".join(["site", "freq_A", "freq_C", "freq_G", "freq_T", "freq_gap"]))
out_variableSites.write("\n")

for line in variable_sites:
    out_variableSites.write(line)
    out_variableSites.write("\n")

out_variableSites.close()


#use variable sites only for read analysis

var_list = []

for line in variable_sites:
    site = line.split("\t")[0]
    var_list.append(int(site)-1)

trimmed_seq_list = []

for i in range(0,len(seq_list)):
    seq_list[i] = "".join(seq_list[i])

for seq in seq_list:

    trimmed_seq = trim_seq(seq)

#    trimmed_seq = []

#    for pos in var_list:
#        trimmed_seq.append(seq[pos])

    if 'N' not in trimmed_seq:
        trimmed_seq = "".join(trimmed_seq)
        trimmed_seq_list.append(trimmed_seq)


trimmed_seq_counts = collections.Counter(trimmed_seq_list)
sorted_trimmed_seq_counts = sorted(trimmed_seq_counts.items(), key=operator.itemgetter(1), reverse = True)

candidate_seq_list = []


#filter by read support

for line in sorted_trimmed_seq_counts:

    total_reads = len(trimmed_seq_list)

    if line[1] > int(total_reads)*min_read_perc and line[1] >= 3:
        candidate_seq_list.append(line)


#remove chimeric sequences

real_allele_list = candidate_seq_list[0:2]
chimeric_seq_list = []
read_pile = candidate_seq_list[0:2]

for candidate in candidate_seq_list[2:]:

    candidate_seq = candidate[0]

    a_flag_list = []
    b_flag_list = []

    for read in read_pile:

        read_seq = read[0]

        a_flag = 0
        b_flag = 0

        for i in range(1,len(candidate_seq)):

            a = candidate_seq[:i]
            b = candidate_seq[i:]

            if read_seq.startswith(a):
                a_flag += 1

            if read_seq.endswith(b):
                b_flag += 1

        a_flag_list.append(a_flag)
        b_flag_list.append(b_flag)

    head = sorted(a_flag_list, reverse = True)[0]
    tail = sorted(b_flag_list, reverse = True)[0]

#    if sorted(a_flag_list, reverse = True)[0] + sorted(b_flag_list, reverse = True)[0] >= len(candidate_seq) and \
#       (var_list[sorted(a_flag_list, reverse = True)[0]] - var_list[len(candidate_seq)-sorted(b_flag_list, reverse = True)[0]-1]) > 50:
    if head + tail >= len(candidate_seq) and \
       (var_list[head] - var_list[len(candidate_seq)-tail-1]) > 20:

        chimeric_seq_list.append(candidate)

    else:
        real_allele_list.append(candidate)

    read_pile.append(candidate)


print("Writing allele fasta file")

n = 1
for i in range(0,len(real_allele_list)):

    real_allele_list[i] = list(real_allele_list[i])

    allele = real_allele_list[i]

    read_count = allele[1]
    variable_sites = " ".join(allele[0]).split(" ")

    seq = restore_seq(variable_sites)
    seq = degap(seq)

    allele_name = []

    if args.db != 1:
        allele_name = check_allele(seq)

    if not allele_name:
        allele_id = "_".join([">allele", str(n), str(read_count), "reads"])
    else:
        allele_id = "_".join([">allele", str(n), str(read_count), "reads", allele_name])

    real_allele_list[i].insert(1, allele_id)

    n += 1

    out_alleles.write(allele_id)
    out_alleles.write("\n")
    out_alleles.write(seq)
    out_alleles.write("\n")

out_alleles.close()


out_chimera_file = out_prefix + ".chimera.fa"
out_chimera = open(out_chimera_file, "w")

n = 1
for chimera in chimeric_seq_list:

    read_count = chimera[1]
    variable_sites = " ".join(chimera[0]).split(" ")

    seq_title = "_".join([">chimera", str(n), str(read_count), "reads"])
    n += 1

    seq = restore_seq(variable_sites)
    seq = degap(seq)

    out_chimera.write(seq_title)
    out_chimera.write("\n")
    out_chimera.write(seq)
    out_chimera.write("\n")

out_chimera.close()


if args.evidence == 0:

    print("Writing allele evidence file")

    for seq_id, seq in read_dict.items():

        trimmed_seq = "".join(trim_seq(seq))

        seq_id = seq_id.replace(">", "")

        for i in range(0,len(real_allele_list)):
            if trimmed_seq == real_allele_list[i][0]:
                real_allele_list[i].append(seq_id)

    out_allele_evidence_file = out_prefix + ".allele_evidence.txt"
    out_allele_evidence = open(out_allele_evidence_file, "w")

    out_allele_evidence.write("\t".join(["allele", "read_IDs"]))
    out_allele_evidence.write("\n")

    for allele in real_allele_list:

        allele_id = allele[1].replace(">", "")

        #read_count = allele[2]

        read_ids = ",".join(allele[3:])

        evidence_line = "\t".join([allele_id, read_ids])

        out_allele_evidence.write(evidence_line)
        out_allele_evidence.write("\n")

    out_allele_evidence.close()


print("Run completed")
