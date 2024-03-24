#! /usr/bin/python3
""" Find and group the same seq under ribosomes """

import sys
import argparse


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta",  type=str, required=False, help="Fasta file with molecular barcodes in names")
    parser.add_argument("-d", "--db",     type=str, required=True,  help="File with molecular barcodes database")
    parser.add_argument("-c", "--consts", type=str, required=True,  help="file with consts (ref_name    5'-const    split   3'-const")
    parser.add_argument("-M", "--max",    type=int, required=True,  help="Maximal length for output")
    #parser.add_argument("-o", "--output", required=False, default="file", help="Output name")

    return parser.parse_args()


def translate(seq, table_num):

    from Bio.Seq import Seq
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[table_num]
    seq = seq[:len(seq)//3*3]
    coding_dna = Seq(seq)
    return coding_dna.translate(table)


def coords(numbers, dlina):

    j = []
    # numbers = [63-x for x in numbers]
    for i in range(-1, dlina):
        count = numbers.count(i)
        j.append(count)

    #z = "\t".join([str(y) for y in j])
    z = [y for y in j]
    return z


def main():

    args = parse_args()

    with open(args.db) as file1:
        stroki1 = file1.readlines()
    if args.fasta:
        with open(args.fasta) as file:
            stroki2 = file.readlines()
    else:
        if not sys.stdin.isatty():
            print("Stdin is empty!", file=sys.stderr)
            return 1
        stroki2 = sys.stdin.readlines()
    with open(args.consts) as file3:
        stroki3 = file3.readlines()

    # create dictionary with molecular barcodes
    stroki1 = [x.rstrip("\n").split("\t") for x in stroki1]
    dict_mb = {x[0]: x[1] for x in stroki1}
    # create dictionary of consts
    stroki3 = [x.rstrip("\n").split("\t") for x in stroki3]
    dict_consts = {x[0]: x[1:] for x in stroki3}

    # extract seq by molecular barcode in names
    dict_ref = {}
    for stroka in stroki2:
        y = stroka.rstrip("\n")
        if stroka[0] == ">":
            # read info from read name
            coord = y.rfind("barcode=")
            mb = y[coord+8:y.rfind(";")]
            sample = y[y.find("sample=")+7:coord-1]
            ref_name = y[y.rfind(";")+5:]
        else:
            seq = y
            dlina = len(seq)
            const = dict_consts[ref_name][1]
            ref = dict_mb[mb] + const
            if ref not in dict_ref:
                prot = str(translate(ref, 11))
                dict_ref[ref] = [ref_name, ref, prot, mb, sample, []]
            new_coord = len(ref) - dlina + 1
            if new_coord < 0:
                dict_ref[ref][-1].append(-1)        # count before AUG codon
            else:
                dict_ref[ref][-1].append(new_coord)
            mb = ""

    ref_names = sorted([x for x in dict_ref])

    for ref in ref_names:
        dict_ref[ref][-1] = coords(sorted(dict_ref[ref][-1]), args.max)
    # summa = sum([sum(dict_ref[x][-1][1:]) for x in ref_names])  # [1:] = exclude count -1 (before ATG)
    for ref in ref_names:
        z = dict_ref[ref]
        # z[-1] = "\t".join([str(x*100/summa) for x in z[-1]])
        z[-1] = "\t".join([str(x) for x in z[-1]])
        print("\t".join(z))

    return 0


if __name__ == '__main__':
    sys.exit(main())
