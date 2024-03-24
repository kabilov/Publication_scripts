#! /usr/bin/python3
""" Script trims poly-T tail by reference which obtained from molecular barcode """

import sys
import argparse


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta", type=str, required=False, help="Fasta file")
    parser.add_argument("-d", "--dict",  type=str, required=True, help="Dictionary fasta file (barcode   seq ref_name)")
    parser.add_argument("-e", "--error", type=int, required=False, default=0, help="Number of errors allowed (mismatch, deletion, insertion")
    parser.add_argument("-m", "--max",   type=int, required=False, default=0, help="Maximal length for output")
    parser.add_argument("-b", "--barcode", action="store_true", required=False, help="Add molecular barcode to seq name")
    parser.add_argument("-i", "--ignore", action="store_true", required=False, help="Remove reads with mismatches relative reference")

    parser.add_argument("-n", "--consts", type=str, required=True, help="file with consts (ref_name    5'-const    3'-const")
    parser.add_argument("-c", "--const5", type=str, required=False, help="5'-const seq")
    parser.add_argument("-s", "--split",  type=str, required=True, nargs=True, help="Seqs for splitting")
    parser.add_argument("-o", "--const3", type=str, required=False, help="3'-const seq")

    return parser.parse_args()


def split_seq(const, seq, errors):

    import regex
    from difflib import SequenceMatcher as SeqM

    #seqs = regex.split("({}){{e<={}}}".format(const, errors), seq, flags=regex.BESTMATCH) #.fuzzy_counts(1, 1, 1)
    seqs = regex.split(f"({const}){{e<={errors}}}", seq, flags=regex.BESTMATCH) #.fuzzy_counts(1, 1, 1)

    if "" in seqs:
        seqs.remove("")
    if len(seqs) < 3:
        return False

    # if number of splitted seqs more than 3 (short constant region was found in random regions (barcode2) also)
    if len(seqs) > 3:
        for j in range(-1,-len(seqs),-1):
            dl = sum([len(x) for x in seqs[j:]])
            if dl > 20:
                homolog = SeqM(None, seqs[j], const).ratio()
                if homolog > 0.8:
                    break
                else:
                    return False
            elif 20 > dl > 10:
                # check next seq (-1) for homology with const
                homolog = SeqM(None, seqs[j-1], const).ratio()
                if homolog > 0.8 and abs(len(const) - len(seqs[j-1])) < 2:
                    j = j - 1
                    break

        # join fragmented region
        seqs2 = ["".join(seqs[:j]), seqs[j], "".join(seqs[j+1:])]
        seqs = seqs2

    return seqs

def find_homology(ref, seq):
    # function finds homology between reference and sequence from 3'-end to 5'-end

    # http://python-lab.blogspot.com/2012/05/difflib.html
    import difflib

    # limit length of ref from 3'-side for exception homology with barcode1
    ref_cut = ref[-len(seq):]
    # if "TTTTT" in seq:
    #     ref_cut = ref_cut[seq.find("TTTTT"):]
    diff = difflib.ndiff(seq, ref_cut)
    newseq = [x.strip() for x in reversed(list(diff))]

    def find_score(bukva):
        if len(bukva) == 1:
            return 2
        else:
            return -2

    scores = []; seq2 = ""; score = 0
    for bukva in newseq:
        score += find_score(bukva)
        if len(bukva) != 1:
            if bukva[0] == "-":
                bukva = bukva[2]
            else:
                bukva = ""
        seq2 = bukva + seq2

        scores.append([seq2, score])

    # find max from the end of scores
    homology, max_score = max(scores[::-1], key=lambda x:x[1])

    match = difflib.SequenceMatcher(None, ref, homology).get_matching_blocks()
    match_len = sum([x[2] for x in match])
    if len(homology) - match_len > 2 or len(homology) < 5:
        return False, homology

    # check percent T in trimmed end
    trimmed_end = seq[:seq.find(homology)]
    if len(trimmed_end) == 0:
        return False, homology
    else:
        # trimmed_end2 = trimmed_end
        # if len(trimmed_end) > 15:
        #     trimmed_end2 = trimmed_end[-15:]
        if len(ref_cut) > len(homology):
            trimmed_end2 = trimmed_end[-5:]
            T_perc = 100 * trimmed_end2.count("T") / len(trimmed_end2)
            # if T_perc < 70 and trimmed_end2[-5:] != "TTGGG":
            if T_perc < 70:
                return False, homology
        elif len(homology) > len(ref_cut):
            return False, homology

    return True, homology

def find_homology_simple(ref, seq):
    # check homology between two seqs continuously

    ref2 = reversed(ref)
    seq2 = reversed(seq)
    i = 0; j = 0; m = 0
    for k,n in zip(ref2, seq2):
        m += 1
        if k == n:
            i += 1; j += 1
        else:
            if j == 0:
                m = m - 2; break
            j = 0
    if i > 0:
        seq3 = seq[-m:]
    else:
        return False

    return seq3


def main():

    args = parse_args()

    if args.fasta:
        file = open(args.fasta); stroki = file.readlines(); file.close()
    else:
        stroki = sys.stdin.readlines()

    # Read dictionary of molecular barcodes
    file = open(args.dict); stroki2 = file.readlines(); file.close()
    stroki2 = [x.rstrip("\n").split("\t") for x in stroki2]
    dict_barcodes = {x[0]:x[1:] for x in stroki2}

    # Read separate consts or from dictionary
    if args.consts:
        if args.const5 or args.const3:
            print("Arguments --const5 and --const3 are not incompatible with --consts !", file=sys.stderr); exit(1)

        # Read file with consts
        # ref_name    5'-const    split   3'-const
        file = open(args.consts); stroki3 = file.readlines(); file.close()
        stroki3 = [x.rstrip("\n").split("\t") for x in stroki3]
        dict_refs = {x[0]: x[1:] for x in stroki3}
    else:
        if not args.const5:
            print("Const5 is absent!", file=sys.stderr); exit(1)
        elif not args.split:
            print("Split is absent!", file=sys.stderr); exit(1)
        elif not args.const3:
            print("Const3 is absent!", file=sys.stderr); exit(1)
        else:
            const5 = args.const5
            const3 = args.const3
    bool = False
    stat = [0, 0, 0, 0, 0]
    dict_ref_stat = {}
    stroki.append(">")
    name = stroki[0].rstrip("\n"); seq = ""
    for stroka in stroki[1:]:
        stroka = stroka.rstrip("\n")
        if stroka[0] == ">":
            stat[0] += 1    # stat1: number of reads in general
            for split in args.split:
                seqs = split_seq(split, seq, args.error)
                if seqs:
                    stat[1] += 1    # stat2: number of splitted reads
                    # check length of barcode2
                    barcode2 = seqs[2]
                    if not 20 >= len(barcode2) >= 9:
                        continue
                    stat[2] += 1    # stat3: number of reads with right length of barcode2
                    # check barcode2 in dictionary
                    if barcode2 in dict_barcodes:
                        stat[3] += 1    # stat4: number of reads with barcode from dictionary
                        if args.consts:
                            ref_name = dict_barcodes[barcode2][1]
                            const5 = dict_refs[ref_name][0]
                            const3 = dict_refs[ref_name][1]
                            # stats by refs
                            if ref_name in dict_ref_stat:
                                dict_ref_stat[ref_name][0] += 1
                            else:
                                dict_ref_stat[ref_name] = [1, 0, []]
                        if split not in const3:
                            continue    # Rare event when ref type is wrong
                        ref = const5 + dict_barcodes[barcode2][0] + const3
                        trimmed = find_homology(ref, seqs[0]+seqs[1])
                        if trimmed[0]:
                            stat[4] += 1    # stat5: number of trimmed reads
                            dict_ref_stat[ref_name][1] += 1
                            dict_ref_stat[ref_name][2].append(seqs[2])
                            if args.barcode:
                                name += "barcode=" + seqs[2]
                            if args.consts:
                                name += ";ref=" + ref_name
                            if bool:
                                print("the same seq is splitted different const twice!!!", file=sys.stderr)
                                exit(1)
                            else:
                                print(name); print(trimmed[1] + seqs[1][:-6])
                            bool = True
            bool = False
            name = stroka; seq = ""
        else:
            seq += stroka

    print(f"Total={stat[0]}\tSplitted={stat[1]}\tWith_norm_len_barcode2={stat[2]}\tSeq_with_barcodes_in_dict={stat[3]}\tTrimmed={stat[4]}\n", file=sys.stderr)
    if dict_ref_stat:
        print("Ref\tWith_right_barcode2\tSplitted\tUnique", file=sys.stderr)
        ref_stat = [[x] + dict_ref_stat[x][:2] + [len(set(dict_ref_stat[x][2]))] for x in dict_ref_stat]
        [print("\t".join([str(x) for x in y]), file=sys.stderr) for y in ref_stat]

    return 0


if __name__ == '__main__':
    sys.exit(main())
