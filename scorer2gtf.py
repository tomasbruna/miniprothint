#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Convert the output of the miniprot boundary scorer to a miniprot-
# like gtf. Warning: the CDS scores will not perfectly match the
# native miniprot.gtf as the miniprot boundary scorer has a
# different way of penalizing gaps and frameshifts. However, the
# differences are usually small. The gene scores will match
# perfectly -- they are parsed directly form the miniprot output.
# ==============================================================


import argparse
import csv
import re


def extractAttribute(row, feature):
    regex = feature + '=([^;]+)'
    return re.search(regex, row[8]).groups()[0]


def loadStopCodons(scorerFile):
    allStops = {}
    validStops = set()
    for row in csv.reader(open(scorerFile), delimiter='\t'):
        if row[2] == "stop_codon":
            parent = extractAttribute(row, "Parent")
            prot = extractAttribute(row, "prot")
            allStops[f'{parent}_{prot}'] = (int(row[3]), int(row[4]))

            proteinEnd = extractAttribute(row, "proteinEnd")
            if proteinEnd == "1":
                validStops.add(f'{parent}_{prot}')

    return allStops, validStops


def convert(scorerFile, stopsInCDS):
    allStops, validStops = loadStopCodons(scorerFile)
    for row in csv.reader(open(scorerFile), delimiter='\t'):
        if row[2] == "mRNA":
            ID = extractAttribute(row, "ID")
            prot = extractAttribute(row, "prot")
            row[2] = "transcript"
            row[8] = f'transcript_id "{ID}_{prot}"; gene_id "{ID}_{prot}";'

            if f'{ID}_{prot}' in validStops:
                if row[6] == '+':
                    row[4] = str(int(row[4]) + 3)
                else:
                    row[3] = str(int(row[3]) - 3)

            gene = row.copy()
            gene[2] = "gene"
            gene[8] = f'gene_id "{ID}_{prot}";'
            print("\t".join(gene))
            print("\t".join(row))
        elif row[2] == "CDS":
            score = extractAttribute(row, "eScore")
            parent = extractAttribute(row, "Parent")
            prot = extractAttribute(row, "prot")

            if stopsInCDS and f'{parent}_{prot}' in allStops:
                if row[6] == '+':
                    if int(row[4]) + 1 == allStops[f'{parent}_{prot}'][0]:
                        row[4] = str(int(row[4]) + 3)
                else:
                    if int(row[3]) - 1 == allStops[f'{parent}_{prot}'][1]:
                        row[3] = str(int(row[3]) - 3)

            row[5] = score
            row[8] = f'transcript_id "{parent}_{prot}"; gene_id "{parent}_{prot}";'
            exon = row.copy()
            exon[2] = "exon"
            exon[7] = "."
            print("\t".join(exon))
            print("\t".join(row))


def main():
    args = parseCmd()
    convert(args.scorerFile, args.stopsInCDS)


def parseCmd():

    parser = argparse.ArgumentParser(description='Convert the output of the \
        miniprot boundary scorer to a miniprot-like gtf. Warning: the CDS\
        scores will not perfectly match the native miniprot.gtf as the\
        miniprot boundary scorer has a different way of penalizing gaps and \
        frameshifts. However, the differences are usually small. The gene\
        scores will match perfectly -- they are parsed directly form the \
        miniprot output.')

    parser.add_argument('scorerFile', metavar='miniprot_parsed.gff', type=str,
                        help='The input gff to collapse.')

    parser.add_argument('--stopsInCDS',  default=False, action='store_true',
                        help='Extend terminal CDS to include stop codons.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
