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
    validStops = set()
    for row in csv.reader(open(scorerFile), delimiter='\t'):
        if row[2] == "stop_codon":
            proteinEnd = extractAttribute(row, "proteinEnd")
            if proteinEnd == "1":
                parent = extractAttribute(row, "Parent")
                validStops.add(parent)
    return validStops


def convert(scorerFile):
    validStops = loadStopCodons(scorerFile)
    for row in csv.reader(open(scorerFile), delimiter='\t'):
        if row[2] == "mRNA":
            ID = extractAttribute(row, "ID")
            row[2] = "transcript"
            row[8] = f'transcript_id "{ID}"; gene_id "{ID}";'

            if ID in validStops:
                if row[6] == '+':
                    row[4] = str(int(row[4]) + 3)
                else:
                    row[3] = str(int(row[3]) - 3)

            gene = row.copy()
            gene[2] = "gene"
            gene[8] = f'gene_id "{ID}";'
            print("\t".join(gene))
            print("\t".join(row))
        elif row[2] == "CDS":
            score = extractAttribute(row, "eScore")
            parent = extractAttribute(row, "Parent")
            row[5] = score
            row[8] = f'transcript_id "{parent}"; gene_id "{parent}";'
            exon = row.copy()
            exon[2] = "exon"
            exon[7] = "."
            print("\t".join(exon))
            print("\t".join(row))


def main():
    args = parseCmd()
    convert(args.scorerFile)


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

    return parser.parse_args()


if __name__ == '__main__':
    main()
