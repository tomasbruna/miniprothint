#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Collapse introns, starts, and stops and compute their coverage/depth.
# ==============================================================


import argparse
import csv
import re


def signature(row):
    return f'{row[0]}_{row[2]}_{row[3]}_{row[4]}_{row[6]}'


def extractAttribute(row, feature):
    regex = feature + '=([^;]+)'
    return re.search(regex, row[8]).groups()[0]


class Feature:
    def __init__(self, row):
        self.row = row
        self.count = 1
        self.prots = [extractAttribute(row, "prot")]
        if row[2] != "cds":
            self.alScore = float(extractAttribute(row, "al_score"))
        if row[2] == "intron":
            self.spliceSites = extractAttribute(row, "splice_sites")
        row[8] = ""

    def add(self, row):
        if row[2] != "cds":
            self.alScore = max(float(extractAttribute(row, "al_score")),
                               self.alScore)
        self.prots.append(extractAttribute(row, "prot"))
        self.count += 1

    def print(self, printProts):
        self.row[5] = str(self.count)
        if self.row[2].lower() != "cds":
            if self.alScore == 0:
                self.alScore = "0"
            self.row[8] += f'al_score={self.alScore};'
        else:
            self.row[2] = "CDS"

        if self.row[2] == "intron":
            self.row[8] += f' splice_sites={self.spliceSites};'

        if printProts:
            self.row[8] += " prots=" + ",".join(self.prots) + ';'

        if self.row[8] == "":
            self.row[8] = "."

        return "\t".join(self.row)


class Codon():
    def __init__(self, arg):
        self.arg = arg


def loadData(inputFile):
    features = {}
    for row in csv.reader(open(inputFile), delimiter='\t'):
        if len(row) != 9:
            continue

        row[2] = row[2].lower()
        if row[2] != "intron" and row[2] != "start_codon" and \
           row[2] != "cds" and row[2] != "start" and \
           row[2] != "stop_codon" and row[2] != "stop":
            continue

        if signature(row) not in features:
            features[signature(row)] = Feature(row)
        else:
            features[signature(row)].add(row)

    return features


def printCollapsed(features, printProts, outputFile=None, append=False):
    if outputFile:
        if append:
            output = open(outputFile, "a")
        else:
            output = open(outputFile, "w")

    for f in features.values():
        if outputFile:
            output.write(f.print(printProts) + "\n")
        else:
            print(f.print(printProts))
    if outputFile:
        output.close()


def collapse(inputFile, printProts=True, outputFile=None, append=False):
    features = loadData(inputFile)
    printCollapsed(features, printProts, outputFile, append)


def main():
    args = parseCmd()
    collapse(args.input, not args.dontPrintProteins)


def parseCmd():

    parser = argparse.ArgumentParser(description='Collapse introns, starts,\
        and stops and compute their coverage/depth')

    parser.add_argument('input', metavar='input.gff', type=str,
                        help='The input gff to collapse.')

    parser.add_argument('--dontPrintProteins', action='store_true',
                        help='Do not print source proteins for each feature.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
