#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Select the best alignments from overlapping miniprot results.
# ==============================================================


import argparse
import re
import sys
import os
import csv
import AlignmentCluster


class Exon():
    def __init__(self, row, parent):
        self.contig = row[0]
        self.start = int(row[3])
        self.end = int(row[4])
        self.strand = row[6]
        self.parent = parent

    def print(self, Parent):
        print("\t".join([self.contig, "miniprot", "CDS", str(self.start),
                         str(self.end), ".", self.strand, ".",
                         f'Parent={Parent};']))

    def printGff2(self, match):
        print("\t".join([self.contig, "miniprot", "HSP", str(self.start),
                         str(self.end), ".", self.strand, ".",
                         f'Match {match}']))

    def __lt__(self, other):
        if self.contig != other.contig:
            return self.contig < other.contig
        elif self.strand != other.strand:
            return self.strand < other.strand
        elif self.start != other.start:
            return self.start < other.start
        else:
            return self.end < other.end


class Alignment():
    def __init__(self, ID, coverage, target):
        self.ID = ID
        self.start = sys.maxsize
        self.end = -sys.maxsize
        self.cluster = None
        self.CDSlen = 0
        self.exons = []
        self.used = False
        self.selected = False
        self.selectedCount = 0
        self.seed = False
        self.subLocus = False
        self.coverage = coverage
        self.identity = 0
        self.target = target.split()[0]

    def addExon(self, exon):
        self.start = min(self.start, exon.start)
        self.end = max(self.end, exon.end)
        self.CDSlen += exon.end - exon.start + 1
        self.exons.append(exon)

    def addScore(self, score):
        self.score = float(score)

    def getContig(self):
        return self.exons[0].contig

    def getStrand(self):
        return self.exons[0].strand

    def print(self, clusterId='', exons=True):
        print("\t".join([self.getContig(), "miniprot", "mRNA", str(self.start),
                         str(self.end), str(self.score), self.getStrand(),
                         ".", f'ID={self.ID}; score={str(self.score)};'
                         f'qcov={str(round(self.coverage, 4))};'
                         f'identity={str(round(self.identity, 4))};'
                         f'cluster={str(clusterId)}']))
        if exons:
            for exon in self.exons:
                exon.print(self.ID)

    def printGff2(self):
        print("\t".join([self.getContig(), "miniprot", "match",
                         str(self.start), str(self.end), str(int(self.score)),
                         self.getStrand(), ".",
                         f'Match {self.target}_{self.ID};'
                         f'coverage {str(round(self.coverage, 2))};'
                         f'subjectName {self.target};'
                         f'Alias {self.target}'
                         ]))
        for exon in self.exons:
            exon.printGff2(f'{self.target}_{self.ID}')

    def printLocus(self, outFh, addQcov=False):
        toPrint = [self.getContig(), self.getStrand(), str(self.start),
                   str(self.end), self.target]
        if addQcov:
            toPrint.append(str(round(self.coverage, 2)))
        outFh.write("\t".join(toPrint) + "\n")

    def getCDSOverlap(self, other):
        if self.start > other.end or self.end < other.start:
            return 0
        i = 0
        j = 0
        overlap = 0
        while i < len(self.exons) and j < len(other.exons):
            exon1 = self.exons[i]
            exon2 = other.exons[j]
            if exon1.end < exon2.start:
                i += 1
            elif exon2.end < exon1.start:
                j += 1
            else:
                overlap += min(exon1.end, exon2.end) - \
                    max(exon1.start, exon2.start) + 1
                if exon1.end < exon2.end:
                    i += 1
                else:
                    j += 1
        return overlap / self.CDSlen

    def __lt__(self, other):
        return self.score > other.score


def extractAttributeGtf(text, att):
    regex = att + ' "([^"]+)"'
    result = re.search(regex, text)
    if result:
        return result.groups()[0]
    else:
        return None


def extractAttributeGff(text, att):
    regex = att + '=([^;]+)'
    result = re.search(regex, text)
    if result:
        return result.groups()[0]
    else:
        return None


def getRootCluster(clusterTree, i):
    """Return the root cluster of a given cluster.

    Args:
        clusterTree: Cluster pointers
        i: Index of the cluster of interest
    """
    while clusterTree[i] != i:
        i = clusterTree[i]
    return i


def clusterAlignments(allExons, alignments):
    """Cluster overlapping seeds. Only CDS-level overlaps are considered.
    """

    clusterId = -1
    prevContig = ""
    prevStrand = ""
    currentClusterEnd = 0
    clusterTree = []
    ID2cluster = {}
    clusters = {}

    # Cluster by exons. Do not build objects for the temp clusters
    for exon in allExons:
        if prevContig != exon.contig or exon.start > currentClusterEnd or \
           prevStrand != exon.strand:
            clusterId += 1
            clusterTree.append(clusterId)
            currentClusterEnd = exon.end
        else:
            if exon.end > currentClusterEnd:
                currentClusterEnd = exon.end

        if exon.parent not in ID2cluster:
            ID2cluster[exon.parent] = clusterId
        elif ID2cluster[exon.parent] != clusterId:
            parentCluster = getRootCluster(clusterTree, clusterId)
            clusterTree[parentCluster] = getRootCluster(clusterTree,
                                                        ID2cluster[exon.parent]
                                                        )
        prevContig = exon.contig
        prevStrand = exon.strand

    # Assign alignments to the final clusters
    for alignment in alignments.values():
        rootClusterId = getRootCluster(clusterTree, ID2cluster[alignment.ID])
        if rootClusterId not in clusters:
            clusters[rootClusterId] = AlignmentCluster.Cluster(rootClusterId)
        clusters[rootClusterId].addAlignment(alignment)

    return clusters


def getCoverageFromPAF(row):
    # These are bed-like coordinates
    return (int(row[4]) - int(row[3])) / int(row[2])


def loadGff(miniprot):
    alignments = {}
    allExons = []
    pafDetected = False
    coverage = -1
    for row in csv.reader(open(miniprot), delimiter='\t'):
        if row[0][0] == "#":
            if row[0] == "##PAF":
                pafDetected = True
                coverage = getCoverageFromPAF(row)
            else:
                continue

        if row[2] == 'mRNA':
            ID = extractAttributeGff(row[8], "ID")
        elif row[2] == 'CDS':
            ID = extractAttributeGff(row[8], "Parent")
        else:
            continue

        if ID not in alignments:
            if pafDetected:
                target = extractAttributeGff(row[8], "Target")
            else:
                target = extractAttributeGff(row[8], "prot")
            alignments[ID] = Alignment(ID, coverage, target)

        if row[2] == 'mRNA':
            alignments[ID].addScore(row[5])
            alignments[ID].identity = float(extractAttributeGff(row[8],
                                                                "Identity"))
            if not pafDetected:
                alignments[ID].coverage = float(extractAttributeGff(row[8],
                                                                    "qcov"))
        else:
            exon = Exon(row, ID)
            alignments[ID].addExon(exon)
            allExons.append(exon)

    return allExons, alignments


def loadAlignments(miniprot):

    ext = os.path.splitext(miniprot)[1]

    if ext == ".gff" or ext == ".gff3":
        allExons, alignments = loadGff(miniprot)
    else:
        sys.exit(f'error: Unexpected file extension: {ext}')

    allExons.sort()
    # We could sort this later when cycling through alignments, but sorting
    # here makes the code more predictable (worth being a bit slower)
    for alignment in alignments.values():
        alignment.exons.sort()

    return allExons, alignments


def printSelected(miniprot, selected):
    selected = set(selected)
    for row in csv.reader(open(miniprot), delimiter='\t'):
        if row[0][0] == "#":
            continue

        if row[2] == 'mRNA':
            ID = extractAttributeGff(row[8], "ID")
        else:
            ID = extractAttributeGff(row[8], "Parent")

        if ID in selected:
            print("\t".join(row))


def main():
    args = parseCmd()
    allExons, alignments = loadAlignments(args.miniprot)

    clusters = clusterAlignments(allExons, alignments)

    selected = []
    for cluster in clusters.values():
        s = cluster.splitByBestAlignments(args.minOverlap4SeedChildren,
                                          args.minScoreFraction,
                                          args.minSeedCoverage,
                                          args.topNperSeed,
                                          args.maxSubFraction,
                                          args.minSubCoverage)
        selected += s

    printSelected(args.miniprot, selected)


def parseCmd():

    parser = argparse.ArgumentParser(description='Select the best\
        alignments from overlapping miniprot results.',
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    parser.add_argument('miniprot', metavar='miniprot.gtf/gff', type=str,
                        help='Raw miniprot alignments in a gff format \
        produced by miniprot or miniprot boundary scorer. If in native \
        miniprot format, each alignment must be preceded by the PAF line')

    parser.add_argument('--minSeedCoverage', type=float, default=0,
                        help='Minimum query coveragy for an alignment to be\
                        selected as a locus seed.')

    parser.add_argument('--minOverlap4SeedChildren', type=float, default=0.01,
                        help='Min CDS overlap of an alignment by a locus seed\
                        to be considered the seed\'s child')

    parser.add_argument('--minScoreFraction', type=float, default=0.9,
                        help='Print seed children in a locus only if their \
                        score is > seed.score * minScoreFraction')

    parser.add_argument('--topNperSeed', type=int, default=10,
                        help='Print a maximum of topNperSeed seed children per\
                        locus.')

    parser.add_argument('--maxSubFraction', type=float, default=0.8,
                        help='Maximum CDS coverage of the parent seed by an \
                        alignment in order for the alignment to be considered\
                        for spawning a subseed.')

    parser.add_argument('--minSubCoverage', type=float, default=0.9,
                        help='Spawn a subseed only if its alignment query \
                        coverage is >= minSubCovoverage. Further, the subseed\
                        needs to have has better average alignment identity \
                        than the parent.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
