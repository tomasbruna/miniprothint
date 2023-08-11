#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Class for representing alignment clusters
#
# ==============================================================


import sys
import math


class Border():

    def __init__(self, coordinate, start):
        self.coordinate = coordinate
        self.start = start

    def print(self):
        print(self.coordinate, self.start)

    def __lt__(self, other):
        return self.coordinate < other.coordinate


class Block():

    def __init__(self, start, end, coverage):
        self.start = start
        self.end = end
        self.coverage = coverage
        self.state = ""

    def print(self, contig="-"):
        print("\t".join([contig, str(self.start), str(self.end),
                         f'Coverage={self.coverage}']))


class Cluster():
    def __init__(self, clusterId):
        self.clusterId = clusterId
        self.start = sys.maxsize
        self.end = -sys.maxsize
        self.alignments = []
        self.borders = []

    def addAlignment(self, alignment):
        self.start = min(self.start, alignment.start)
        self.end = max(self.end, alignment.end)
        self.alignments.append(alignment)

    def getContig(self):
        return self.alignments[0].getContig()

    def print(self):
        print("\t".join([self.getContig(), str(self.start), str(self.end),
                         f'ID={self.clusterId}']))

    def printAlignments(self):
        for alignment in self.alignments:
            alignment.print(clusterId=self.clusterId)

    def printSubClusters(self):
        for subCluster in self.subClusters:
            print("\t".join([self.getContig(), str(subCluster.start),
                             str(subCluster.end),
                             f'ID={self.clusterId}_{subCluster.clusterId}']))

    def defineBorders(self):
        self.borders = []
        self.CDSborders = []
        for al in self.alignments:
            self.borders.append(Border(al.start, True))
            self.borders.append(Border(al.end, False))
            for exon in al.exons:
                self.CDSborders.append(Border(exon.start, True))
                self.CDSborders.append(Border(exon.end, False))
        self.borders.sort()
        self.CDSborders.sort()

    def computeMeanCDSCoverage(self):
        coverage = 0
        prevCoordinate = self.CDSborders[0].coordinate
        meanCoverage = 0
        CDSArea = 0

        for border in self.CDSborders:
            if border.coordinate != prevCoordinate:
                if coverage != 0:
                    CDSArea += border.coordinate - prevCoordinate
                    meanCoverage += (border.coordinate - prevCoordinate) * \
                        coverage

            if border.start:
                coverage += 1
            else:
                coverage -= 1

            prevCoordinate = border.coordinate

        self.meanCDSCoverage = meanCoverage / CDSArea

    def computeCoverage(self):
        coverage = 0
        self.blocks = []
        self.maxCoverage = 0
        self.meanCoverage = 0
        prevCoordinate = self.borders[0].coordinate
        blockLen = self.borders[-1].coordinate - self.borders[0].coordinate

        for border in self.borders:
            if border.coordinate != prevCoordinate:
                self.meanCoverage += ((border.coordinate - prevCoordinate)
                                      / blockLen) * coverage
                if len(self.blocks) != 0 and \
                   self.blocks[-1].coverage == coverage:
                    self.blocks[-1].end = border.coordinate
                else:
                    self.blocks.append(Block(prevCoordinate, border.coordinate,
                                             coverage))

                    if coverage > self.maxCoverage:
                        self.maxCoverage = coverage

            if border.start:
                coverage += 1
            else:
                coverage -= 1

            prevCoordinate = border.coordinate

    def detectBridges(self, baseline, bridgeEnterT, bridgeExitT):
        state = "start"
        bridgeEnter = int(bridgeEnterT * baseline)
        bridgeExit = math.ceil(bridgeExitT * baseline)
        self.subClusters = []

        for block in self.blocks:
            if state == "start":
                # Always start the first subcluster from the beginning because
                # a bridge needs to connect two subclusters
                self.subClusters.append(Cluster(len(self.subClusters) + 1))
                self.subClusters[-1].start = block.start
                if block.coverage <= bridgeEnter:
                    state = "edge"
                else:
                    state = "subCluster"
            elif state == "subCluster":
                if block.coverage <= bridgeEnter:
                    state = "bridge"
                    self.subClusters[-1].end = block.start - 1
            elif state == "bridge":
                if block.coverage >= bridgeExit:
                    state = "subCluster"
                    self.subClusters.append(Cluster(len(self.subClusters) + 1))
                    self.subClusters[-1].start = block.start
            elif state == "edge":
                if block.coverage >= bridgeExit:
                    state = "subCluster"

            block.state = state

        # Close the last subcluster. If a bridge was open at the end, make it
        # a part of the subcluster. This mirrors the opening behavior.
        self.subClusters[-1].end = self.blocks[-1].end

    def splitBridges(self, bridgeEnterT, bridgeExitT):
        self.defineBorders()
        self.computeCoverage()
        self.computeMeanCDSCoverage()
        self.detectBridges(self.meanCDSCoverage, bridgeEnterT, bridgeExitT)
        self.printSubClusters()

    def getNextSeed(self, minSeedCov):
        while self.lastSeed < len(self.alignments):
            if self.alignments[self.lastSeed].used:
                self.lastSeed += 1
            else:
                if self.alignments[self.lastSeed].coverage < minSeedCov:
                    self.alignments[self.lastSeed].used = True
                    self.lastSeed += 1
                else:
                    self.alignments[self.lastSeed].seed = True
                    self.alignments[self.lastSeed].selected = True
                    self.alignments[self.lastSeed].used = True
                    return self.alignments[self.lastSeed]
        return None

    def processOverlappingAlignmens(self, seed, minOverlap, minScoreFraction,
                                    topNperSeed, maxSubFrac, minSubCov):
        # maxSubFrac is the maximum coverage of the parent seed by a
        # potential subseed in order for subseed to be used. Setting
        # to 0 turns subseeed off completely.
        # Further, subseed is only selected when it has better average
        # alignment identity than the parent and query coverage >= minSubCov.
        for i in range(self.lastSeed + 1, len(self.alignments)):
            alignment = self.alignments[i]
            if not alignment.used:
                if alignment.getCDSOverlap(seed) > minOverlap:
                    if alignment.score > seed.score * minScoreFraction and \
                            seed.selectedCount < topNperSeed:
                        alignment.selected = True
                        alignment.used = True
                        seed.selectedCount += 1
                    elif seed.getCDSOverlap(alignment) <= maxSubFrac and \
                            alignment.coverage >= minSubCov and \
                            alignment.identity >= seed.identity:
                        # If not selected as a seed rep., give it a chance
                        # to be selected as a sublocus seed later.
                        # A potential issue to consider: All sublocus reps.
                        # need to satisfy the sublocus seed conditions.
                        # Could be resolved by multiple iterations for each
                        # subseed if this really turns out to matter.
                        alignment.subLocus = True
                    else:
                        alignment.used = True

    def splitByBestAlignments(self, minOverlap=0.01,
                              minScoreFraction=0.90, minSeedCov=0,
                              topNperSeed=10, maxSubFrac=0.8, minSubCov=0.9):
        self.alignments.sort()
        self.lastSeed = 0

        seed = self.getNextSeed(minSeedCov)
        while seed is not None:
            self.processOverlappingAlignmens(seed, minOverlap,
                                             minScoreFraction, topNperSeed,
                                             maxSubFrac, minSubCov)
            seed = self.getNextSeed(minSeedCov)

        selected = []
        for alignment in self.alignments:
            if alignment.selected:
                selected.append(alignment.ID)
        return selected
