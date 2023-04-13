#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Visualize miniprothint intron scores
# ==============================================================


import argparse
import csv
import re
import sys
import random
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def getSignature(row):
    return row[0] + "_" + row[3] + "_" + row[4] + "_" + row[6]


def extractFeature(text, feature):
    regex = feature + '=([^;]+)'
    search = re.search(regex, text)
    if search:
        return search.groups()[0]
    else:
        return None


def loadAnnotation(annotFile):
    annot = set()
    for row in csv.reader(open(annotFile), delimiter='\t'):
        if (row[2].lower() == "intron"):
            annot.add(getSignature(row))
    return annot


def plotScores(annot, inputFile, outputFile, args):
    TP = 0
    FP = 0
    allX = []
    allY = []
    colors = []
    maxX = 0
    maxY = 0

    for row in csv.reader(open(inputFile), delimiter='\t'):
        if (row[2].lower() != "intron"):
            continue

        signature = getSignature(row)
        color = 'purple'
        if signature in annot:
            TP += 1
            color = 'green'
        else:
            FP += 1

        x = float(extractFeature(row[8], "al_score")) + \
            random.uniform(-0.01, 0.01)
        y = float(row[5]) + random.uniform(-0.5, 0.5)

        if args.ylim != -1:
            if y > args.ylim:
                continue

        if x > maxX:
            maxX = x

        if y > maxY:
            maxY = y

        allX.append(x)
        allY.append(y)
        colors.append(color)

    plt.scatter(x=allX, y=allY, marker='o', s=0.1, color=colors,
                alpha=args.opacity, clip_on=False)

    lstyle = '--'
    lsize = 2
    plt.plot([0.25, maxX], [3.5, 3.5], color='tab:red', ls=lstyle,
             linewidth=lsize)
    plt.plot([0.25, 0.25], [3.5, maxY], color='tab:red', ls=lstyle,
             linewidth=lsize)
    plt.plot([0.1, 0.1], [0, maxY], color='b', ls=':', linewidth=lsize)

    # Legend
    yMargin = 0
    yShift = -0.048
    plt.text(0.05, 1 + yMargin, "   TP (" + str(TP) + ")" +
             "\n   FP (" + str(FP) + ")",
             va="top", transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', edgecolor='black',
             boxstyle='square,pad=0.75'))
    plt.plot(0.05, .989, '.', ms=10, color='green', clip_on=False,
             transform=plt.gca().transAxes, zorder=4)
    plt.plot(0.05, .989 + yShift, '.', ms=10, color='purple', clip_on=False,
             transform=plt.gca().transAxes, zorder=4)

    plt.xlabel("Intron borders alignment (IBA)")
    plt.ylabel("Intron mapping coverage (IMC)")
    if args.logYScale:
        plt.yscale('log')
    plt.box(on=None)
    if args.ylim != -1:
        plt.ylim(0, args.ylim + 0.5)
    plt.xlim(-0.01, maxX)
    plt.savefig(outputFile)


def main():
    args = parseCmd()
    with open("scatter.sh", "w") as file:
        file.write("#!/usr/bin/env bash\n")
        file.write(" ".join(sys.argv) + "\n")
    annot = loadAnnotation(args.annotation)
    plotScores(annot, args.input, args.output, args)


def parseCmd():

    parser = argparse.ArgumentParser(description='Visualize miniprothint\
                                                  intron scores')

    parser.add_argument('input', metavar='miniprothint.gff', type=str,
                        help='Introns scored by miniprot boundary \
                              scorer/miniprohtint')

    parser.add_argument('annotation', metavar='introns_annot.gff', type=str,
                        help='Annotated introns.')

    parser.add_argument('output', type=str,
                        help='Output figure.')

    parser.add_argument('--opacity', type=float, default=1,
                        help='Opacity of individual points. Default = 1')

    parser.add_argument('--ylim', type=float, default=-1,
                        help='Limit y axis by this number. Default = No limit')

    parser.add_argument('--logYScale',  default=False, action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()
