#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# ==============================================================


import argparse
import sys
import time
import subprocess
import collapseGff
import tempfile
import shutil
import os


tempFiles = []
workDir = ''
binDir = ''


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('[' + time.ctime() + '] error: exited due to an ' +
                 'error in command: ' + cmd)


def callScript(name, args):
    systemCall(binDir + '/' + name + ' ' + args)


def temp(prefix, suffix):
    if not os.path.isdir(workDir + "/tmp"):
        os.mkdir(workDir + "/tmp")
    tmp = tempfile.NamedTemporaryFile(delete=False, dir=workDir + "/tmp",
                                      prefix=prefix, suffix=suffix)
    tempFiles.append(tmp.name)
    return tmp


def cleanup():
    for file in tempFiles:
        os.remove(file)
    shutil.rmtree(workDir + "/tmp")


def setup(args):
    global workDir
    global binDir
    workDir = args.workdir
    if not os.path.isdir(workDir):
        os.mkdir(workDir)
    binDir = os.path.abspath(os.path.dirname(__file__))


def processMiniprotOutput(miniprot):
    processIntrons(miniprot)
    processStarts(miniprot)
    processStops(miniprot)
    callScript('print_high_confidence.py',
               f'{workDir}/miniprothint.gff > {workDir}/hc.gff')


def processIntrons(miniprot):
    intronsAll = temp('intronsAll', '.gff')
    systemCall(f'grep intron {miniprot} > {intronsAll.name}')
    introns01 = temp('introns01', '.gff')
    callScript('print_high_confidence.py',
               f'{intronsAll.name} --intronCoverage 0 --intronAlignment 0.1 '
               f'--addAllSpliceSites > {introns01.name}')
    collapseGff.collapse(introns01.name, outputFile=f'{workDir}/miniprothint.gff')


def processStops(miniprot):
    stopsAll = temp('stopsAllEnd', '.gff')
    systemCall(f'grep stop_codon {miniprot} | grep proteinEnd=1 > {stopsAll.name}')
    stopsPositive = temp('stopsPositive', '.gff')
    callScript('print_high_confidence.py',
               f'{stopsAll.name} --stopCoverage 0 --stopAlignment 0.01 '
               f'> {stopsPositive.name}')
    collapseGff.collapse(stopsPositive.name,
                         outputFile=f'{workDir}/miniprothint.gff',
                         append=True)


def processStarts(miniprot):
    startsAll = temp('startsAll', '.gff')
    systemCall(f'grep start_codon {miniprot} > {startsAll.name}')
    startsPositive = temp('startsPositive', '.gff')
    callScript('print_high_confidence.py',
               f'{startsAll.name} --startCoverage 0 --startAlignment 0.01 '
               f'> {startsPositive.name}')
    startsCollapsed = temp('startsCollapsed', '.gff')
    startsCollapsedS = temp('startsCollapsedSorted', '.gff')
    collapseGff.collapse(startsPositive.name, outputFile=startsCollapsed.name)
    systemCall(f'sort -k1,1 -k4,4n -k5,5n {startsCollapsed.name} > '
               f'{startsCollapsedS.name}')

    cds = temp('cds', '.gff')
    cdsC = temp('cdsCollapsed', '.gff')
    cdsSupported = temp('cdsCollapsed', '.gff')
    cdsSupportedS = temp('cdsCollapsed', '.gff')
    systemCall(f'grep CDS {miniprot} > {cds.name}')
    collapseGff.collapse(cds.name, printProts=False, outputFile=cdsC.name)

    # This is crucial as there is so much noise in the CDS alignments.
    # Without this step, almost no starts are left with a larger database.
    callScript('cds_with_upstream_support.py',
               f'{cdsC.name} {startsCollapsedS.name} '
               f'{workDir}/miniprothint.gff > {cdsSupported.name}')
    systemCall(f'sort -k1,1 -k4,4n -k5,5n {cdsSupported.name} > '
               f'{cdsSupportedS.name}')

    callScript('count_cds_overlaps.py',
               f'{startsCollapsedS.name} {cdsSupportedS.name} >> '
               f'{workDir}/miniprothint.gff')


def main():
    args = parseCmd()
    setup(args)
    processMiniprotOutput(args.miniprot)
    if (not args.nocleanup):
        cleanup()


def parseCmd():

    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    parser.add_argument('miniprot', metavar='miniprot_scored.gff', type=str,
                        help='Miniprot output scored by the miniprot\
                              boundary scorer.')

    parser.add_argument('--workdir', type=str, default='.',
                        help='Keep all the temporary files.')

    parser.add_argument('--nocleanup', action='store_true',
                        help='Keep all the temporary files.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
