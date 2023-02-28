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
workDir = "."


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('[' + time.ctime() + '] error: exited due to an ' +
                 'error in command: ' + cmd)


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
    workDir = args.workdir
    if not os.path.isdir(workDir):
        os.mkdir(workDir)


def processMiniprotoutput(miniprot):
    collapsed = temp("collapsed", ".gff")
    collapseGff.collapse(miniprot, outputFile=collapsed.name)


def main():
    args = parseCmd()
    setup(args)
    processMiniprotoutput(args.miniprot)
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
