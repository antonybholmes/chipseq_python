# -*- coding: utf-8 -*-
"""
Annotate chipseeqer peaks

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import re
from typing import List
import numpy as np
import os
import argparse


import pychipseq.genomic
import pychipseq.human.peaks


parser = argparse.ArgumentParser(prog='annotate_chipseeqer_peaks.py')
parser.add_argument('-m', '--matches', nargs=1)
parser.add_argument('-g', '--genome', nargs=1, default=['hg19'])
parser.add_argument('--prom5', nargs=1, default=[2000])
parser.add_argument('--prom3', nargs=1, default=[1000])
parser.add_argument('-k', '--skip', nargs=1, default=[False])
args = parser.parse_args()

print(args)

matches = set(args.matches[0].split(','))
genome = args.genome[0]
prom_ext_5p = int(args.prom5[0])
prom_ext_3p = int(args.prom3[0])
skipIfExists = args.skip[0]

# exit()

#file = sys.argv[1]
#out = sys.argv[2]
#prom_ext_5p = int(sys.argv[3])
#prom_ext_3p = int(sys.argv[4])

print(f'----------', file=sys.stderr)
print(f'Annotation', file=sys.stderr)
print(f'----------', file=sys.stderr)

peak_annotation = pychipseq.human.peaks.AnnotatePeak(
    "Peak", prom_ext_5p, prom_ext_3p)


for (dir, dirs, files) in os.walk('.'):
    found = False

    for m in matches:
        if m in dir:
            found = True
            break

    if not found:
        continue

    for file in files:
        if 'TF_target' in file:
            out = file.replace('TF_targets', 'Peaks').replace('txt', 'tsv').replace(
                '.tsv', f'_5p{prom_ext_5p}_3p{prom_ext_3p}_{genome}.tsv')

            print(file=sys.stderr)
            print(f'dir: {dir}', file=sys.stderr)
            print(f'in: {file}', file=sys.stderr)
            print(f'out: {out}', file=sys.stderr)
            file = f'{dir}/{file}'
            out = f'{dir}/{out}'

            if skipIfExists and os.path.exists(out):
                print(f'file exists, skipping', file=sys.stderr)
            else:
                print(f'working...', file=sys.stderr)
                df = peak_annotation.parse(file)
                df.to_csv(out, sep='\t', header=True, index=False)

            print(f'finished', file=sys.stderr)
            

    # pass

print(f'----', file=sys.stderr)

exit(0)

print("Annotating " + file, file=sys.stderr)

df = peak_annotation.parse(file)

df.to_csv(out, sep='\t', header=True, index=False)

exit(0)

#peaksPerGene = pychipseq.peaks.PeaksPerGeneAnnotation()


#f = open(file, "r")

# skip header if not a chipseeqer file
# if "TF_targets" not in file:
#    f.readline()

# header
#header = []
# header.append(LOCATION_HEADING)
# header.append(PVALUE_HEADING)
##sys.stdout.write("\tMax Height (reads)")
# header.append(SCORE_HEADING)
# header.append(WIDTH_HEADING)

# header.extend(peak_annotation.get_header())

# header.extend(peaksPerGene.get_names())

#print(header, file=sys.stderr)
# print("\t".join(header))

#locations = []
#annotations = []

for line in f:
    tokens = line.strip().split("\t")

    if pychipseq.genomic.is_location(tokens[0]):
        location = pychipseq.genomic.parse_location(tokens[0])
    else:
        location = pychipseq.genomic.Location(
            tokens[0], int(tokens[1]), int(tokens[2]))

    #chr = tokens[0]

    # Skip chrM
    if "chrM" in location.chr:
        continue

    locations.append(location)

    #start = int(tokens[1])
    #end = int(tokens[2])
    width = location.end - location.start + 1

    if "TF_targets" in file:
        # deal with infs
        if re.match(r'.*[Ii]nf.*', tokens[3]):
            p = -1000
        else:
            p = float(tokens[3])

        score = float(tokens[4])
    else:
        p = 0
        score = 0

    #max_height = int(tokens[6])

    # location = lib.genomic.Location(chr, start, end) #chr + ":" + str(start) + "-" + str(end)

    # Write the initial portion
    annotation = [str(location), str(p), str(score), str(width)]

    row_map = {LOCATION_HEADING: str(
        location), PVALUE_HEADING: p, SCORE_HEADING: score, WIDTH_HEADING: width}

    # Add the gene annotation
    annotation.extend(peak_annotation.annotate(location, row_map=row_map))

    annotations.append(annotation)

f.close()

peaksPerGene.annotate(header, annotations)

#
# Output table
#
for i in range(0, len(locations)):
    print("\t".join(annotations[i]))
