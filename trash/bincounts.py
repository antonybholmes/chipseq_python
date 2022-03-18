#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:56:33 2018

@author: antony
"""

import sys
import libseq

SAMTOOLS='/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin/samtools'

bam = sys.argv[1]
genome = sys.argv[2]
power = int(sys.argv[3])
mode = sys.argv[4]

bc = libseq.BinCountWriter(bam, genome, power, mode=mode, samtools=SAMTOOLS)
bc.write_all()
