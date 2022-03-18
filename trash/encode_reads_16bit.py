# -*- coding: utf-8 -*-
"""
Encode read counts per base in 2 bytes

@author: Antony Holmes
"""

import sys

import lib.encode

chr_size_file = sys.argv[1]
file = sys.argv[2]
chr = sys.argv[3]
read_length = int(sys.argv[4])
window = int(sys.argv[5])

#lib.encode.encode_sam_16bit(chr_size_file, file, chr, read_length, window)
lib.encode.encode_sam(chr_size_file, file, chr, read_length, window)
