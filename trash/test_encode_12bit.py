# -*- coding: utf-8 -*-
"""
Encode read counts per base in 2 bytes

@author: Antony Holmes
"""

import sys

import lib.encode

counts = [0, 4000, 256]
bytes = lib.encode.encode_counts(counts, "chr1", 40, 12)



sys.stderr.write(' '.join('{:08b}'.format(x) for x in bytes) + "\n")

ret = lib.encode.decode_12bit(bytes[1:len(bytes)], 1, 2)

sys.stderr.write("0 offset " + ' '.join([str(x) for x in ret]) + "\n")
