#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 12:24:32 2018

@author: antony
"""

import sys
sys.path.append('/ifs/scratch/cancer/Lab_RDF/ngs/scripts/python/lib/libdna')
sys.path.append('/ifs/scratch/cancer/Lab_RDF/ngs/scripts/python/lib/libbam')
import argparse
import libdna
import libbam
import re

SAMTOOLS='/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-1.8/bin/samtools'
#SAMTOOLS='/ifs/scratch/cancer/Lab_RDF/ngs/tools/samtools-0.1.19/samtools'
DNA_DIR='/ifs/scratch/cancer/Lab_RDF/ngs/genomes/ucsc/2bitext/hg19'

#
# Print header
#

def merge(bam, chr, out='', max_len=1000, samtools=SAMTOOLS, dna_dir=DNA_DIR):
    read_map = {} #'collections.defaultdict(pysam.AlignedSegment)

    sam = libbam.BamReader(bam, samtools=samtools)
    
    if out == '':
        out = re.sub(r'\.bam','.merged.bam', bam)
    
    print('Writing to {} using {} with max length {}...'.format(out, samtools, max_len))
    
    
    out_bam = libbam.BamWriter(out, samtools=samtools)
    
    out_bam.write_header(sam)
    
    
    #out = pysam.AlignmentFile('merged.sam', 'w', template=samfile)
    
    c = 1
    
    pc = 0
    pc_ins = 0
    pc_trunc = 0
    
    dna = libdna.CachedDNA2Bit(dna_dir)
    
    for read in sam.reads(chr):
        #for read in samfile:
            
        if not read.is_paired or not read.is_proper_pair:
            continue
        
        if read.qname in read_map:
            # found a pair so check the distance
            
            r1 = read_map[read.qname]
            r2 = read
            
            if r1.pos > r2.pos:
                # r1 must start before r2, if not, swap them
                r3 = r2
                r2 = r1
                r1 = r3
            
            # merge the sequences together
            seq = dna.merge_read_pair_seq(r1, r2)
            
            if len(seq) > len(r1.seq) + len(r2.seq):
                pc_ins += 1
            else:
                pc_trunc += 1
        
            # remove from map so we don't waste memory storing already
            # processed pairs
            del read_map[r1.qname]
            
            # Use a default flag of 0 to say the merged read is on
            # the forward strand, is unpaired and maps
            flag = 0
            
            # make a new read that is the merge of the two
            read = libbam.SamRead(r1.qname,
                 flag,
                 r1.rname,
                 r1.pos,
                 r1.mapq,
                 '{}M'.format(len(seq)),
                 '*',
                 0,
                 0,
                 seq,
                 '*',
                 tags=['NM:i:0'])
            
            
            #print('\t'.join([r1.qname, str(r1.flag), r1.chr, str(r1.pos),  str(r1.mapq), '{}M'.format(len(seq)), '*', '0', '0', seq, '*', 'NM:i:0']))
            
            # Write to bam file
            
            if read.length <= max_len:
                out_bam.write(read)
        
            pc += 1
        else:
            # Store until we find its partner
            read_map[read.qname] = read
        
        if c % 100000 == 0:
            print('Processed {} reads, inserts={}, truncations={}, cache={}...'.format(c, pc_ins, pc_trunc, len(read_map)))
      
        c += 1
    
    out_bam.close()
    #samfile.close()    
    #out.close()




if __name__ == "__main__":
    #bam = '/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/data/samples/hg19/rdf/elodie/LY1_H3K27ACAM_EP001/genome/hg19/hisat2/LY1_H3K27ACAM_EP001_hg19.fixmate.sorted.markdup.bam' #sys.argv[1]
    #bam = '/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/data/samples/hg19/rdf/elodie/LY1_Input_EP003/genome/hg19/hisat2/LY1_Input_EP003_hg19.fixmate.sorted.markdup.bam' #sys.argv[1]
    #bam = '/ifs/scratch/cancer/Lab_RDF/ngs/ChIP_seq/data/samples/hg19/rdf/elodie/LY1_H3K27ACDiag_EP002/genome/hg19/hisat2/LY1_H3K27ACDiag_EP002_hg19.fixmate.sorted.markdup.bam'
    #bam = sys.argv[1]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs=1)
    parser.add_argument('chr', nargs=1)
    parser.add_argument('--samtools', nargs='?', default=SAMTOOLS)
    parser.add_argument('--out', nargs='?', default='')
    parser.add_argument('--max_len', nargs='?', type=int, default=1000)
    args = parser.parse_args()
    
    merge(args.file[0], args.chr[0], out=args.out, max_len=args.max_len, samtools=args.samtools)
    
