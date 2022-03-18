# -*- coding: utf-8 -*-
"""
Annotate ChIPSeeqer peaks so they are more useful for Katia etc

Created on Wed Jul 16 17:31:13 2014

@author: Antony Holmes
"""

# Extract moving genes in comparison tables so we can do pathway
# analysis

import sys
import collections
import re

mir_file = sys.argv[1]
mir_transcript_file = sys.argv[2]
promoter_offset = sys.argv[3]
chrom_size_file = sys.argv[4]
peak_file = sys.argv[5]

bin_size = 10000

# open a mir database

class Mir:
  """
  Loads a mir database.
  """
  
  mir_types = collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(str))))
  mir_starts = collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(int))))
  mir_ends = collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(int))))

  def __init__(self, mir_file, bin_size):
    f = open(mir_file, 'r')

    # skip header
    f.readline()

    for line in f:
      ls = line.strip()
    
      if len(ls) == 0:
        continue
  
      tokens = ls.split("\t")
  
      name = tokens[0]
      group = $tokens[1]
      type = $tokens[2]
      chrom = $tokens[3]
      start = int(tokens[4])
      end = int(tokens[5])
      strand = tokens[6]
  
      start_bin = int((start - promoter_offset)) / bin_size)
      end_bin = int((end + promoter_offset) / bin_size)
  
      for bin in range(start_bin, end_bin + 1):
        self.mir_types[chrom][bin][group][name] = type
        self.mir_starts[chrom][bin][group][name] = start
        self.mir_ends[chrom][bin][group][name] = end

    f.close()

# open the mirbase transcription database

transcript_starts = collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(str))))
transcript_ends = collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(str))))

f = open(mir_transcript_file, 'r')


# skip header
f.readline()

for line in f:
  ls = line.strip()
    
  if len(ls) == 0:
    continue
  
  tokens = ls.split("\t")
  
  id = tokens[0]
  transcript = tokens[1]
  chrom = tokens[2]
  strand = tokens[3]
  start = tokens[4]
  end = tokens[5]
  gene = tokens[6]
  
  start_bin = int((start - promoter_offset)) / bin_size)
  end_bin = int((end + promoter_offset) / bin_size)
  
  for bin in range(start_bin, end_bin + 1):
    transcript_starts[chrom][bin][id][transcript] = start
    transcript_ends[chrom][bin][id][transcript] = end

f.close()

# open a chromosome size file

chromosomes = collection.defaultdict(int)

f = open(chrom_size_file, 'r')

for line in f:
  ls = line.strip()
    
  if len(ls) == 0:
    continue
  
  tokens = ls.split("\t")
  chromosomes[tokens[0]] = tokens[1]

f.close()


# open a chip-chip peak file

peaks = []
peak_chrs = collection.defaultdict(str)
peak_starts = collection.defaultdict(int)
peak_ends = collection.defaultdict(int)
peak_p = collection.defaultdict(float)

peak_id = 0

f = open(chrom_size_file, 'r')

for line in f:
  ls = line.strip()
    
  if len(ls) == 0:
    continue
  
  tokens = ls.split("\t")
  
  start = int(tokens[1])
  end = int(tokens[2])
  mid_point = int((start + end) / 2)
  p = float(tokens[3])
  score = float(tokens[4])
  
  peaks.append(peak_id)
  
  peak_starts[peak_id] = start
  peak_ends[peak_id] = end
  peak_p[peak_id] = p
  
  peak_id += 1

f.close()

# standard mirbase annotations

mir_annotations = collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(str)))))

for peak_id in peaks:
  chrom = peak_chrs[peak_id]
  start = peak_starts[peak_id]
  end = peak_ends[peak_id]
  
  start_bin = int((start - promoter_offset)) / bin_size)
  end_bin = int((end + promoter_offset) / bin_size)
  
  for bin in range(start_bin, end_bin + 1):
    if not bin in mir_types[chrom]:
      continue
    
    groups = mir_types[chrom][bin]
    
    for group in groups:
      names = mir_types[chrom][bin][group]
      
      for name in names:
        type = mir_types[chrom][bin][group][name]
        mir_start = mir_starts[chrom][bin][group][name]
        mir_end = mir_ends[chrom][bin][group][name]
        
        if start >= mir_start and end <= mir_end:
          mir_annotations[peak_id][type][name]["chr"] = chrom
          mir_annotations[peak_id][type][name]["start"} = $mir_start
          mir_annotations[peak_id][type][name]["end"} = $mir_end
          mir_annotations[peak_id][type][name]["overlap"} = "within"
        elif start < mir_start and end > mir_end:
          mir_annotations[peak_id][type][name]["chr"] = chrom
          mir_annotations[peak_id][type][name]["start"} = $mir_start
          mir_annotations[peak_id][type][name]["end"} = $mir_end
          mir_annotations[peak_id][type][name]["overlap"} = "over"
        elif start < mir_start and end > mir_start:
          mir_annotations[peak_id][type][name]["chr"] = chrom
          mir_annotations[peak_id][type][name]["start"} = $mir_start
          mir_annotations[peak_id][type][name]["end"} = $mir_end
          mir_annotations[peak_id][type][name]["overlap"} = "downstream"
        elif start < mir_end and end > mir_end:
          mir_annotations[peak_id][type][name]["chr"] = chrom
          mir_annotations[peak_id][type][name]["start"} = $mir_start
          mir_annotations[peak_id][type][name]["end"} = $mir_end
          mir_annotations[peak_id][type][name]["overlap"} = "upstream"
        else:
          pass


# transcript annotation

# standard mirbase annotations

mir_transcript_annotations = collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(lambda: collection.defaultdict(str))))

for peak_id in peaks:
  chrom = peak_chrs[peak_id]
  start = peak_starts[peak_id]
  end = peak_ends[peak_id]
  
  start_bin = int((start - promoter_offset)) / bin_size)
  end_bin = int((end + promoter_offset) / bin_size)
  
  for bin in range(start_bin, end_bin + 1):
    if not bin in mir_types[chrom]:
      continue
    
    groups = mir_types[chrom][bin]
    
    for group in groups:
      names = mir_types[chrom][bin][group]

foreach my $peak_id (sort(keys(%peaks))) {
	my $chr = $peaks{$peak_id}{"chr"};
	my $start = $peaks{$peak_id}{"start"};
	my $end = $peaks{$peak_id}{"end"};

	foreach my $strand (sort(keys(%{$transcripts{$chr}}))) {
		foreach my $mir (sort(keys(%{$transcripts{$chr}{$strand}}))) {
			foreach my $transcript (sort(keys(%{$transcripts{$chr}{$strand}{$mir}}))) {
				
				# test whether there is overlap
				
				my %mir_annotation = %{$transcripts{$chr}{$strand}{$mir}{$transcript}};
				
				my $mir_start = $mir_annotation{"start"};
				my $mir_end = $mir_annotation{"end"};
				
				if ($strand eq "-") {
					# coordinates in forward space
					#$mir_start = $chromosomes{$chr} - $mir_start;
					#$mir_end = $chromosomes{$chr} - $mir_end;
					
					#my $t = $mir_start;
					#$mir_start = $mir_end;
					#$mir_end = $t;
				}

				if ($start >= $mir_start && $end <= $mir_end) {
					$mir_transcript_annotations{$peak_id}{$mir}{"chr"} = $chr;
					$mir_transcript_annotations{$peak_id}{$mir}{"start"} = $mir_start;
					$mir_transcript_annotations{$peak_id}{$mir}{"end"} = $mir_end;
					$mir_transcript_annotations{$peak_id}{$mir}{"overlap"} = "within";
				} elsif ($start < $mir_start && $end > $mir_end) {
					$mir_transcript_annotations{$peak_id}{$mir}{"chr"} = $chr;
					$mir_transcript_annotations{$peak_id}{$mir}{"start"} = $mir_start;
					$mir_transcript_annotations{$peak_id}{$mir}{"end"} = $mir_end;
					$mir_transcript_annotations{$peak_id}{$mir}{"overlap"} = "over";
				} elsif ($start < $mir_start && $end > $mir_start) {
                            $mir_transcript_annotations{$peak_id}{$mir}{"chr"} = $chr;
					$mir_transcript_annotations{$peak_id}{$mir}{"start"} = $mir_start;
					$mir_transcript_annotations{$peak_id}{$mir}{"end"} = $mir_end;
					$mir_transcript_annotations{$peak_id}{$mir}{"overlap"} = $strand eq "+" ? "upstream" : "downstream";
				} elsif ($start < $mir_end && $end > $mir_end) {
                            $mir_transcript_annotations{$peak_id}{$mir}{"chr"} = $chr;
					$mir_transcript_annotations{$peak_id}{$mir}{"start"} = $mir_start;
					$mir_transcript_annotations{$peak_id}{$mir}{"end"} = $mir_end;
					$mir_transcript_annotations{$peak_id}{$mir}{"overlap"} = $strand eq "+" ? "downstream" : "upstream";
				} else {
					# do nothing
				}
			}
		}
	}
}






# print them out

print "peak_id\tchr\tstart\tend\tlocation\tp-value\ttype\tmir\tmir_start\tmir_end\n";

foreach my $peak_id (sort(keys(%peaks))) {
  my $chr = $peaks{$peak_id}{"chr"};
	my $start = $peaks{$peak_id}{"start"};
	my $end = $peaks{$peak_id}{"end"};
	my $strand = $peaks{$peak_id}{"strand"};
	my $p = $peaks{$peak_id}{"p"};
  my $location = $chr . ":" . $start . "-" . $end;
  
  my $found = 0;
  
  #if (exists($mir_annotations{$peak_id}{"miRNA_primary_transcript"})) {
    #foreach my $mir (sort(keys(%{$mir_annotations{$peak_id}{"miRNA_primary_transcript"}}))) {
      ##my $chr = $mir_annotations{$peak_id}{"miRNA_primary_transcript"}{$name}{"chr"};
      ##my $start = $mir_annotations{$peak_id}{"miRNA_primary_transcript"}{$name}{"start"};
      ##my $end = $mir_annotations{$peak_id}{"miRNA_primary_transcript"}{$name}{"end"};
      ##my $overlap = $mir_annotations{$peak_id}{"miRNA_primary_transcript"}{$name}{"overlap"};
      
      #my $mir_start = $mir_annotations{$peak_id}{"miRNA_primary_transcript"}{$mir}{"start"};
      #my $mir_end = $mir_annotations{$peak_id}{"miRNA_primary_transcript"}{$mir}{"end"};
      #my $mir_overlap = $mir_annotations{$peak_id}{"miRNA_primary_transcript"}{$mir}{"overlap"};
      
      #print "p_" . $peak_id . "\t" . $p . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . $location . "\tgenomic_position\t" . $mir .  "\t" . $mir_start . "\t" . $mir_end . "\n";
    #}
  
		#$found = 1;
	#}
  
  if (exists($mir_annotations{$peak_id})) {
    foreach my $type (sort(keys(%{$mir_annotations{$peak_id}}))) {
      foreach my $mir (sort(keys(%{$mir_annotations{$peak_id}{$type}}))) {
        my $mir_start = $mir_annotations{$peak_id}{$type}{$mir}{"start"};
        my $mir_end = $mir_annotations{$peak_id}{$type}{$mir}{"end"};
        my $mir_overlap = $mir_annotations{$peak_id}{$type}{$mir}{"overlap"};
      
        print "p_", $peak_id, "\t", $chr, "\t", $start, "\t", $end, "\t", $location, "\t", $p, "\t", $type, "\t", $mir, "\t", $mir_overlap, "\t", $mir_start, "\t", $mir_end, "\n";
      }
		}
    
    $found = 1;
	}
  
	#if (exists($mir_annotations{$peak_id}{"miRNA"})) {
    #foreach my $mir (sort(keys(%{$mir_annotations{$peak_id}{"miRNA"}}))) {
      ##my $chr = $mir_annotations{$peak_id}{"miRNA"}{$name}{"chr"};
      #my $mir_start = $mir_annotations{$peak_id}{"miRNA"}{$mir}{"start"};
      #my $mir_end = $mir_annotations{$peak_id}{"miRNA"}{$mir}{"end"};
      #my $mir_overlap = $mir_annotations{$peak_id}{"miRNA"}{$mir}{"overlap"};
      
      #print "p_" . $peak_id . "\t" . $p . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . $location . "\tgenomic_position\t" . $mir .  "\t" . $mir_start . "\t" . $mir_end . "\n";
    #}
		
    #$found = 1;
	#}
  
  
  if (exists($mir_transcript_annotations{$peak_id})) {
    foreach my $mir (sort(keys(%{$mir_transcript_annotations{$peak_id}}))) {
      #my $chr = $mir_annotations{$peak_id}{"miRNA"}{$name}{"chr"};
      my $mir_start = $mir_transcript_annotations{$peak_id}{$mir}{"start"};
      my $mir_end = $mir_transcript_annotations{$peak_id}{$mir}{"end"};
      my $mir_overlap = $mir_transcript_annotations{$peak_id}{$mir}{"overlap"};
      
      print "p_", $peak_id, "\t", $chr, "\t", $start, "\t", $end, "\t", $location, "\t", $p, "\ttranscript\t",$mir, "\t", $mir_start, "\t", $mir_end, "\n";
    }
		
    $found = 1;
	}
  
  
  if ($found == 0) {
    print "p_", $peak_id, "\t", $chr, "\t", $start, "\t", $end, "\t", $location, "\t", $p, "\tn/a\tn/a\tn/a\tn/a\n";
  }
}
		
		
