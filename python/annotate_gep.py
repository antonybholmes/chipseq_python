import sys

import lib.expression


file = sys.argv[1]
col = int(sys.argv[2])

expression_list = []
expression_list_headers = []

expression_list_headers.append("GEP Affy GCvsN")
expression_list_headers.append("GEP Affy GCvsM")
expression_list.append(lib.expression.AffyGeneCBvsNExpression())
expression_list.append(lib.expression.AffyGeneCBvsMExpression())

expression_list_headers.append("GEP Affy DZvsLZ")
expression_list.append(lib.expression.AffyGeneDZvsLZExpression())

# RNA-seq
expression_list_headers.append("RNA-seq GCvsN")
expression_list_headers.append("RNA-seq GCvsM")
expression_list.append(lib.expression.RnaSeqGeneCBvsNExpression())
expression_list.append(lib.expression.RnaSeqGeneCBvsMExpression())

# RNA-seq from deseq2
expression_list_headers.append("RNA-seq GCvsN DeSeq2")
expression_list_headers.append("RNA-seq GCvsM DeSeq2")
expression_list.append(lib.expression.RnaSeqGeneCBvsNDeSeq2Expression())
expression_list.append(lib.expression.RnaSeqGeneCBvsMDeSeq2Expression())

f = open(file, 'r')

header = f.readline().strip().split("\t")

for h in expression_list_headers:
  header.append(h)
  
sys.stdout.write("\t".join(header) + "\n")

for line in f:
  tokens = line.strip().split("\t")
  
  gene = tokens[col]
  
  for expression in expression_list:
    tokens.append(expression.get_expression(gene))
    
  sys.stdout.write("\t".join(tokens) + "\n")
  
f.close()
