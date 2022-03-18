# -*- coding: utf-8 -*-
"""
Create an xlsx version of a genes file

Created on Tue Sep 23 09:46:59 2014

@author: Antony Holmes
"""

import sys
import re

import xlsxwriter

import excel

def create_xlsx(text_file, xlsx_file):
  sys.stderr.write("Creating " + xlsx_file + " from " + text_file + "...\n")
  workbook = xlsxwriter.Workbook(xlsx_file)
  worksheet = workbook.add_worksheet()

  left_font = workbook.add_format()
  left_font.set_font_size(11);
  left_font.set_font_name("Arial")
  left_font.set_align("left")
  left_font.set_align("vcenter")
  left_font.set_text_wrap()
  
  right_font = workbook.add_format();
  right_font.set_font_size(11);
  right_font.set_font_name("Arial")
  right_font.set_align("right")
  right_font.set_align("vcenter")
  right_font.set_text_wrap()

  center_font = workbook.add_format();
  center_font.set_font_size(11);
  center_font.set_font_name("Arial")
  center_font.set_align("center")
  center_font.set_align("vcenter") 
  center_font.set_text_wrap()
  
  header_font = workbook.add_format();
  header_font.set_font_size(11);
  header_font.set_font_name("Arial")
  header_font.set_bold()
  header_font.set_align("center")
  header_font.set_align("vcenter")
  header_font.set_text_wrap()
  
  r = 0

  f = open(text_file, 'r')
  
  header = f.readline().strip().split("\t")
  
  
  
  columns = len(header)
  
  left_indices = set()
  left_indices.add(excel.find_first_column("Region Relative To Gene", header))
  left_indices.add(excel.find_first_column("Region Relative To miR", header))
  left_indices.add(excel.find_first_column("Region Genomic Locations (hg19)", header))
  
  
			
  right_indices = set()
  right_indices.add(excel.find_first_column("Region TSS Distance", header))
  right_indices.add(excel.find_first_column("Region TSS Closest Distance", header))
  right_indices.add(excel.find_first_column("Region miR Start Closest Distance", header))
  right_indices.add(excel.find_first_column("Region miR Start Distance", header))
  right_indices.add(excel.find_first_column("Best P-value (ChIPseeqer)", header))

  #sys.stderr.write(";".join([str(d) for d in sorted(right_indices)]) + "\n")

  for i in range(0, columns):
    worksheet.write_string(r, i, header[i], header_font)

  r += 1
  
  for line in f:
    line = line.strip()
    
    if len(line) == 0:
      continue
    
    tokens = line.split("\t")
    
    for i in range(0, columns):
      if i in left_indices:
        font = left_font
      elif i in right_indices:
        font = right_font
      else:
        font = center_font
      
      worksheet.write_string(r, i, tokens[i], font)
    
    r += 1
    
  f.close()
  
  #worksheet.set_column(0, columns - 1, 80)
  
  used = set()  
  
  used.add(excel.set_column_width(worksheet, 0, excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, 1, excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, 2, excel.SHORT_COLUMN_WIDTH))

  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Count", header), excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Genomic Locations (hg19)", header), excel.LONG_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Relative To Gene", header), excel.MEDIUM_COLUMN_WIDTH))
  
  for i in range(0, columns):
    if i not in used:
      excel.set_column_width(worksheet, i, excel.DEFAULT_COLUMN_WIDTH)
  
  workbook.close()
  
  
file = sys.argv[1]

xlsx_file = re.sub(r'\.[^\.]+$', '.xlsx', file)

create_xlsx(file, xlsx_file)
