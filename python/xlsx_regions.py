# -*- coding: utf-8 -*-
"""
Create an excel file of the regions

Created on Tue Sep 23 13:34:02 2014

@author: Antony Holmes
"""

import sys

import xlsxwriter

import excel

def create_xlsx(text_file, xlsx_file, locked):
  sys.stderr.write("Creating " + xlsx_file + " from " + text_file + "...\n")
  
  header = excel.get_header(text_file)
  
  left_indices = set()
  left_indices.add(0)
  
  # Ensure the sample columns are left aligned
  c = excel.find_first_column("Sample", header)
  
  for i in range(0, 3):
    if "Sample" not in header[c]:
      break
    
    left_indices.add(c)
    
    c += 1

  left_indices.add(excel.find_first_column("Entrez ID", header))
  left_indices.add(excel.find_first_column("Ensembl Gene ID", header))
  left_indices.add(excel.find_first_column("Entrez Transcript ID", header))
  left_indices.add(excel.find_first_column("RefSeq", header))
  left_indices.add(excel.find_first_column("Gene Symbol", header))
  left_indices.add(excel.find_first_column("miR Symbol", header))
  
  left_indices.add(excel.find_first_column("Region Relative To Gene", header))
  left_indices.add(excel.find_first_column("Region Relative To miR", header))
  left_indices.add(excel.find_first_column("Region Relative To Closest Gene", header))
			
  right_indices = set()
  right_indices.add(excel.find_first_column("Region TSS Distance", header))
  right_indices.add(excel.find_first_column("Region TSS Closest Distance", header))
  right_indices.add(excel.find_first_column("Region miR Start Closest Distance", header))
  right_indices.add(excel.find_first_column("Region miR Start Distance", header))
  right_indices.add(excel.find_first_column("Peak / Overlap Width", header))
  right_indices.add(excel.find_first_column("Best P-value (ChIPseeqer)", header))
  
  center_indices = set()
  
  workbook = xlsxwriter.Workbook(xlsx_file)

  worksheet = excel.create_worksheet(text_file, workbook, left_indices, right_indices, center_indices, locked)
  
  used = set()  
  
  used.add(excel.set_column_width(worksheet, 0, excel.MEDIUM_COLUMN_WIDTH))
  
  c = excel.find_first_column("Sample", header)
  
  for i in range(0, 3):
    if "Sample" not in header[c]:
      break
    
    used.add(excel.set_column_width(worksheet, c, excel.MEDIUM_COLUMN_WIDTH))
    
    c += 1
  
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Entrez ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Ensembl Gene ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Ensembl Transcript ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("RefSeq", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Gene Symbol", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("miR Symbol", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Closest Ensembl Gene ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Closest Ensembl Transcript ID", header), excel.MEDIUM_COLUMN_WIDTH))


  used.add(excel.set_column_width(worksheet, excel.find_first_column("Peak / Overlap Width", header), excel.SHORT_COLUMN_WIDTH))
  #used.add(excel.set_column_width(worksheet, excel.find_first_column("Number Of Overlapping Peaks", header), excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region TSS Closest Distance", header), excel.SHORT_COLUMN_WIDTH))
  
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Count", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Relative To Closest Gene", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Relative To Gene", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Relative To miR", header), excel.MEDIUM_COLUMN_WIDTH))
  
  
  
  excel.set_default_widths(worksheet, header, used)

  workbook.close()
  
  
file = sys.argv[1]
locked = len(sys.argv) > 2 and sys.argv[2] == "locked"

xlsx_file = excel.create_xlsx_file(file, locked)

create_xlsx(file, xlsx_file, locked)
