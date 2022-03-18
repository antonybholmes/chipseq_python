# -*- coding: utf-8 -*-
"""
Create an xlsx version of a genes file

Created on Tue Sep 23 09:46:59 2014

@author: Antony Holmes
"""

import sys

import xlsxwriter

import excel

def create_xlsx(text_file, xlsx_file, locked):
  sys.stderr.write("Creating " + xlsx_file + " from " + text_file + "...\n")
  
  header = excel.get_header(text_file)
  
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
  
  center_indices = set()
  center_indices.add(0)
  center_indices.add(1)
  center_indices.add(2)
  
  workbook = xlsxwriter.Workbook(xlsx_file)

  worksheet = excel.create_worksheet(text_file, workbook, left_indices, right_indices, center_indices, locked)  

  
  #worksheet.set_column(0, columns - 1, 80)
  
  used = set()  
  
  used.add(excel.set_column_width(worksheet, 0, excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, 1, excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, 2, excel.SHORT_COLUMN_WIDTH))

  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Count", header), excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Genomic Locations (hg19)", header), excel.LONG_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region Relative To Gene", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Region TSS Closest Distance", header), excel.SHORT_COLUMN_WIDTH))
  
  excel.set_default_widths(worksheet, header, used)

  workbook.close()
  
  
file = sys.argv[1]
locked = len(sys.argv) > 2 and sys.argv[2] == "locked"

xlsx_file = excel.create_xlsx_file(file, locked)

create_xlsx(file, xlsx_file, locked)
