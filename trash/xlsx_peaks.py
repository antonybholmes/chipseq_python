# -*- coding: utf-8 -*-
"""
Create an excel file of the peaks

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
  left_indices.add(4)
  left_indices.add(5)
  
  left_indices.add(excel.find_first_column("Entrez ID", header))
  left_indices.add(excel.find_first_column("Ensembl Gene ID", header))
  left_indices.add(excel.find_first_column("Ensembl Transcript ID", header))
  left_indices.add(excel.find_first_column("RefSeq", header))
  left_indices.add(excel.find_first_column("Gene Symbol", header))
  left_indices.add(excel.find_first_column("miR Symbol", header))
  
  left_indices.add(excel.find_first_column("Peak Relative To Gene", header))
  left_indices.add(excel.find_first_column("Peak Relative To miR", header))
  left_indices.add(excel.find_first_column("Peak Relative To Closest Gene", header))
			
  right_indices = set()
  right_indices.add(excel.find_first_column("Peak TSS Distance", header))
  right_indices.add(excel.find_first_column("Peak TSS Closest Distance", header))
  right_indices.add(excel.find_first_column("Peak miR Start Closest Distance", header))
  right_indices.add(excel.find_first_column("Peak miR Start Distance", header))
  right_indices.add(excel.find_first_column("Peak Size", header))
  right_indices.add(excel.find_first_column("P-value (ChIPseeqer)", header))
  right_indices.add(excel.find_first_column("Max Height (ChIPseeqer)", header))

  center_indices = set()  
  
  workbook = xlsxwriter.Workbook(xlsx_file)

  worksheet = excel.create_worksheet(text_file, workbook, left_indices, right_indices, center_indices, locked)  

  
  #worksheet.set_column(0, columns - 1, 80)
  
  used = set()  
  
  used.add(excel.set_column_width(worksheet, 0, excel.MEDIUM_COLUMN_WIDTH))
  
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Entrez ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Ensembl Gene ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Ensembl Transcript ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("RefSeq", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Gene Symbol", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("miR Symbol", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Closest Ensembl Gene ID", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Closest Ensembl Transcript ID", header), excel.MEDIUM_COLUMN_WIDTH))

  used.add(excel.set_column_width(worksheet, excel.find_first_column("Max Height (reads)", header), excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Peak Width", header), excel.SHORT_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Peak TSS Closest Distance", header), excel.SHORT_COLUMN_WIDTH))
 
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Peak Relative To Closest Gene", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Peak Relative To Gene", header), excel.MEDIUM_COLUMN_WIDTH))
  used.add(excel.set_column_width(worksheet, excel.find_first_column("Peak Relative To miR", header), excel.MEDIUM_COLUMN_WIDTH))
  
  excel.set_default_widths(worksheet, header, used)

  workbook.close()
  
  
file = sys.argv[1]
locked = len(sys.argv) > 2 and sys.argv[2] == "locked"

xlsx_file = excel.create_xlsx_file(file, locked)

create_xlsx(file, xlsx_file, locked)
