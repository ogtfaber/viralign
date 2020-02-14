# pipreqs . --force

import sys
import csv

# sequence name, sequence length, # mapped read-segments and # unmapped read-segments
idxstats_filename = 'reads-aligned.virus.sorted.idxstats'
sqlite_file = 'U-RVDBv17.0.sqlite'

data_directory = '/data'
ref_directory = '/ref'
# data_directory = '/Users/onnofaber/Projects/rm/viralign/sample'
# ref_directory = '/Users/onnofaber/Projects/rm/viralign/ref'

reader = csv.reader(open(f'{data_directory}/{idxstats_filename}.csv'), delimiter="\t")

header = [
  "sequence_name",
  "sequence_length",
  "mapped_readsegs",
  "unmapped_readsegs",
]

# import operator
# sortedlist = sorted(reader, key=operator.itemgetter(2), reverse=True)
sortedlist = sorted(reader, key=lambda row: int(row[header.index('mapped_readsegs')]), reverse=True)

# normalize and calculate totals

total_mapped_readsegs = 0
total_unmapped_readsegs = 0
for index in range(len(sortedlist)):
  row = sortedlist[index]
  sequence_length = int(row[header.index('sequence_length')])
  mapped_readsegs = int(row[header.index('mapped_readsegs')])
  unmapped_readsegs = int(row[header.index('unmapped_readsegs')])
  row[header.index('sequence_length')] = sequence_length
  row[header.index('mapped_readsegs')] = mapped_readsegs
  row[header.index('unmapped_readsegs')] = unmapped_readsegs
  total_mapped_readsegs += mapped_readsegs
  total_unmapped_readsegs += unmapped_readsegs
  sortedlist[index] = row

total_readsegs = total_mapped_readsegs + total_unmapped_readsegs

# mapped_reads/(reference genome size * totally mapped reads)
print('total_mapped_readsegs', total_mapped_readsegs)
print('total_unmapped_readsegs', total_unmapped_readsegs)
print('total_readsegs', total_readsegs)
print('% of mapped readsegs', (total_mapped_readsegs / total_readsegs) * 100)

def has_reads(row):
  return row[header.index('sequence_length')] > 0 and row[header.index('mapped_readsegs')] > 0

sortedlist_filtered = list(filter(has_reads, sortedlist))

header2 = header + [
  "rpkm_norm",
]
for index in range(len(sortedlist_filtered)):
  row = sortedlist_filtered[index]
  mapped_readsegs = row[header.index('mapped_readsegs')]
  sequence_length = row[header.index('sequence_length')]
  size_kb = sequence_length / 10**3
  rpkm = mapped_readsegs / (size_kb * (total_mapped_readsegs / 10**6))
  row2 = row + [rpkm]
  sortedlist_filtered[index] = row2

import sqlite3
conn = sqlite3.connect(f'{ref_directory}/{sqlite_file}.db')

def get_columns_annotations():
  cursor=conn.cursor()
  cursor.execute('SELECT * FROM rvdb LIMIT 1')
  column_names = list(map(lambda x: x[0], cursor.description))
  return column_names

def get_locus_annotations(locus_name):
  cursor=conn.cursor()
  cursor.execute('SELECT * FROM rvdb WHERE accs=?', [locus_name])
  result = cursor.fetchone()
  return result

header3 = header2 + get_columns_annotations()

for index in range(len(sortedlist_filtered)):
  row = sortedlist_filtered[index]
  locus_name = row[0].split('|')[2]
  annotations = get_locus_annotations(locus_name)
  row2 = row + list(annotations)
  sortedlist_filtered[index] = row2

with open(f'{data_directory}/{idxstats_filename}.processed.csv', 'w', newline='') as myfile:
  wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
  wr.writerow(header3)
  for index in range(len(sortedlist_filtered)):
    row = sortedlist_filtered[index]
    wr.writerow(row)
