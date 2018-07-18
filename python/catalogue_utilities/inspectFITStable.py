# CE Rusu, July 18 2018
# Inspect the structure and content of a FITS table. This mimics the UNIX !head/!tail command for an ASCII file
# run as python inspectFITStable.py [table.fits] [head/tail] [number of head lines]

import sys
import fitsio # https://github.com/esheldon/fitsio

file = str(sys.argv[1])
mode = str(sys.argv[2])
number = int(str(sys.argv[3]))
table = fitsio.FITS(file)
columns = table[1].get_colnames()

if mode == 'head':
    read = table[1][0 : number]
if mode == 'tail':
    read = table[1][int(table[1].get_nrows()) - number : int(table[1].get_nrows())]

print(table)
print "Number of rows: ", table[1].get_nrows()

strcolumns = ""
for i in range(len(columns)):
    strcolumns = strcolumns + columns[i] + "\t"
print strcolumns

for i in range(number):
    strrows = ""
    for j in range(len(columns)):
        strrows = strrows + "%.5e \t" % read[i][j])
    print strrows
