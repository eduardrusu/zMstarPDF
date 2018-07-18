# CE Rusu, July 18 2018
# Inspect the structure and content of a FITS table. This mimics the !head command for an ASCII file
# run as python inspectFITStable.py [table.fits] [number of head lines]

import sys
import fitsio # https://github.com/esheldon/fitsio

file = str(sys.argv[1])
head = int(str(sys.argv[2]))
table = fitsio.FITS(file)
columns = table[1].get_colnames()
read = table[1][0 : head - 1]

print(table)

strcolumns = ""
for i in range(len(columns)):
    strcolumns = strcolumns + columns[i] + "\t"
print strcolumns

for i in range(head - 1):
    strrows = ""
    for j in range(len(columns)):
        strrows = strrows + str(read[i][j]) + "\t"
    print strrows
