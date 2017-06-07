import csv
import sys
import numpy

gal1 = list(csv.reader(open(sys.argv[1], "rb"), delimiter=','))
header1 = numpy.array(gal1[0])
gal2 = list(csv.reader(open(sys.argv[2], "rb"), delimiter=','))
header2 = numpy.array(gal2[0])

for i in range(1,len(gal))
