import numpy
import csv
import sys
import os


gal = list(csv.reader(open(sys.argv[1], "rb"), delimiter=','))
ndata = len(gal)
header = numpy.array(gal[0])
gal[1:] = sorted(gal[1:],key=lambda l:l[numpy.where(header == "image")[0][0]])
path = "Field/"

for line in range(1,ndata):
	imgIndex = numpy.where(header == "image")[0][0]
        fieldName = path+gal[line][imgIndex]
        fileName = gal[line][imgIndex].replace(".gz", "")
	if not(os.path.isfile(fieldName) or os.path.isfile(path+fileName)):
                print("Downloading", line)
		ra = gal[line][numpy.where(header == "ra")[0][0]]
		dec = gal[line][numpy.where(header == "dec")[0][0]]
		run = gal[line][numpy.where(header == "run")[0][0]]
		rerun = gal[line][numpy.where(header == "rerun")[0][0]]
		camcol = gal[line][numpy.where(header == "camcol")[0][0]]
		#field = gal[line][numpy.where(header == "field")[0][0]]
		dr7id = gal[line][numpy.where(header == "dr7objid")[0][0]]
            
		cmd = "wget --inet4-only -r -nd --directory-prefix=Field http://das.sdss.org/raw/"
		cmd += str(run) + "/"
		cmd += str(rerun) + "/corr/"
		cmd += str(camcol) + "/"
		cmd += fileName + ".gz"
                print(cmd)
		pr = os.popen(cmd)
		print(pr.read())
		# unzip the image
		cmd = "gzip -d " + path + fileName + ".gz"
		pr = os.popen(cmd)
		pr.read()
	else:
		print("Found",fieldName,line)
print("Done")
