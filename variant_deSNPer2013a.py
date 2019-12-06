#Written in python 2 in 2013 by Scott Kennedy.

import sys
from collections import defaultdict
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", action="store", type="string", dest="SNPfile", default=None)
(o, args) = parser.parse_args()

SNPfile = open(o.SNPfile, 'r')
SNPList = []

for SNPline in SNPfile :
	try:
		SNPList.append((SNPline.split('\t')[0], SNPline.split('\t')[1]))
	except IndexError:
		print "indexerror", SNPline
		continue
for varLine in sys.stdin :
    if (varLine.strip('\n').split('\t')[0],varLine.strip('\n').split('\t')[1]) not in SNPList and 'N' not in (varLine.strip('\n').split('\t')[2]) and 'N' not in (varLine.strip('\n').split('\t')[3]):
		print(varLine.strip('\n'))

SNPfile.close()
