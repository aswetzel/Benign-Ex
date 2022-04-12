import sys, os
from collections import defaultdict

if len(sys.argv) != 3:
    print 'BREAKDGVINTOCHR.py <INFILE> <OUTPUT DIRECTORY>'
    sys.exit(2)


FILE = open(sys.argv[1].strip(), 'r').readlines()
DDCHR = defaultdict(list)
FOLDER = sys.argv[2].strip()

# Make OUTPUT Folder if it doesn't exist
if os.path.exists(FOLDER):
    pass
else:
    os.mkdir(FOLDER)


# For each line of INPUT file (skip header)
    # Skip lines where variant is an <inversion>
    # Skip lines where chrom is 'N', '<blank>', or contains '_'
HEADER = []
NOT = ["_","N","M"]
FIELDS = defaultdict(lambda : defaultdict(list))
for f in FILE:
    DATA = f.strip().replace('_',' ').split('\t')
    if 'variantaccession' in f.strip():
        HEADER = DATA[:7] + [DATA[8]] + DATA[14:17]

    if 'variantaccession' not in f.strip():
        if DATA[5].strip() != 'inversion': 
            DATA[1] = DATA[1].strip().replace(' ','_')
            if not DATA[1].strip() == "" and not any(x in DATA[1].strip() for x in NOT):
                if int(DATA[3])-int(DATA[2]) >= 2:
                    STATUS = DATA[5].strip()
                    TMP = DATA[:5] + [STATUS] + [DATA[6]] + [DATA[8]] + DATA[14:17]
                    for i, j in zip(HEADER, TMP): FIELDS[DATA[1].strip()][i.strip()].append(j.strip())

# Specify order of dictionary items
# Write data to individual CHR files in sorted order
ORDER = [ 'variantaccession', 'chr', 'start', 'end', 'varianttype', 'variantsubtype', 'reference', 'method', 'samplesize', 'observedgains', 'observedlosses']
for k in FIELDS.keys():
    OFILE = open(sys.argv[2].strip() + '/chr' + k.strip() + '_' + sys.argv[1].strip().split('/')[-1].replace('.txt', '_PROCESSED.txt'), 'w')

    TMP = [ FIELDS[k.strip()][y] for y in ORDER ]
    TMP = [ [int(x[2]), int(x[3]), '\t'.join(x).strip()] for x in zip(*TMP) ]
    TMP.sort(key=lambda x: x[1])
    TMP.sort(key=lambda x: x[0])
    for tmp in TMP:
        if len(tmp[2].split('\t')) == 11:
            OFILE.write(tmp[2].strip().replace(' ','_') + '\n')
    OFILE.close()
