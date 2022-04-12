import sys, os
from collections import defaultdict

if len(sys.argv) != 5:
    print 'BREAKINTOCHR.py <INFILE> <OUTPUT DIRECTORY> <HEADER; Y or N> <DATABASE; DGV or CLINGEN>'
    sys.exit(2)

DDCHR = defaultdict(lambda : defaultdict(list))
FIELDS = defaultdict(lambda : defaultdict(list))
ALLFIELDS = defaultdict(lambda : defaultdict(list))
FOLDER = sys.argv[2].strip()
FILE = open(sys.argv[1].strip(), 'r').readlines()


# Make OUTPUT Folder if it doesn't exist
if os.path.exists(FOLDER):
    pass
else:
    os.mkdir(FOLDER)

#######################
# FORMATTING DGV DATA #
#######################
if sys.argv[4] == "DGV":

 # Get Header or Supply Default Header
 if sys.argv[3] == 'Y':
    for f in FILE:
        HEADER = f.strip().split('\t')
        if len(HEADER) == 20:
            print '\nERROR: MULTIPLE INPUT FILE FORMAT ERRORS FOUND. INPUT FILE SHOULD CONTAIN 11 COLUMNS.'
            print 'IF USING DOWNLOADED DGV DATABASE PLEASE USE COMMAND:', '\n', '  ', 'python FORMATDGV.py', sys.argv[1].strip(), sys.argv[2].strip()
            sys.exit()

        elif len(HEADER) != 11: 
            print 'ERROR: MULTIPLE INPUT FILE FORMAT ERRORS FOUND. INPUT FILE SHOULD CONTAIN 11 COLUMNS. PLEASE REVIEW README FILE.'
            sys.exit()
        FILE = FILE[1:]
        break

 elif sys.argv[3] == 'N': 
    HEADER = [ 'variantid', 'chr', 'start', 'end', 'varianttype', 'variantsubtype', 'reference', 'method', 'samplesize', 'observedgains', 'observedlosses']
    for each in FILE[0].strip().split('\t'):
        if each in HEADER:
            print '\nERROR: POTENTIAL HEADER DETECTED BUT USER INDICATES NO HEADER IN FILE.', '\n', FILE[0].strip()
            sys.exit()
 ERR=[]
 ctr=0
 NOT = ["_","N","M"]

 # For Each Line of File
 for f in FILE:
    ctr=ctr+1
    DATA = f.strip().split('\t')
    try: DATA[5].strip()
    except: 
        print "\nERROR: PLEASE REMOVE EMPTY LINE AT END OF FILE."
        sys.exit()

    # Skip lines where variant is an <inversion>
    # Skip lines where chrom is 'N', <blank>, or contains '_'
    if len(DATA) == 11 and DATA[5].strip() != 'inversion':
        if not DATA[1].strip() == "" and not any(x in DATA[1].strip() for x in NOT):
            STATUS = DATA[5].strip()
            TMP = DATA[:5] + [STATUS] + DATA[6:]
            for i, j in zip(HEADER, TMP): ALLFIELDS[DATA[1].strip()][i.strip()].append(j.strip())
            for k, l in zip(HEADER, TMP): FIELDS[DATA[1].strip()][k.strip()].append(l.strip())

        # Add lines that don't meet criteria to counter      
        else:
            STATUS = DATA[5].strip()
            TMP = DATA[:5] + [STATUS] + DATA[6:]
            for i, j in zip(HEADER, TMP): ALLFIELDS[DATA[1].strip()][i.strip()].append(j.strip())

    # Add to 'ERR' array if missing data
    # Add to ALLFIELDS dictionary to count where missing data is
    elif len(DATA) != 11:
        ERR.append(ctr)
        STATUS = DATA[5].strip()
        TMP = DATA[:5] + [STATUS] + DATA[6:]
        for i, j in zip(HEADER, TMP): ALLFIELDS[DATA[1].strip()][i.strip()].append(j.strip())

    # Adds lines that don't meet crtieria to counter for 
    else:
        STATUS = DATA[5].strip()
        TMP = DATA[:5] + [STATUS] + DATA[6:]
        for i, j in zip(HEADER, TMP): ALLFIELDS[DATA[1].strip()][i.strip()].append(j.strip())

# Error Handler for Missing Data
 SUM = defaultdict(list)
 if len(ERR) > 1:
    print 'ERROR: MULTIPLE INPUT FILE FORMAT ERRORS FOUND. PLEASE REVIEW README FILE. ' + repr(len(ERR)) + ' LINES WITH INCORRECT FORMATTING WERE SKIPPED'
    print 'Number of Items Found for Each Column:'

    for k in ALLFIELDS.keys():
        for y in HEADER: 
            SUM[y].append(len(ALLFIELDS[k][y]))
    for z in HEADER:
        print '\t', sum(SUM[z]), '\t', z 
    print 'TOTAL LINES IN FILE:', len(FILE)
    sys.exit()

 elif len(DATA) == 1:
    print 'ERROR: INPUT FILE FORMAT INCORRECT AT LINE ' + repr(ctr) + '. PLEASE REVIEW README FILE.'
    print 'LINE WITH INCORRECT FORMAT WAS SKIPPED'
    print f.strip()
    sys.exit()



 # Specify order of dictionary items
 # Write data to individual CHR files in sorted order
 for k in FIELDS.keys():
  try:
    OFILE2 = open(sys.argv[2].strip() + '/chr' + k.strip() + '_' + sys.argv[1].strip().split('/')[-1].replace('.txt', '_PROCESSED.txt'), 'w')
    TMP = [ FIELDS[k.strip()][y] for y in HEADER ]
    TMP = [ [int(x[2]), int(x[3]), '\t'.join(x).strip().replace('_',' ')] for x in zip(*TMP) ]
    TMP.sort(key=lambda x: x[1])
    TMP.sort(key=lambda x: x[0])
    for tmp in TMP: OFILE2.write(tmp[2].strip() + '\n')
    OFILE2.close()
  except ValueError:
    if sys.argv[3] == "N": 
        print '\nERROR: POTENTIAL HEADER DETECTED BUT USER INDICATES NO HEADER IN FILE.', '\n', TMP
        sys.exit()

###########################
# FORMATTING CLINGEN DATA #
###########################
elif sys.argv[4] == "CLINGEN":
 NOT = ["_","N"]

 BENIGN = ["B","Benign"]
 USLB = ["USLB","LB","Likely benign","LikelyBenign","Uncertain significance: likely benign","Uncertain significance:likely benign","Uncertain significance: likely benign"]
 VUS = ["VUS","Uncertain","Uncertain Significance","Variant of unknown significance","Variant of unclear clinical significance"]

 GAIN = ["Gain","Duplication"]
 LOSS = ["Loss","Deletion"]

 for f in FILE:
    HEADER = f.strip().split('\t')
    if len(HEADER) != 6: 
        print 'ERROR: MULTIPLE INPUT FILE FORMAT ERRORS FOUND. INPUT FILE SHOULD CONTAIN 6 COLUMNS. PLEASE REVIEW README FILE.'
        sys.exit()
    break

 # For each line of file
 # Standardize Classification { Error Handling if non-acceptable classification provided }
 if sys.argv[3] == 'Y': FILE = FILE[1:]
 ctr=0
 for f in FILE:
    ctr=ctr+1
    DATA = f.strip().split('\t')
    if not DATA[0].strip() == "" and not any(x in DATA[0].strip() for x in NOT):

      if DATA[5].strip().lower() in (gain.lower() for gain in GAIN):
        if DATA[4].strip().lower().replace('_',' ') in (name.lower() for name in BENIGN): 
            NEWDATA = DATA[0:4] + ['Benign'] + [DATA[5].lower()]
            DDCHR[DATA[0].strip()][DATA[-1].strip().lower()].append('\t'.join(NEWDATA))
        elif DATA[4].strip().lower().replace('_',' ') in (name.lower() for name in USLB):
            NEWDATA = DATA[0:4] + ['Uncertain significance: likely benign'] + [DATA[5].lower()]
            DDCHR[DATA[0].strip()][DATA[-1].strip().lower()].append('\t'.join(NEWDATA))
        elif DATA[4].strip().lower().replace('_',' ') in (name.lower() for name in VUS):
            NEWDATA = DATA[0:4] + ['Uncertain_Significance'] + [DATA[5].lower()]
            DDCHR[DATA[0].strip()][DATA[-1].strip().lower()].append('\t'.join(NEWDATA))
        else: 
            print 'ERROR: INCORRECT CLASSIFICATION...', DATA[4], 'AT LINE ' + repr(ctr) + '. PLEASE REVIEW README FILE.'
            sys.exit()

      elif DATA[5].strip().lower() in (loss.lower() for loss in LOSS):
        if DATA[4].strip().lower().replace('_',' ') in (name.lower() for name in BENIGN): 
            NEWDATA = DATA[0:4] + ['Benign'] + [DATA[5].lower()]
            DDCHR[DATA[0].strip()][DATA[-1].strip().lower()].append('\t'.join(NEWDATA))
        elif DATA[4].strip().lower().replace('_',' ') in (name.lower() for name in USLB):
            NEWDATA = DATA[0:4] + ['Uncertain significance: likely benign'] + [DATA[5].lower()]
            DDCHR[DATA[0].strip()][DATA[-1].strip().lower()].append('\t'.join(NEWDATA))
        elif DATA[4].strip().lower().replace('_',' ') in (name.lower() for name in VUS):
            NEWDATA = DATA[0:4] + ['Uncertain_Significance'] + [DATA[5].lower()]
            DDCHR[DATA[0].strip()][DATA[-1].strip().lower()].append('\t'.join(NEWDATA))
        else: 
            print 'ERROR: INCORRECT CLASSIFICATION...', DATA[4], 'AT LINE ' + repr(ctr) + '. PLEASE REVIEW README FILE.'
            sys.exit()

      elif DATA[5].strip().lower() == "direction": 
        print '\nERROR: POTENTIAL HEADER DETECTED BUT USER INDICATES NO HEADER IN FILE.', '\n', f.strip()
        sys.exit()

      else: 
        print 'ERROR: INCORRECT DIRECTION...', DATA[5], 'AT LINE ' + repr(ctr) + '. PLEASE REVIEW README FILE.'
        sys.exit()

 for k in DDCHR.keys():
    GFILE = open(sys.argv[2].strip() +'/' + k.strip() + '_' + sys.argv[1].strip().split('/')[-1].replace('.txt', '_PROCESSED_GAIN.txt'), 'w')      
    TMP = [ [int(x.strip().split('\t')[1]),int(x.strip().split('\t')[2]), x.strip().split('\t')[3], x.strip()] for x in DDCHR[k.strip()]['gain'] ]
    TMP.sort(key=lambda x: x[2])
    TMP.sort(key=lambda x: x[1])
    TMP.sort(key=lambda x: x[0])
    for tmp in TMP: GFILE.write(tmp[3].strip() + '\n')
    GFILE.close()

    LFILE = open(sys.argv[2].strip() + '/' + k.strip() + '_' + sys.argv[1].strip().split('/')[-1].replace('.txt', '_PROCESSED_LOSS.txt'), 'w')
    TMP = [ [int(x.strip().split('\t')[1]),int(x.strip().split('\t')[2]), x.strip().split('\t')[3], x.strip()] for x in DDCHR[k.strip()]['loss'] ]
    TMP.sort(key=lambda x: x[2])
    TMP.sort(key=lambda x: x[1])
    TMP.sort(key=lambda x: x[0])
    for tmp in TMP: LFILE.write(tmp[3].strip() + '\n')
    LFILE.close()

