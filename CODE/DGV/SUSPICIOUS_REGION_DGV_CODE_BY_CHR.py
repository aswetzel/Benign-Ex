#load required python modules
from __future__ import division 
import sys, glob, os, networkx
from collections import defaultdict
from progress.bar import Bar
from intervaltree import Interval, IntervalTree

# Error Handler if Incorrect Number of Arguments
if len(sys.argv) != 3:
    print "SUSPICIOUS_REGION_DGV_CODE.py <RAWFILES PATH> <PARAMETERFILES PATH>"
    sys.exit(1)

GAINKEY = ['gain', 'duplication', 'gain+loss']
LOSSKEY = ['loss', 'deletion', 'gain+loss']

#print "Check for exact phrashe"
#print "compre heatmaps to see if significant difference in which is optimal"
#sys.exit()

try:

 ALLFILES = glob.glob(sys.argv[2].strip() + '/CHR*/*.txt')
 ALLFILES = [ ls for ls in ALLFILES if not 'ADDITIONALREGIONS' in ls if not 'FINALREGIONS' in ls ] # Remove ADDITIONALREGIONS files from list set if they exist
 SUM=len(ALLFILES)*3

 # Function to add dgv data to DICT_DGV
    # For each DGV chr file (DGV_SOURCE/HGXX_XXXX)
        # Open chr file; put each tab-delimted line in array 'DATA'
        # Loop over YEARS
            # If year filter APPLIED: If study year >= {YEAR FILTER} then create dictionary key 'chr#' and associate dgv data as 'value' 
            # If year filter NOT applied: create dictionary key 'chr#' and associate dgv data as 'value' { value = ID, chr, start, end, dir, n, #gain, #loss } 
 def RAWFILE(d):
    FILE = open(d.strip(), 'r')
    for line in FILE:
        DATA = filter(None, line.strip().split('\t'))
        iYEAR = DATA[6].strip().replace(' ','_')
        if 'gnomAD' in iYEAR: yr='2019'
        elif '1000_Genomes_Consortium_Pilot_Project' in iYEAR: yr='2010'
        elif '1000_Genomes_Consortium_Phase_1' in iYEAR: yr='2012'
        elif '1000_Genomes_Consortium_Phase_3' in iYEAR: yr='2015'
        else: yr = iYEAR.split('_')[-1].strip().replace('b','')

        for x in YEARS:
            if ( (x.isdigit()) and (int(yr) >= int(x)) ): 
                CHR = DATA[1].strip()
                DICT_DGV[x].append(DATA[0:4] + [DATA[5]] + DATA[8:])
            elif not x.isdigit():
                CHR = DATA[1].strip()
                DICT_DGV[x].append(DATA[0:4] + [DATA[5]] + DATA[8:])

        # Add dictionary entry which appends the variant, type, lesion, study, method, sample size, gain freq, loss freq
        ele = filter(None, line.strip().split('\t')) # split each line put into elements of list
        v = ele[2] + '_' + ele[3]                 # start_stop
        middle = 'Type:' + ele[4] + ' ' +\
             'Lesion:' + ele[5] + ' ' +\
             'Study:' + ele[6] + ' ' +\
             'Method:' + ele[7].replace(' ','_') + ' ' +\
             'Sample_size:' + ele[8]          # Middle = Type:<> & Lesion:<> & Study:<> & Method:<> & SampleSize: <>
        dict_all[v].append('Name:' + ele[0] + ' ' + middle + ' ' + 'Freq_Gain:' + repr(int(ele[9]) / int(ele[8])) + ' ' + 'Freq_Loss:' + repr(int(ele[10]) / int(ele[8])))
    FILE.close()

 # Create function to lookup how many dgv entries overlap the 500bp interval
 def compare_benign(lk, lo, listn, counter, DDCB):
    region = lk[0] + '-' + lk[1]
    for a in lo:
        if len(range(max(int(lk[0]), int(a[0])), min(int(lk[1]), int(a[1])) + 1)) > 1:
            DDCB[region].append([int(a[0]), int(a[1])])
            counter += 1
    listn.append([region, counter])
    return listn

 # Function to Determine which GAIN/LOSS Intervals to Add as Reciprocal Events
 def PROCESSINPUT(pf):
        global COUNT_TOTAL
        COUNT_TOTAL += 1 # Keep track of the total number of non 'additional regions' files processed
        READPFILE = open(pf.strip(), 'r').readlines() # read file
        GAINLIST = []
        CHRN = 0
        LOSSLIST = []
        DDO = defaultdict(set)
        OUTPUT = []
        ALLDATA = defaultdict(str)
        FILENAME = pf.split('/')[-1].strip() # Get Filename
        CHRN = FILENAME.split('CHR')[1].split('_')[0].strip() # Get Chromosome
        EVENTSGL = int(FILENAME.split('_')[1].strip())  # Get Min Event Freq
        FREQ = float(FILENAME.split('_')[4]) # Get Min Freq
        if pf.strip().split('_')[-1].replace('.txt','').isdigit(): PF_YEAR = pf.strip().split('_')[-1].replace('.txt','') # Get YEAR for year filter for file (if applied)
        else: PF_YEAR = 'NA'
        RAWD = DICT_DGV[PF_YEAR] # Get Dictionary Entries for that chromosome
        DICT_DGV_CHR = defaultdict(set)

        # For each dictionary entry : Associate dgv variant name with coordiantes
        for rawd in RAWD:
            DICT_DGV_CHR[rawd[0].strip()].add(':'.join(rawd[2:4]))

        # If the FILENAME has 'GAIN' THEN for each LINE <pfi in READPFILE>:
        if 'GAIN' in FILENAME:
            for pfi in READPFILE:
                DATAPFI = pfi.strip().split('\t') # Split data into array
                ATLEASTBL = []
                CALLED_COORD = []

                # Add line (Gain Benign Region) to <GAIN> DD_BEXREGIONS with FILENAME as key; LINE as VALUE
                DD_BEXREGIONS[pf.strip()].append(pfi.strip())

                # If UNCALLED column exists THEN put CALLED data into CALLED array; Variant name & Type into UNCALLED array
                    # If type is 'Loss' or 'Gain+Loss' get COORDINATES from DICT_DGV_CHR; add COORDINATES to ATLEASTBL array
                    # If variant length >= length of benign region then add to FILTERATLEASTBL array
                if DATAPFI[-1].strip().startswith('UNCALLED: '): 
                    CALLED = [ call.replace('CALLED: ','').split(' ') for call in DATAPFI[-2].split(';') ]
                    UNCALLED = [ [ uncall.replace('UNCALLED: ','').split(' ')[0].split(':')[1].strip(), uncall.replace('UNCALLED: ','').split(' ')[2].split(':')[1].strip() ] for uncall in DATAPFI[-1].split(';') ]
                    for atbl in UNCALLED: 
                        ATLEASTBL.extend( i for i in list(DICT_DGV_CHR[atbl[0].strip()]) if atbl[1].lower() in LOSSKEY )

                # Else put data into CALLED; UNCALLED array is blank
                else:
                    CALLED = [ call.replace('CALLED: ', '').split(' ') for call in DATAPFI[-1].split(';') ]
                    UNCALLED = []

                # If CALLED Gain Freq >= 1% then add to RESULTG array
                RESULTG = [ '||'.join(rg).replace('CALLED:||', '') for rg in CALLED if float(rg[-2].split(':')[1]) * 100 >= 1 ]

                if len(RESULTG) > 0:
                   # For called dgv entry; get coordinates; Associate Coordiantes with ID
                   GCALLED = []
                   CALLEDIDG = [ rg.split('||')[0].replace('Name:','') for rg in RESULTG ]

                   for gr in RESULTG:
                      for ddcn in list(DICT_DGV_CHR[filter(None, gr.split('||'))[0].replace('Name:', '').strip()]):
                         COORD = [int(ddcn.split(':')[0]), int(ddcn.split(':')[1])]
                         CALLED_COORD.append(COORD)
                         GCALLED.append([gr.split('||')[0].replace('Name:',''), COORD])

                   if len(CALLED_COORD) > 1:
                      # For each pair of CALLED_COORD(n,n+1); Get Max Start & Min End; If CNVs Overlap (Range >1); If CNVS bookended (Range=1); If CNVs DO NOT Overlap (Range=0)
                      # If CNVs Overlap Save to CHECKOVERLAP array
                      # Get Unique Coordinates in RESULT_OVERL
                      CHECKOVERLAP = [ ['_'.join([repr(cco[0]), repr(cco[1])]), '_'.join([repr(rd[0]), repr(rd[1])])] for cco, rd in zip(CALLED_COORD, CALLED_COORD[1:]) if len(range(max(cco[0], rd[0]), min(cco[1], rd[1]) + 1)) > 1 ]
                      RESULT_OVERL = list(set(sum(CHECKOVERLAP, [])))

                      # If Overlapping CNVs with >1% AF were Found:
                      if len(CHECKOVERLAP) > 2:
                   
                         # Sort all START/END Coordinates; Get Regions Between Each Start/End Coordinate
                         S = [ int(RO.split('_')[0]) for RO in RESULT_OVERL ]
                         E = [ int(RO.split('_')[1]) for RO in RESULT_OVERL ]
                         ALL = list(set(S+E))
                         ALL.sort(key=int)
                         all_windows = []
                         i=0; j=1
                         while i < len(ALL)-1:
                            while j <= len(ALL)-1:
                               all_windows.append([int(ALL[i]), int(ALL[j])])
                               i = i + 1; j = j + 1

                         # Find out how many >1% AF CNVs overlap each region
                         combine = [ [RO.split('_')[0], RO.split('_')[1]] for RO in RESULT_OVERL ]
                         targets_count = []; countert = 0; DDCB = defaultdict(list)
                         for each in all_windows:
                            compare_benign([repr(each[0]), repr(each[1])], combine, targets_count, countert, DDCB)

                         # If Number of OVERLAPPING CALLED GAIN variants with >1% Frequency is 2x Min Number of EVENTS
                         for y in targets_count: 
                            if y[1] >= 2* EVENTSGL:

                               # Check if At Least ONE UNCALLED LOSS REGION OVERLAPS
                               for fatbl in ATLEASTBL:
                                  if len(range(int(max(y[0].split('-')[0], fatbl.split(':')[0])), int(min(y[0].split('-')[1], fatbl.split(':')[1])) + 1)) > 1: 

                                     # Get DGV info for supporting GAIN variants
                                     CALLEDINDEX = []
                                     for ddcb in DDCB[y[0]]: CALLEDINDEX.extend([ item for item in dict_all['_'.join(map(str,ddcb))] if item.split(' ')[0].replace('Name:','') in CALLEDIDG ])

                                     # Add line (Gain Benign SubRegion) to <LOSS> DD_ADDITIONALREGIONS with FILEANME as key; DGV Supporting Info as VALUE
                                     CALLDATA = "chr"+CHRN +'\t'+ y[0].split('-')[0] +'\t'+ y[0].split('-')[1] +'\t'+ 'CALLED: ' + ';'.join(CALLEDINDEX) + '\n'
                                     DD_ADDITIONALREGIONS[pf.strip().replace('GAIN', 'LOSS')].append(CALLDATA)
                                     break
            return DD_ADDITIONALREGIONS

        if 'LOSS' in FILENAME:
            ctr=0
            for ifp in READPFILE:
                ctr+=1
                DATAIFP = ifp.strip().split('\t')
                ATLEASTLBL = []
                CALLED_COORDL = []
                DD_BEXREGIONS[pf.strip()].append(ifp.strip())
                if DATAIFP[-1].strip().startswith('UNCALLED: '): 
                    CALLEDL = [ lcall.replace('CALLED: ','').split(' ') for lcall in DATAIFP[-2].split(';') ]
                    UNCALLEDL = [ [ luncall.replace('UNCALLED: ','').split(' ')[0].split(':')[1].strip(), luncall.replace('UNCALLED: ','').split(' ')[2].split(':')[1].strip() ] for luncall in DATAIFP[-1].split(';') ]
                    for atl in UNCALLEDL:
                       ATLEASTLBL.extend( i for i in list(DICT_DGV_CHR[atl[0].strip()]) if atl[1].lower() in GAINKEY )
                else:
                    CALLEDL = [ lcall.replace('CALLED: ', '').split(' ') for lcall in DATAIFP[-1].split(';') ]
                    UNCALLED = []
                RESULTL = [ '||'.join(rl).replace('CALLED:||', '') for rl in CALLEDL if float(rl[-1].split(':')[1]) * 100 >= 1 ]
                if len(RESULTL) > 0:
                   LCALLED = []
                   CALLEDIDL = [ rl.split('||')[0].replace('Name:','') for rl in RESULTL ]
                   for rl in RESULTL:
                      for ddcnl in list(DICT_DGV_CHR[filter(None, rl.split('||'))[0].replace('Name:', '').strip()]):
                         COORDL = [int(ddcnl.split(':')[0]), int(ddcnl.split(':')[1])]
                         CALLED_COORDL.append(COORDL)
                         LCALLED.append([rl.split('||')[0].replace('Name:',''), COORDL])
                   if len(CALLED_COORDL) > 1:
                      CHECKOVERLAPL = [ ['_'.join([repr(cco[0]), repr(cco[1])]), '_'.join([repr(rd[0]), repr(rd[1])])] for cco, rd in zip(CALLED_COORDL, CALLED_COORDL[1:]) if len(range(max(cco[0], rd[0]), min(cco[1], rd[1]) + 1)) > 1 ]
                      RESULT_OVERLL = list(set(sum(CHECKOVERLAPL, [])))
                      if len(CHECKOVERLAPL) > 2:
                         S = [ int(RO.split('_')[0]) for RO in RESULT_OVERLL ]
                         E = [ int(RO.split('_')[1]) for RO in RESULT_OVERLL ]
                         ALL = list(set(S+E))
                         ALL.sort(key=int)
                         all_windowsl = []
                         i=0; j=1
                         while i < len(ALL)-1:
                            while j <= len(ALL)-1:
                               all_windowsl.append([int(ALL[i]), int(ALL[j])])
                               i = i + 1; j = j + 1
                         combinel = [ [RO.split('_')[0], RO.split('_')[1]] for RO in RESULT_OVERLL ]
                         targets_countl = []; countertl = 0; DDCBL = defaultdict(list)
                         for each in all_windowsl:
                            compare_benign([repr(each[0]), repr(each[1])], combinel, targets_countl, countertl, DDCBL)
                         for y in targets_countl:
                            if y[1] >= 2* EVENTSGL:
                               for fatl in ATLEASTLBL:
                                  if len(range(int(max(y[0].split('-')[0], fatl.split(':')[0])), int(min(y[0].split('-')[1], fatl.split(':')[1])) + 1)) > 1: 
                                     CALLEDINDEXL = []
                                     for ddcbl in DDCBL[y[0]]: CALLEDINDEXL.extend([ iteml for iteml in dict_all['_'.join(map(str,ddcbl))] if iteml.split(' ')[0].replace('Name:','') in CALLEDIDL ])
                                     CALLDATAL = "chr"+CHRN +'\t'+ y[0].split('-')[0] +'\t'+ y[0].split('-')[1] +'\t'+ 'CALLED: ' + ';'.join(CALLEDINDEXL) + '\n'
                                     DD_ADDITIONALREGIONS[pf.strip().replace('LOSS', 'GAIN')].append(CALLDATAL)
                                     break
            return DD_ADDITIONALREGIONS

 def ARFILES(AR):
    NEWOFILE = open(AR.strip().replace('.txt', '_FINALREGIONS.txt'), 'w') # Open File for Writing
    CHR = AR.strip().split('/')[-1].split('_')[0].lower().replace('x','X').replace('y','Y')
    tree=IntervalTree()

    # Add DD_BEXREGIONS to IntervalTree
    for ddbex in DD_BEXREGIONS[AR.strip()]: tree.add(Interval(int(ddbex.split('\t')[1]),int(ddbex.split('\t')[2])))

    # Add DD_ADDITIONALREGIONS to IntervalTree
    for ddar in DD_ADDITIONALREGIONS[AR.strip()]: tree.add(Interval(int(ddar.split('\t')[1]),int(ddar.split('\t')[2])))

    # Merge Regions
    tree.merge_overlaps(strict=False)

    # Output Final Regions to File
    for item in sorted(tree):
        NEWOFILE.write( '\t'.join([ CHR, repr(item[0]), repr(item[1]) ]) +'\n' )


 # Get list of CHR Folders to Process
 CHRFOLDERS = glob.glob(sys.argv[2].strip() + '/CHR*')
 PROGRESS = Bar('Processing', max=SUM, suffix='%(percent)d%%')
 ctr=0

 for CHR in CHRFOLDERS: 
    CHR = CHR.split('/')[-1]
    COUNT_TOTAL = 0
 
    # Create list of years + 'NA' value supplied to Benign-Ex
    PARAMETERFILES = glob.glob(sys.argv[2].strip() + '/' + CHR + '/*.txt')
    PARAMETERFILES = [ ls for ls in PARAMETERFILES if not 'ADDITIONALREGIONS' in ls if not 'FINALREGIONS' in ls ] # Remove ADDITIONALREGIONS and FINALREGIONS files from list set if they exist
    if len(PARAMETERFILES) == 0: continue

    # Get Year Information
    DICT_DGV = defaultdict(list); dict_all = defaultdict(list)
    YEARS = []
    for i in PARAMETERFILES:
       if i.strip().split('/')[-1].split('_')[-1].replace('.txt','').isdigit(): 
           YEARS.append(i.strip().split('/')[-1].split('_')[-1].replace('.txt',''))
       else: YEARS.append('NA')
    YEARS = list(set(YEARS))

    # Get Raw DGV Data for Current Chromosome    
    d = glob.glob(sys.argv[1].strip() + '/'+ CHR.replace('CHR','chr') + '_*')[0]
    RAWFILE(d)

    # Process Each OUTPUT File (Current Chromosome)
    DD_ADDITIONALREGIONS = defaultdict(list); DD_BEXREGIONS = defaultdict(list)
    for pf in PARAMETERFILES:
        PROCESSINPUT(pf)
        PROGRESS.next()
        ctr+=1
    
    # Write DATA to TMP File
    for addreg in DD_ADDITIONALREGIONS.keys():
        NEWOFILE = open(addreg.replace('.txt', '_ADDITIONALREGIONS.txt'), 'w')
        ADDRESULT = DD_ADDITIONALREGIONS[addreg.strip()]
        for adres in ADDRESULT:
            NEWOFILE.write(adres.strip() + '\n')
        NEWOFILE.close()
        PROGRESS.next()
        ctr+=1

    # Get all Benign Intervals from Original and Additional Regions Files; Merge Benign Intervals
    for pf in PARAMETERFILES:
        ARFILES(pf)
        PROGRESS.next()
        ctr+=1

 PROGRESS.finish()
 print SUM, ctr

except KeyboardInterrupt:
 print ''
 sys.exit(1)