## load required python modules
from __future__ import division
import os, sys
from collections import defaultdict

if len(sys.argv) != 3:
    print "RANK_GENOME_ONE_OVERLAP.py <PARAMETERFILES PATH> <DGV or CLINGEN>"
    sys.exit(1)

try:
 PATH = os.path.dirname(os.path.realpath(__file__)) # Folder where this file is located
 OFILE = open(PATH.replace('CODE/GENERIC', '').strip() + '/' + sys.argv[1].strip() + '/' + sys.argv[2].strip() + '_RANK_LIST_ONE_OVERLAP_SEQUENCE_BP.txt', 'w') # RANK_LIST_ONE_JACCARD_SEQUENCE_BP File
 IFILE = open(PATH.replace('CODE/GENERIC', '').strip() + '/' + sys.argv[1].strip() + '/' + sys.argv[2].strip() + '_ONE_OVERLAP_BY_CHR.txt') # ONE_OVERLAP File

 dict_LIST = defaultdict(list)
 PARAMTERS_DATA = []
 OFILE.write('PATH_LIST\tPARAMETERS\tDIRECTION\tBENIGN_SC\tPATHOGENIC_SC\tINTERSECTION\tUNION\tMINIMUM\tONE_OVERLAP\tRANK_BENIGN_SC\tRANK_ONE_OVERLAP\tAVG_RANKED_SCORES\n')

 # For each line of <ONE_OVERLAP_BY_CHR>: 
 # Skip First Line; Add to dict_LIST with Pathogenic List Name as KEY; line as VALUE
 for k in IFILE:
    if k.strip().split('vs')[0].strip().split('_')[1] == "PATHOGENIC": next
    else:
        list_name = '_'.join(k.strip().split('vs')[0].strip().split('_')[1:])
        dict_LIST[list_name].append(k.strip())
 IFILE.close()


 # For Chosen PATHOGENIC LIST
 # For each CHR BENIGN PATHOGENIC COMBINATION
    # Isolate PARAMETER Info
    # Add to DICT_PARAMTERS with PARAMETER as Key; LINE as value
 # For each PARAMETER
    # Calculate GENOMEWIDE: BenignSC, PathogenicSC, Intersection, Union, Minimum, 1-OC, 1-OC (path)
    # Add to DICT_PARAMTERS with PARAMETERS as Key; GENOMEWIDE calculations as value
 def CHR_COMBINED(l1):
    DICT_PARAMTERS = defaultdict(list)
    for i in l1:
        PSET = ''
        if sys.argv[2].strip() == 'DGV':
            PSET = '_'.join(i.strip().split('\t')[0].split('vs')[1].replace('_ADDITIONALREGIONS.txt', '').split('_')[1:])
        elif sys.argv[2].strip() == 'CLINGEN':
            PSET = '_'.join(i.strip().split('\t')[0].split('vs')[1].replace('.txt', '').split('_')[1:])
        DICT_PARAMTERS[PSET.strip()].append(i.strip().split('\t')[1:])

    for t in DICT_PARAMTERS:
        SUM_SEQ_CONTENT_B = int(sum([ float(i[0]) for i in DICT_PARAMTERS[t.strip()] ])) # Calculate GENOMEWIDE BenignSC
        SUM_SEQ_CONTENT_P = int(sum([ float(i[1]) for i in DICT_PARAMTERS[t.strip()] ])) # Calculate GENOMEWIDE PathogenicSC
        SUM_INTERSECTION = int(sum([ float(i[2]) for i in DICT_PARAMTERS[t.strip()] ])) # Calculate GENOMEWIDE INTERSECTION
        SUM_UNION = int(sum([ float(i[3]) for i in DICT_PARAMTERS[t.strip()] ])) # Calculate GENOMEWIDE UNION
        SUM_MINIMUM = min(SUM_SEQ_CONTENT_B, SUM_SEQ_CONTENT_P)
        if SUM_MINIMUM > 0: 
            ONE_OVERLAP = 1-(SUM_INTERSECTION / SUM_MINIMUM)
        else: ONE_OVERLAP = 0
        DICT_PARAMTERS[t.strip()] = [repr(SUM_SEQ_CONTENT_B), repr(SUM_SEQ_CONTENT_P), repr(SUM_INTERSECTION), repr(SUM_UNION), repr(min(SUM_SEQ_CONTENT_B,SUM_SEQ_CONTENT_P)), repr(ONE_OVERLAP) ]
    return DICT_PARAMTERS


 # For PATHOGENIC LIST, GENOMEWIDE DATA (For each Parameter Set)
 def COMPUTE_SCORES(d, l1):
    D0 = [ k[0] for k in l1 ] # PARAMETER SET
    D1 = [ int(k[1]) for k in l1 ] # BENIGN SC
    D2 = [ int(k[2]) for k in l1 ] # PATHOGENIC SC
    D3 = [ int(k[3]) for k in l1 ] # INTERSECTION
    D4 = [ int(k[4]) for k in l1 ] # UNION
    D5 = [ int(k[5]) for k in l1 ] # MINIMUM
    D6 = [ float(k[6].strip()) for k in l1 ] # 1-OC
    D7 = [ k[7] for k in l1 ] # DIRECTION

    var = []
    MAX_SEQUENCE_CONTENT = max([ i - min(D1) for i in D1 ])   # MAX VALUE FOR: [BENIGN SC-min(BENIGN SC)]
    if repr(MAX_SEQUENCE_CONTENT) == '0': MAX_SEQUENCE_CONTENT = 1 # If MAX VALUE for: [BENIGN SC-min(BENIGN SC)] = 0; Reset to 1
    DARBRO_RANK_SEQUENCE_CONTENT = [ (i - min(D1)) / MAX_SEQUENCE_CONTENT for i in D1 ]

    MAX_ONE_OC = max([ i - min(D6) for i in D6 ])  # MAX VALUE FOR: [1-OC-min(1-OC)]
    if repr(MAX_ONE_OC) == '0.0': MAX_ONE_OC = 1 # If MAX VALUE for: [1-OC-min(1-OC)] = 0; Reset to 1
    DARBRO_RANK_ONE_OC = [ (i - min(D6)) / MAX_ONE_OC for i in D6 ]

    AVG_RANKED_SCORES = [ (i + j) / 2 for i, j in zip(DARBRO_RANK_ONE_OC, DARBRO_RANK_SEQUENCE_CONTENT) ]

    for d0, d7, d1, d2, d3, d4, d5, d6, A, B, C in zip(D0, D7, D1, D2, D3, D4, D5, D6, DARBRO_RANK_SEQUENCE_CONTENT, DARBRO_RANK_ONE_OC, AVG_RANKED_SCORES):
        var.append(d +'\t'+ d0 +'\t'+ d7 +'\t'+ repr(d1) +'\t'+ repr(d2) +'\t'+ repr(d3) +'\t'+ repr(d4) +'\t'+ repr(d5) +'\t'+ repr(d6) +'\t'+ repr(A) +'\t'+ repr(B) +'\t'+ repr(C))
    return var


 # For each Pathogenic List in dict_LIST
    # Calculate GENOMEWIDE
    # For each PARAMETER SET
        # If PARAMETER contains 'LOSS' or 'DELS': add to <DELS> array
        # If PARAMETER contains 'GAIN' or 'DUPS': add to <DUPS> array
    # If DUPS and DELS are NOT EMPTY: COMPUTE SCORES & Write to File
    # If DUPS or DELS are EMPTY: Compute SCORES & Write to Files for NON_EMPTY
 for e in dict_LIST.keys():
    DELS = []; DUPS = []
    PARAMETERS_DATA = CHR_COMBINED(dict_LIST[e])

    for h in PARAMETERS_DATA.keys():
        if 'LOSS' in h.strip():
            DATA_L = PARAMETERS_DATA[h.strip()]
            DELS.append([h.strip().replace('_LOSS', ''), DATA_L[0], DATA_L[1], DATA_L[2], DATA_L[3], DATA_L[4], DATA_L[5], 'LOSS'])
        elif 'DELS' in h.strip():
            DATA_L = PARAMETERS_DATA[h.strip()]
            DELS.append([h.strip().replace('_DELS', ''), DATA_L[0], DATA_L[1], DATA_L[2], DATA_L[3], DATA_L[4], DATA_L[5], 'LOSS'])
        elif 'DUPS' in h.strip():
            DATA_G = PARAMETERS_DATA[h.strip()]
            DUPS.append([h.strip().replace('_DUPS', ''), DATA_G[0], DATA_G[1], DATA_G[2], DATA_G[3], DATA_G[4], DATA_G[5], 'GAIN'])
        elif 'GAIN' in h.strip():
            DATA_G = PARAMETERS_DATA[h.strip()]
            DUPS.append([h.strip().replace('_GAIN', ''), DATA_G[0], DATA_G[1], DATA_G[2], DATA_G[3], DATA_G[4], DATA_G[5], 'GAIN'])
    if len(DUPS) > 0 and len(DELS) > 0:
        INTRAEVENT_DUPS = sorted(COMPUTE_SCORES(e.strip(), DUPS))
        INTRAEVENT_DELS = sorted(COMPUTE_SCORES(e.strip(), DELS))
        OFILE.write('\n'.join(INTRAEVENT_DUPS).strip() + '\n')
        OFILE.write('\n'.join(INTRAEVENT_DELS).strip() + '\n')
    elif len(DUPS) > 0:
        INTRAEVENT_DUPS = sorted(COMPUTE_SCORES(e.strip(), DUPS))
        OFILE.write('\n'.join(INTRAEVENT_DUPS).strip() + '\n')
    elif len(DELS) > 0:
        INTRAEVENT_DELS = sorted(COMPUTE_SCORES(e.strip(), DELS))
        OFILE.write('\n'.join(INTRAEVENT_DELS).strip() + '\n')
 OFILE.close()

except KeyboardInterrupt:
 print ''
 sys.exit(1)
