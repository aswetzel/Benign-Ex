## load required python modules
from __future__ import division
import sys,os,glob
from progress.bar import Bar
from collections import defaultdict
from intervaltree import Interval, IntervalTree

if len(sys.argv) != 4:
    print "ONE_OVERLAP.py <PARAMETERFILES PATH> <DGV or CLINGEN> <GENOME>"
    sys.exit(1)

try:
 PATH=os.path.dirname(os.path.realpath(__file__)); # Folder where this file is located
 PATHFILES1=glob.glob(PATH+'/PATHOGENIC_INTERVALS/'+sys.argv[3].strip()+'/*.BED'); # List of Pathogenic Intervals
 PATHFILES2=glob.glob(PATH+'/PATHOGENIC_INTERVALS/'+sys.argv[3].strip()+'/*.bed'); # List of Pathogenic Intervals
 PATHFILES = PATHFILES1 + PATHFILES2
 CHRFOLDERS=glob.glob(PATH.replace('CODE/GENERIC',sys.argv[1].strip())+'/CHR*'); # List of CHR Folders
 OFILE=open(PATH.replace('CODE/GENERIC',sys.argv[1].strip())+'/'+sys.argv[2].strip()+"_ONE_OVERLAP_BY_CHR.txt","w"); #OUTPUT FILE: ONE_OVERLAP_BY_CHR.txt

 # For each CHR Folder <i>: 
    # DGV ::: Get the list of files that have 'FINALREGIONS.txt' <DATAFILES>
    # CLINGEN ::: Get the list of files <DATAFILES>
 DATAFILES=[]
 for i in CHRFOLDERS:
    if sys.argv[2].strip()=='DGV':
        FILES = glob.glob(i+'/CHR*FINALREGIONS*.txt')
        for f in FILES[:]:
           if os.path.isfile(f.replace('.txt','.tmp')) == True: FILES.remove(f)
        DATAFILES = DATAFILES + FILES
    elif sys.argv[2].strip()=='CLINGEN':
        FILES=glob.glob(i+'/CHR*.txt')
        DATAFILES = DATAFILES + FILES

 # Total Number of Things to Process:
    # len(DATAFILES) = Number of TOTAL PARAMETER FILES
    # len(PATHFILES)*len(CHRFILES) = Number of total CHRxPATH combinations (24 chromosomes; 1-22,X,Y)
 SUM = len(PATHFILES)*24 + len(DATAFILES) + 2
 PROGRESS = Bar('Processing', max=SUM, suffix='%(percent)d%%')
 ctr=0

 # Function to calculate Benign SC, Pathogenic SC, Union, Intersection, Minimum, Jaccard
    # Open .tmp file for writing
    # For each <Pathogenic Interval Set>: {If 'FINALREGIONS' file is not empty} {If interval >1bp in lenth} Add to PathogenicTREE. Add to UnionTREE
        # For each benign region <t> in DICT_FR: Loop through PTREE. Add to BenignTree. Add to UnionTree. If BENIGN intersects PATHOGENIC add benign to INTN
        # Calculate Intersection
        # Merge Intervals in BTREE, PTREE, UTREE. Calculate Benign SC, Pathogenic SC, Union, Minimum
 def PROCESSINPUT(h):
    OFILEH=open(h.replace('.txt','.tmp'),'w')
    CHR_NUMBER =  h.strip().split('/')[-2]

    for p in [match for match in DICT_PDATA.keys() if CHR_NUMBER+"_" in match]:
        UTREE=IntervalTree()
        BDIFF=[]; PDIFF=[]; UDIFF=[]; IDIFF=[]
        BEN = h.split('/')[-1]
        PATH = p.split('/')[-1]

        BTREE = DICT_BDATA[BEN]
        PTREE = DICT_PDATA[PATH]
        UTREE = BTREE | PTREE
        UTREE.merge_overlaps(strict=False)

        for b in BTREE: BDIFF.append(int(b[1])-int(b[0])) # Benign Interval Length
        for p in PTREE: PDIFF.append(int(p[1])-int(p[0])) # Pathogenic Interval Length
        for u in UTREE: UDIFF.append(int(u[1])-int(u[0])) # Union

        INTN=[]
        ITREE=IntervalTree()
        if len(PTREE) > 0: 
            for btree in BTREE:
                if len(PTREE.overlap(btree[0],btree[1])) > 0:
                    for i in PTREE.overlap(btree[0],btree[1]):
                        ITREE.add(Interval(max(btree[0],i[0]),min(btree[1],i[1])))
            ITREE.merge_overlaps(strict=False) 
            for i in ITREE: IDIFF.append(int(i[1])-int(i[0])) # Intersection
        else: IDIFF=[]

        MINIMUM=[min(sum(BDIFF),sum(PDIFF)) if min(sum(BDIFF),sum(PDIFF)) >0 else 1][0]
        # CHR_PATH_BENIGN + BenignSC + PathSC + Intersection + Union + Minimum
        OFILEH.write(PATH.replace(".BED","").replace(".bed","")+'vs'+BEN.replace("_FINALREGIONS.txt","").replace(".txt","")+\
            '\t'+repr(sum(BDIFF))+\
            '\t'+repr(sum(PDIFF))+\
            '\t'+repr(sum(IDIFF))+\
            '\t'+repr(sum(UDIFF))+\
            '\t'+repr(MINIMUM)+\
            '\n')
    OFILEH.close()

 # For each Pathogenic Inverval Set <e>:
    # Put data into dictionary <DICT_CHRFILES> with filename as key; interval as value
    # Store Intervals in LIST1: Embed dictionary with chr as key; interval as value within DICT_CHRFILES dictionary 
        #{ 'file1': {chr1:[int1, int2,...], chr2[...], ...}, 'file2': {chr1:[int1, int2,...], chr2[...], ...}} 
 DICT_PDATA=defaultdict(list)
 for e in PATHFILES:
    DICT_PFILES=defaultdict(list)
    LIST1=[fi.strip().split('\t') for fi in set(open(e.strip(),"r").readlines())] # Set removes duplicate entries to save time
    for lis1 in LIST1:
        KEY = lis1[0].upper().strip() +"_"+ e.strip().split('/')[-1]
        DICT_PFILES[KEY].append(lis1)        
    for p in DICT_PFILES.keys():
        PTREE = IntervalTree()
        for u in DICT_PFILES[p]: PTREE.add(Interval(int(u[1]),int(u[2])))
        PTREE.merge_overlaps(strict=False)
        DICT_PDATA[p] = PTREE       
        PROGRESS.next()
        ctr+=1
 while ctr < len(PATHFILES)*24: 
    ctr+=1
    PROGRESS.next()

 # For each <DATAFILE> : put data from 'FINALLREGIONS' into dictionary <DICT_FR> with filename as key; 
 # Embed dictionary with chr as key; interval as value within DICT_FR dictionary 
    #{ 'file1': {chr1:[int1, int2,...], chr2[...], ...}, 'file2': {chr1:[int1, int2,...], chr2[...], ...}} 
 DICT_BDATA=defaultdict(list); BTREE=IntervalTree() 
 for j in DATAFILES:
    LIST1=[fi.strip().split('\t') for fi in set(open(j.strip(),"r").readlines())]
    for lis1 in LIST1: BTREE.add(Interval(int(lis1[1]),int(lis1[2])))
    BTREE.merge_overlaps(strict=False)
    DICT_BDATA[j.strip().split('/')[-1]] = BTREE    
    PROCESSINPUT(j.strip())
    DICT_BDATA=defaultdict(list); BTREE=IntervalTree() 
    PROGRESS.next()
    ctr+=1

 # Get List of files with .tmp: Save file contents in dictionary <LKEY> with filename as key; lines as value THEN remove from list of File
 LKEY=defaultdict(list);
 if sys.argv[2].strip()=='DGV':
    for addr in glob.glob(PATH.replace('CODE/GENERIC',sys.argv[1].strip())+'/CHR*'+'/CHR*FINALREGIONS*.tmp'):
        LKEY[addr]=open(addr.strip(),'r').readlines();
        os.remove(addr.strip());
 elif sys.argv[2].strip()=='CLINGEN':
    for addr in glob.glob(PATH.replace('CODE/GENERIC',sys.argv[1].strip())+'/CHR*'+'/*.tmp'):
        LKEY[addr]=open(addr.strip(),'r').readlines();
        os.remove(addr.strip());
 PROGRESS.next()

 # Write <LKEY> Dictionary to file
 OFILE.write('CHR_PATHOGENICvsBENIGN'+'\t'+'BENIGN_SC'+'\t'+'PATHOGENIC_SC'+'\t'+'INTERSECTION'+'\t'+'UNION'+'\t'+'MINIMUM'+'\n')
 for oput in sorted(LKEY.keys()):
     for elkey in LKEY[oput.strip()]:
          OFILE.write(elkey.strip()+'\n');
 OFILE.close();
 PROGRESS.next()
 PROGRESS.finish()

except KeyboardInterrupt:
 print ''
 sys.exit(1)
