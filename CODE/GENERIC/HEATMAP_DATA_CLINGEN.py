# load required python modules
from collections import defaultdict
import re,sys,os;
import glob;

if len(sys.argv) != 3:
    print "HEATMAP_DATA_CLINGEN.py <PARAMETERFILES PATH> CLINGEN"
    sys.exit(1)

try:

 # If GAIN DATA Folder Exists; Deleted it and Make a New One; Else Create GAIN DATA FOLDER
 PATH=os.path.dirname(os.path.realpath(__file__));
 try:
    if os.path.isdir(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/'):
        COMMAND='rm -r '+PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/';
        os.system(COMMAND);    
        os.makedirs(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/');
    else:
        os.makedirs(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/');
 except OSError:
    if not os.path.isdir(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()):
        raise;

 # If LOSS DATA Folder Exists; Deleted it and Make a New One; Else Create LOSS DATA FOLDER
 try:
    if os.path.isdir(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/'):
        COMMAND='rm -r '+PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/';
        os.system(COMMAND);
        os.makedirs(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/');
    else:
        os.makedirs(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/');
 except OSError:
    if not os.path.isdir(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()):
        raise;


 # Loop through DGV_RANK_LIST_ONE_JACCARD_SEQUENCE_BP.txt <iFile>
 # For each line of <iFile> add the final metric value <AVG(RANK BENIGN; RANK 1-OVERLAP)> to the GAIN/LOSS array
 iFile=open(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_RANK_LIST_ONE_OVERLAP_SEQUENCE_BP.txt',"r");
 GAIN=defaultdict(list); LOSS=defaultdict(list); GAIN_DICT=defaultdict(list); LOSS_DICT=defaultdict(list);

 for i in iFile:
    DATA = i.strip().split('\t')
    PATH_LIST = DATA[0].strip()
    BEN_LIST = DATA[1].strip()
    STATUS = DATA[2].strip()
    FINAL_METRIC = DATA[-1].strip()
    if STATUS == "GAIN": GAIN[PATH_LIST].append([ BEN_LIST, FINAL_METRIC ])
    if STATUS == "LOSS": LOSS[PATH_LIST].append([ BEN_LIST, FINAL_METRIC ])
 iFile.close()

 # Create File for each Pathogenic List with Benign Parameter Settings & AVG(Ranks) Value in appropriate GAIN/LOSS folders
 # Create Master GAIN/LOSS File with all Benign Parameter Settigns + Avg(Ranks) Value Pairs
 CAT_GAIN=open(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/'+'COMBINED_DATA_GAIN.txt',"w");
 BEN_SORTED = sorted([ x[0] for x in GAIN[PATH_LIST] ], key=lambda item: (int(item.partition('_')[0]) if item[0].isdigit() else float('inf'), item))

 for j in sorted(GAIN.keys()):
    FILENAME=PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/'+j.strip()+".txt";
    oFile=open(FILENAME,"w");

    for ben in BEN_SORTED:
        for u in GAIN[j.strip()]: 
            if u[0] == ben : 
                oFile.write(u[0]+'\t'+u[1]+'\n');
                CAT_GAIN.write(j+'\t'+u[0]+'\t'+u[1]+'\n');
                GAIN_DICT[ben].append(u[1])
    oFile.close();
 CAT_GAIN.close()

 CAT_LOSS=open(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/'+'COMBINED_DATA_LOSS.txt',"w");
 for j in sorted(LOSS.keys()):
    FILENAME=PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/'+j.strip()+".txt";
    oFile=open(FILENAME,"w");

    for ben in BEN_SORTED:
        for u in LOSS[j.strip()]: 
            if u[0] == ben : 
                oFile.write(u[0]+'\t'+u[1]+'\n');
                CAT_LOSS.write(j+'\t'+u[0]+'\t'+u[1]+'\n');
                LOSS_DICT[ben].append(u[1])
    oFile.close();
 CAT_LOSS.close()

except KeyboardInterrupt:
 print ''
 sys.exit(1)
