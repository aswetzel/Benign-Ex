# load required python modules
import sys, os
from collections import defaultdict

if len(sys.argv) != 3:
    print "HEATMAP_DATA_DGV.py <PARAMETERFILES PATH> DGV"
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
    PATH_LIST = DATA[0].strip().split('.')[0]
    BEN_LIST = DATA[1].strip()
    STATUS = DATA[2].strip()
    FINAL_METRIC = DATA[-1].strip()
    if STATUS == "GAIN": GAIN[PATH_LIST].append([ BEN_LIST, FINAL_METRIC ])
    if STATUS == "LOSS": LOSS[PATH_LIST].append([ BEN_LIST, FINAL_METRIC ])
 iFile.close()

 # Determine how many heatmaps need to be created
 BEN_SORTED = sorted([ x[0] for x in GAIN[PATH_LIST] ], key=lambda item: (int(item.partition('_')[0]) if item[0].isdigit() else float('inf'), item))
 SS= []; YR=[]; MTHD=[]; MAPS=[]
 for ben in BEN_SORTED:
    DATA = ben.strip().split("_")
    SS.append(ben.strip().split('_')[1])
    MTHD.append(ben.strip().split('_')[3])
    if len(DATA) == 4: YR.append('NA')
    else: YR.append(ben.strip().split('_')[4])

 SS = list(set(SS))
 MTHD = list(set(MTHD))
 YR = list(set(YR))
 for i in SS:
  for j in MTHD:
   for k in YR:
    MAPS.append('_'.join([i,j,k]))

 # Create File for each Pathogenic List with Benign Parameter Settings & AVG(Ranks) Value in appropriate GAIN/LOSS folders
 # Create Master GAIN/LOSS File with all Benign Parameter Settigns + Avg(Ranks) Value Pairs
 ALL_GAIN=open(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/'+'COMBINED_DATA_GAIN.txt',"w");
 ALL_LOSS=open(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/'+'COMBINED_DATA_LOSS.txt',"w");

 # Create Individual COMBINED_DATA_GAIN Files for Each Heatmap
 for j in sorted(GAIN.keys()):
    FILENAME=PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/'+j.strip()+".txt";
    oFile=open(FILENAME,"a");
    if len(MAPS) >1:
     for i in MAPS:
        CAT_GAIN=open(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_GAIN_DATA/'+'COMBINED_DATA_GAIN_'+i.strip()+'.txt',"a");
        for ben in BEN_SORTED:
            DATA = ben.strip().split('_')
            if len(DATA) == 4: ben = '_'.join([ben,'NA'])
            ben2 = '_'.join([ben.strip().split('_')[1],ben.strip().split('_')[3],ben.strip().split('_')[4]])
            if i == ben2:
                for u in GAIN[j.strip()]:
                    if len(u[0].strip().split('_')) == 4: u[0] = '_'.join([u[0],'NA'])
                    u2 = '_'.join([u[0].strip().split('_')[1],u[0].strip().split('_')[3],u[0].strip().split('_')[4]])       
                    if u[0] == ben:
                       oFile.write(u[0]+'\t'+u[1]+'\n');
                       ALL_GAIN.write(j+'\t'+u[0]+'\t'+u[1]+'\n');
                       CAT_GAIN.write(j+'\t'+u[0]+'\t'+u[1]+'\n')
    elif len(MAPS) == 1:
     for i in MAPS:
        for ben in BEN_SORTED:
            DATA = ben.strip().split('_')
            if len(DATA) == 4: ben = '_'.join([ben,'NA'])
            ben2 = '_'.join([ben.strip().split('_')[1],ben.strip().split('_')[3],ben.strip().split('_')[4]])
            if i == ben2:
                for u in GAIN[j.strip()]:
                    if len(u[0].strip().split('_')) == 4: u[0] = '_'.join([u[0],'NA'])
                    u2 = '_'.join([u[0].strip().split('_')[1],u[0].strip().split('_')[3],u[0].strip().split('_')[4]])       
                    if u[0] == ben:
                       oFile.write(u[0]+'\t'+u[1]+'\n');
                       ALL_GAIN.write(j+'\t'+u[0]+'\t'+u[1]+'\n');

 ALL_GAIN.close(); oFile.close();
 try:
  CAT_GAIN.close()
 except NameError:
  pass

 # Create Individual COMBINED_DATA_LOSS Files for Each Heatmap
 for j in sorted(LOSS.keys()):
    FILENAME=PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/'+j.strip()+".txt";
    oFile=open(FILENAME,"a");
    if len(MAPS) >1:
     for i in MAPS:
        CAT_LOSS=open(PATH.replace('CODE/GENERIC','').strip()+'/'+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_LOSS_DATA/'+'COMBINED_DATA_LOSS_'+i.strip()+'.txt',"a");
        for ben in BEN_SORTED:
            DATA = ben.strip().split('_')
            if len(DATA) == 4: ben = '_'.join([ben,'NA'])
            ben2 = '_'.join([ben.strip().split('_')[1],ben.strip().split('_')[3],ben.strip().split('_')[4]])
            if i == ben2:
                for u in LOSS[j.strip()]:
                    if len(u[0].strip().split('_')) == 4: u[0] = '_'.join([u[0],'NA'])
                    u2 = '_'.join([u[0].strip().split('_')[1],u[0].strip().split('_')[3],u[0].strip().split('_')[4]])       
                    if u[0] == ben:
                       oFile.write(u[0]+'\t'+u[1]+'\n');
                       ALL_LOSS.write(j+'\t'+u[0]+'\t'+u[1]+'\n');
                       CAT_LOSS.write(j+'\t'+u[0]+'\t'+u[1]+'\n')
    elif len(MAPS) == 1:
     for i in MAPS:
        for ben in BEN_SORTED:
            DATA = ben.strip().split('_')
            if len(DATA) == 4: ben = '_'.join([ben,'NA'])
            ben2 = '_'.join([ben.strip().split('_')[1],ben.strip().split('_')[3],ben.strip().split('_')[4]])
            if i == ben2:
                for u in LOSS[j.strip()]:
                    if len(u[0].strip().split('_')) == 4: u[0] = '_'.join([u[0],'NA'])
                    u2 = '_'.join([u[0].strip().split('_')[1],u[0].strip().split('_')[3],u[0].strip().split('_')[4]])       
                    if u[0] == ben:
                       oFile.write(u[0]+'\t'+u[1]+'\n');
                       ALL_LOSS.write(j+'\t'+u[0]+'\t'+u[1]+'\n');
 ALL_LOSS.close(); oFile.close();
 try:
  CAT_LOSS.close()
 except NameError:
  pass

except KeyboardInterrupt:
 print ''
 sys.exit(1)
