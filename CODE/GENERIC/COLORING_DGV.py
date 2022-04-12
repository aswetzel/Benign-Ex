import glob
import os
import sys
import re

if len(sys.argv) != 4:
    print "COLORING_DGV.py <PARAMETERFILES PATH> DGV <R Command>"
    sys.exit(1)

try:
 PATH=os.path.dirname(os.path.realpath(__file__));
 FILES=glob.glob(PATH.replace('CODE/GENERIC','').strip()+sys.argv[1].strip()+'/'+sys.argv[2].strip()+'_*_DATA'+'/'+'*COMBINED_DATA*');
 RVERSION = sys.argv[3]
 GENOME = re.findall('HG[0-9]{2}', FILES[0])

 # For each COMBINED_DATA File
 # Open HEATMAP_R_CODE.R and Adjust Script Contents Appropriately
 for rf in FILES:
    RCODE=open(PATH+'/HEATMAP_R_CODE_DGV.R','r');
    RCODED=[]; TAG=[];
    EXT = '_'.join(rf.strip().split('/')[-1].split('.')[0].split('_')[3:])

    for coder in RCODE:
        rcode=coder.strip();

        # Change where to read the file from
        if 'COMBINED_DATA <- read.csv' in rcode.strip():
            LINERCODE = rcode.strip().split("'")
            RCODED.append(LINERCODE[0].strip()+"'"+rf.strip()+"'"+LINERCODE[2].strip())

        # Change where to output the MATRIX file
        elif 'write.csv(MATRIX' in rcode.strip():
            LINERCODE = rcode.strip().split("'")
            RCODED.append(LINERCODE[0].strip()+"'"+rf.strip().replace("COMBINED_DATA","MATRIX").replace(".txt",".csv")+"'"+LINERCODE[2].strip())

        # Change when to output the LIST_SCORES file
        elif 'write.csv(BestParam' in rcode.strip():
            LINERCODE = rcode.strip().replace('#','').split("'")
            if rf.strip().split('/')[-1].split('_')[-1] in ('GAIN.txt', 'LOSS.txt'):
                RCODED.append(LINERCODE[0].strip()+"'"+rf.strip().replace("COMBINED_DATA","LIST_SCORES").replace(".txt",".csv")+"'"+LINERCODE[2].strip())
            else: 
                RCODED.append('#'+LINERCODE[0].strip()+"'"+rf.strip().replace("COMBINED_DATA","LIST_SCORES").replace(".txt",".csv")+"'"+LINERCODE[2].strip())

        # Change where to output the heatmap PDF
        elif 'pdf(' in rcode.strip():
            LINERCODE = rcode.strip().split("'")
            RCODED.append(LINERCODE[0].strip()+"'"+rf.strip().replace("COMBINED_DATA",sys.argv[2]).replace(".txt",".pdf")+"'"+LINERCODE[2].strip())

        # Change plot title (DB_GAIN or DB_LOSS)
        elif 'main=' in rcode.strip():
            LINERCODE = rcode.strip().split("'")
 
            if len(rf.strip().split('/')[-1].split('.')[0].split('_')) > 3:
                 RCODED.append(LINERCODE[0].strip()+"'"+sys.argv[2]+'_'+'_'.join(rf.strip().split('/')[-1].split('.')[0].split('_')[-4:])+" ("+GENOME[0]+"): '"+LINERCODE[2].strip())
            else: 
                 RCODED.append(LINERCODE[0].strip()+"'"+sys.argv[2]+' '+rf.strip().split('/')[-1].split('_')[-1].split('.')[0]+" ("+GENOME[0]+"): '"+LINERCODE[2].strip())

        else: RCODED.append(rcode.strip())
    RCODE.close();

    # Write Changes to File
    OFILE=open(PATH+'/HEATMAP_R_CODE_DGV.R','w');
    for rcodl in RCODED:
        OFILE.write(rcodl.strip()+'\n');
    OFILE.close();
    
    # Run HEATMAP_R_CODE.R
    COMMAND = RVERSION +' '+ 'CMD BATCH' +' '+ PATH + '/HEATMAP_R_CODE_DGV.R'
    os.system(COMMAND);

 # Remove .Rout file
 PATH_TO_SCRIPT=os.path.abspath(os.path.realpath(__file__));
 L=PATH_TO_SCRIPT.split('/')[1:];
 KL=0;
 for u in range(KL,len(L)):
    if os.path.exists('/'.join(['']+L[0:KL]).strip()+'/HEATMAP_R_CODE_DGV.Rout'):
        COMMAND='rm '+'/'.join(['']+L[0:KL]).strip()+'/HEATMAP_R_CODE_DGV.Rout'
        os.system(COMMAND);
        break;
    else:
        KL+=1;

except KeyboardInterrupt:
 print ''
 sys.exit(1)
