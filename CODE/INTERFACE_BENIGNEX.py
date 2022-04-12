## load required python modules
import sys, os, glob, re, subprocess, importlib, shutil
from collections import defaultdict

if len(sys.argv) < 2:
    print "\nERROR: AUTOMATION FILE NOT PROVIDED. PLEASE REVIEW README FILE."
    print 'INTERFACE_BENIGNEX.py <AUTOMATION FILE>'
    sys.exit(1)
elif len(sys.argv) != 2:
    print "\nERROR: INCORRECT ARGUMENTS\nINTERFACE_BENIGNEX.py <AUTOMATION FILE PATH>"
    sys.exit(1)

try:
 CWF = os.path.dirname(os.path.realpath(__file__))
 PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
 PATH_TO_PARAMETERS = CWF.replace('CODE', 'PARAMETERS').strip()

 try:
    AUTOMATION = open(sys.argv[1].strip(), 'r').readlines() 
    FOLDERNAME = AUTOMATION[0].strip()
    PYVERSION = AUTOMATION[3].strip()
    RVERSION = AUTOMATION[4].strip()
    OFILELIST = AUTOMATION[5].strip().replace(' ','').split(',')
    if ';' in AUTOMATION[1].strip() or ';' in AUTOMATION[2].strip() or ';' in AUTOMATION[5].strip():
        print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. PLEASE USE \',\' NOT \';\'.'
        sys.exit()
    PATHLISTS = [ i.strip().split('/')[-1] for i in glob.glob(PATH_TO_SCRIPT + '/GENERIC/PATHOGENIC_INTERVALS/*/*')]
    for item in OFILELIST: 
        if item not in ('DEFAULT', 'ALL') and item not in PATHLISTS: 
            print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. \nTHIS PARAMETER: [', item, '] IS INVALID, NEED "DEFAULT" OR "ALL" or PATHOGENIC_LIST as VALUE' 
            sys.exit()
 except IndexError:
    print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. PLEASE REVIEW README FILE.'
    sys.exit()



 # Check if Python Version Requirements (>=2.7.13) Met
 print 'RUNNING A DEPENDENCIES CHECK ...'
 PYTHON_VER = sys.version.split('\n')[0].split(' ')[0].strip()  # Get Python Version
 if '.'.join(PYVERSION.split('.')[0:2]).strip() >= '2.7':
    print 'YOUR PYTHON VERSION INSTALLED IS: ', PYTHON_VER
 else:
    print 'PLEASE INSTALL PYTHON VERSION >=2.7.13'
    sys.exit()

 # Check if R Version Requirments (>= 3.6.0) Met
 R = RVERSION
 R_VER = RVERSION + ' --version'  # Get R Version from AUTOMATION File
 from subprocess import check_call, CalledProcessError, Popen, PIPE

 try:
  FNULL = open(os.devnull, 'w')
  check_call(['which', R],stdout=FNULL, stderr=subprocess.STDOUT)
  try:
    SUBPRC = subprocess.Popen(R_VER, stdout=subprocess.PIPE, shell=True).stdout.read()
    print 'YOUR R VERSION INSTALLED IS: ', SUBPRC.strip().split(' ')[2].strip()
    VER = SUBPRC.strip().split(' ')[2].strip()
    if float(VER[0:3]) >= 3.6: 
        pass
    else:
        print '\nSPECIFIED R VERSION [' + VER + '] INCOMPATIBLE. PLEASE INSTALL R VERSION (>=3.6.0)'
        raise SystemExit
  except subprocess.CalledProcessError:
    print 'COULD NOT CHECK R VERSION... PLEASE INSTALL R'
    sys.exit()
 except CalledProcessError:
    print '\nSPECIFIED R VERSION [' + R + '] NOT INSTALLED. PLEASE INSTALL R VERSION (>=3.6.0)'
    sys.exit()

 # RUN DEPENDENCIES :: Check If Required Python Dependencies Installed; Install if Needed
 COMMAND = PYVERSION + ' ' + CWF.strip() + '/PACKAGES/DEPENDENCIES.py ' + RVERSION
 try:
  if os.system(COMMAND) != 0: raise Exception
 except:
  sys.exit()



 # Create "main" function
 para = []
 def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, 'ha', [ # Define the parameter arguments used
         'ifile='])
    except getopt.GetoptError: # If a required parameter is [blank] output usage instructions
        print 'INTERFACE_BENIGNEX.py <AUTOMATION FILE>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':  # If using the help flag (-h) exit & output usage instructions
            print 'INTERFACE_BENIGNEX.py <AUTOMATION FILE> (*)'
            print '* Required\t# Optional'
            sys.exit()
        elif opt in ('-a', '--auto'): # Else ... print each of the arguments to the array 'para'
            autofile = arg
            para.append(autofile)
    return para


 from progress.bar import Bar
 DIRPARA = glob.glob(PATH_TO_PARAMETERS + '/*_PARAMETERS')
 DEFAULT_PARAVER = defaultdict(list)
 DD_DPARA = defaultdict(list)
 DD_DBVER = defaultdict(set)
 PARAFDIR = []

 # Error Handler: AUTOMATION File Provided? If Not; Exit
 try:
    sys.argv[1].strip()

    # Error Handler: AUTOMATION File Format Correct? If Not; Exit
    try:
        AUTOMATION[2].strip()

        # Store Chosen DATA SOURCE(S) from AUTOMATION (DGV CLINGEN USER) in <DEFAULT_PARAVER>. EX: {'DGV': ['HG19_2020', 'HG18_2016'], 'CLINGEN': ['HG19_2020']}
        # Put all possible PARAMETER FILES in <DD_DPARA>; Database <DGV CLINGEN DEFAULT> as Key; File as Value
        # Get PARAMETER FILE(S) that match Chosen DATA SOURCE(S) <PARAFDIR>
        # Get GENOME, YEAR for DATABASE from AUTOMATION File <DD_DBVER>
        # Store BENIGNEX=Y/N and HEATMAP=Y/N in <APPROVALIST>

        if AUTOMATION[2].strip().startswith('BENIGN-EX'):  
            for deflt in AUTOMATION[1].strip().replace(' ','').split(','): DEFAULT_PARAVER[deflt.strip().split('_')[0].strip()].append('_'.join(deflt.strip().split('_')[1:]).strip())
            for dpara in DIRPARA: DD_DPARA[dpara.strip().split('/')[-1].split('_')[0].strip()].append(dpara.strip())
            for autdb in AUTOMATION[1].replace(' ','').strip().split(','):
                for autop in DD_DPARA[autdb.split('_')[0].strip()]:
                   PARAFDIR.append(autop.strip())
            for ver in AUTOMATION[1].replace(' ','').strip().split(','):
                DD_DBVER[ver.split('_')[0].strip()].add('_'.join(ver.split('_')[1:]).strip()) 
            APPROVALIST = [AUTOMATION[2].strip().split(',')[0].strip().split('=')[1].strip(), AUTOMATION[2].strip().split(',')[1].strip().split('=')[1].strip()]
            if len(APPROVALIST) != 2: 
                print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. PLEASE REVIEW README FILE.'
                sys.exit(1)
            if APPROVALIST[0].strip() == 'N' and APPROVALIST[1].strip() == 'N':
                print '\nERROR: AUTOMATION FILE FORMATION INCORRECT. NOTHING TO RUN.'
                sys.exit(1)
        else:    
            print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. PLEASE REVIEW README FILE.'
            sys.exit()
    except IndexError:
        print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. PLEASE REVIEW README FILE.'
        sys.exit()
 except IndexError:
    print '\nERROR: AUTOMATION FILE NOT PROVIDED. PLEASE REVIEW README FILE.'
    print 'INTERFACE_BENIGNEX.py -a <AUTOMATION FILE>'
    sys.exit()



 # If DGV PARAMETERS: Check PARAMETER Format; print if ERROR in Formatting
 # If CLINGEN PARAMETERS: Check CLINGEN Format; print if ERROR in Formatting
 def CHECK_FORMAT(L, FN):
    PASSED = []
    LEN2 = defaultdict(list)
    ctr=0

    for line in L:
        ctr=ctr+1
        PARMDATA = filter(None, line.strip().replace('EVENTS=', '').split(','))
        LEN2[len(PARMDATA)].append(ctr)

    if len(LEN2.keys()) > 1 and FN.upper().strip() == 'DGV': 
        print "\nERROR: THESE PARAMETER LINE(S) DO NOT HAVE ENOUGH PARAMETERS:", LEN2[min(LEN2.keys())],"\nPLEASE CHECK PARAMETER INPUT FORMAT IN THE README FILE."
        sys.exit(1)

    for line in L:
        PARMDATA = filter(None, line.strip().replace('EVENTS=', '').split(','))
        ctr=ctr+1
        if FN.upper().strip() == 'DGV':
            if len(PARMDATA) >= 5:
                if PARMDATA[0].strip().isdigit():
                    if PARMDATA[1].strip().isdigit():
                        if PARMDATA[2].strip() == 'GAIN' or PARMDATA[2].strip() == 'LOSS':
                            if float(PARMDATA[3].strip()) or PARMDATA[3].strip().isdigit():
                                if PARMDATA[4].strip() == 'Y' or PARMDATA[4].strip() == 'N':
                                    if len(PARMDATA) == 5: PASSED.append(PARMDATA)
                                    elif len(PARMDATA) == 6 and PARMDATA[5].strip() == 'NA' or PARMDATA[5].strip().isdigit(): PASSED.append(PARMDATA)
                                    else: print 'LINE', ctr, 'ERROR: THIS PARAMETER: ', PARMDATA[5].strip(), 'IN', PARMDATA, ' IS INVALID, NEED "NA" OR YEAR VALUE' 
                                else: print 'LINE', ctr, 'ERROR: THIS PARAMETER: ', PARMDATA[4].strip(), 'IN', PARMDATA, ' IS INVALID, NEED "Y" OR "N" VALUE'
                            else: print 'LINE', ctr, 'ERROR: THIS PARAMETER: ', PARMDATA[3].strip(), 'IN', PARMDATA, ' IS INVALID, NEED FLOAT OR INTEGER VALUE'
                        else: print 'LINE', ctr, 'ERROR: THIS PARAMETER: ', PARMDATA[2].strip(), 'IN', PARMDATA, ' IS INVALID, NEED  "GAIN" OR "LOSS" VALUE'
                    else: print 'LINE', ctr, 'ERROR: THIS PARAMETER: ', PARMDATA[1].strip(), 'IN', PARMDATA, ' IS INVALID, NEED AN INTEGER VALUE'
                else: print 'LINE', ctr, 'ERROR: THIS PARAMETER: ', PARMDATA[0].strip(), 'IN', PARMDATA, ' IS INVALID, NEED AN INTEGER VALUE'
            else: print 'LINE', ctr, 'ERROR: THIS PARAMETER LINE DOES NOT HAVE ENOUGH PARAMETERS'

        elif FN.upper().strip() == 'CLINGEN':
            CNVTYPE = ['B', 'USLB', 'VUS', 'B/USLB', 'USLB/VUS', 'B/USLB/VUS']
            if len(PARMDATA) == 2:
                if PARMDATA[0].strip().isdigit():
                    if PARMDATA[1].strip().replace('CNVTYPES=', '') in CNVTYPE:
                        PASSED.append([PARMDATA[0].strip(), PARMDATA[1].strip().replace('CNVTYPES=', '')])
                    else: print 'LINE', ctr, 'ERROR: CNVTYPE "' + PARMDATA[1].strip().replace('CNVTYPES=', '') + '" IS INVALID. ACCEPTABLE VALUES: [B, USLB, VUS, B/USLB, USLB/VUS, B/USLB/VUS].'
                else: print 'LINE', ctr, 'ERROR: THIS PARAMETER: "' + PARMDATA[0].strip() + '" IS INVALID, NEED AN INTEGER VALUE'
            else: print 'LINE', ctr, 'ERROR: THIS PARAMETER LINE [' + ','.join(PARMDATA) + '] HAS AN INCORRECT NUMBER OF PARAMETERS'
    return PASSED


 # Create Sorting Function
 _nsre = re.compile('([0-9]+)')
 def natural_sort_key(s):
    return [ (int(text) if text.isdigit() else text.lower()) for text in re.split(_nsre, s) ]

 # For each PARAMETER File; Get DATABASE; Get PARAMETER SETS
 # Check is PARAMETER Format is CORRECT. Exit if 1+ Errors Found
 # Get List of YEAR/GENOME Combinations for DATABASE
 for file in list(set(PARAFDIR)):
    DBNAME = file.strip().split('/')[-1].split('_')[0].strip()
    PARAMS = []
    with open(file.strip()) as p_file:
        for line in p_file:
            PARAMS.append(line.strip())
    if len(PARAMS) >= 1:
        print ''
        PASSED_PARAMS = CHECK_FORMAT(PARAMS, DBNAME)
        if len(PARAMS) != len(PASSED_PARAMS): 
            print '\nERROR: INCORRECT PARAMETER FORMATS IDENTIFIED. PLEASE CHECK PARAMETER INPUT FORMAT IN THE README FILE.'
            sys.exit(0)
        if FOLDERNAME == 'exit' or FOLDERNAME == 'EXIT' or FOLDERNAME == "" or FOLDERNAME == 'NA' or FOLDERNAME is None:
            print '\nERROR: INVALID FOLDERNAME PROVIDED.'
            sys.exit(0)

        # If DATABASE=DGV & 1+ Correctly Formatted PARAMETERS
        if len(PASSED_PARAMS) > 0 and DBNAME == 'DGV':
            print '**************************'
            print 'PROCESSING DGV DATABASE...'
            print '**************************'
            EV_DGV = [ evdgv.strip() for evdgv in list(DD_DBVER[DBNAME.strip()]) if 'default' not in evdgv.strip() ]

            # If BENIGN-EX=N; HEATMAP=Y
            if len(filter(None, EV_DGV)) > 0 and APPROVALIST[0].strip() == 'N' and APPROVALIST[1].strip() == 'Y':
                for V_DGV in EV_DGV:

                    # Compute OVERLAP COEFFICIENT and CREATE HEATMAP
   
                    try:
                       print 'COMPUTING OVERLAP COEFFICIENT FOR...' + DBNAME
                       print 'THIS MIGHT TAKE TIME, PLEASE BE PATIENT...'
                       COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME + ' ' + V_DGV.strip().split('_')[0].strip()
                       if os.system(COMMAND) != 0: raise Exception
                       COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/RANK_GENOME_ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME
                       if os.system(COMMAND) != 0: raise Exception

                       print '\nIDENTIFYING THE OPTIMAL PARAMETER SET FOR ...' + DBNAME
                       PROGRESS = Bar('Processing', max=2, suffix='%(percent)d%%')
                       COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/HEATMAP_DATA_DGV.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME
                       if os.system(COMMAND) != 0: raise Exception
                       PROGRESS.next()
                       COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/COLORING_DGV.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME + ' ' + RVERSION
                       if os.system(COMMAND) != 0: raise Exception
                    except:
                       sys.exit()
                    PROGRESS.next()
  
                    # Get List of OUTPUT Bed files for each Parameter Combination
                    OFILES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV + '/' + '*.bed')

                    # Read LIST_SCORES File into DICT_SCORES
                    DICT_SCORES=defaultdict(list)
                    SCORES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV + '/' + DBNAME + "*DATA" + "/" + "LIST_SCORES*")
                    for s in SCORES:
                        SS = s.strip().split('/')[-1].split('_')[-1].split('.')[0] 
                        DICT_SCORES[SS.strip()]=defaultdict(list)
                        LIST1=[fi.strip().split('\t') for fi in set(open(s.strip(),"r").readlines())]
                        for lis1 in LIST1: DICT_SCORES[SS.strip()][lis1[0].strip().replace('"','').split(',')[0]].append(lis1[0].strip().replace('"','').split(',')[2])

                        # Get List of PathLists to create 'Optimal' Bed files from
                        for item in OFILELIST:
                            if item in ('DEFAULT', 'ALL'): FINAL = "AVERAGE"
                            else: FINAL = item.strip().split('.')[0]

                            A = '_'.join(DICT_SCORES[SS.strip()][FINAL][0].strip().replace('_NA','').split('_')[0:2])
                            B = '_'.join(DICT_SCORES[SS.strip()][FINAL][0].strip().replace('_NA','').split('_')[2:])
                            C = '_'.join([A,SS,B,'BENIGN-EX'])

                            # Copy the appropriate optimal Gain/Loss file to head OUTPUT folder
                            for f in OFILES:
                                F = f.strip().split('/')[-1]
                                ff = f.strip().replace('BENIGN-EX',FINAL.replace('_','').replace('_','')).replace(DBNAME + '_' + V_DGV + '/','')
                                if C in F: shutil.copy(f, os.path.join(ff))
                    PROGRESS.finish()
                    print '\nYOUR OUTPUT IS SUBMITTED HERE: ', '/'.join(PATH_TO_SCRIPT.split('/')[0:-1]).strip()+ '/' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV

            # If BENIGN-EX=Y; HEATMAP=Y/N
            elif len(filter(None, EV_DGV)) > 0 and APPROVALIST[0].strip() == 'Y' and APPROVALIST[1].strip() in ['N','Y']:

                # Get the list of CHR FILES
                # Get the list of OUTPUT CHR DIRECTORIES & Create FOLDERS (If Don't Already Exist)
                for V_DGV in EV_DGV:
                    GENOME = V_DGV.split('_')[0]
                    TDGVSD = glob.glob(PATH_TO_SCRIPT.replace('CODE', DBNAME + '_SOURCE') + '/' + V_DGV + '/chr*')
                    TDGVSD.sort(key=natural_sort_key)
                    for dgv_inpf in TDGVSD:
                        CHR = dgv_inpf.strip().split('/')[-1].split('_')[0].upper().strip()
                        INPUTFILE = dgv_inpf.strip()
                        OUTDIR = '/'.join(PATH_TO_SCRIPT.split('/')[0:-1]).strip() + '/' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV + '/' + CHR
                        try:
                            os.makedirs(OUTDIR)
                        except OSError:
                            if not os.path.isdir(OUTDIR): raise

                        # Create OUTPUT file
                        # Get parameter arguments for DGV_ALGORITHM
                        # Store DGV_ALGORITHM command in <COMMAND> & RUN DGV_ALGORITHM
                        PROGRESS = Bar('Processing', max=len(PASSED_PARAMS), suffix='%(percent)d%%')
                        for args in PASSED_PARAMS:
                            STRSGRA = ''
                            OUTFILE = OUTDIR + '/' + CHR + '_' + '_'.join(args).strip() + '.txt'
                            for esgra in [sgra for sgra in zip(['-e','-n','-t','-f','-m','-y'],args)]: STRSGRA+=' '+esgra[0].strip()+' '+esgra[1].strip()
                            COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/' + DBNAME + '/' + DBNAME + '_ALGORITHM.py -i ' + INPUTFILE + ' -o ' + OUTFILE + STRSGRA
                            try:
				if os.system(COMMAND) != 0: raise Exception
                            except:
                                sys.exit()
                            PROGRESS.next()
                        print '\t' + CHR.strip()
                        PROGRESS.finish()

                    # Process SUSPICIOUS REGIONS
                    print 'PROCESSING SUSPICIOUS REGIONS ...'
                    COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/' + DBNAME + '/SUSPICIOUS_REGION_DGV_CODE.py ' + PATH_TO_SCRIPT.replace('CODE', DBNAME + '_SOURCE') + '/' + V_DGV + ' ' + '/'.join(PATH_TO_SCRIPT.split('/')[0:-1]).strip() + '/' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV
                    try:
	                if os.system(COMMAND) != 0: raise Exception
                    except:
                        sys.exit()

                    # if HEATMAP=N
                    if APPROVALIST[1].strip() == 'N':
                      try:
                        print '\nCOMPUTING OVERLAP COEFFICIENT FOR...' + DBNAME
                        print 'THIS MIGHT TAKE TIME, PLEASE BE PATIENT...'
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME + ' ' + V_DGV.strip().split('_')[0].strip()
                        if os.system(COMMAND) != 0: raise Exception
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/RANK_GENOME_ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME
                        if os.system(COMMAND) != 0: raise Exception
                      except:
                        sys.exit()

                    # if HEATMAP=Y
                    if APPROVALIST[1].strip() == 'Y':
                      try:
                        print '\nCOMPUTING OVERLAP COEFFICIENT FOR...' + DBNAME
                        print 'THIS MIGHT TAKE TIME, PLEASE BE PATIENT...'
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME + ' ' + V_DGV.strip().split('_')[0].strip()
                        if os.system(COMMAND) != 0: raise Exception
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/RANK_GENOME_ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME
                        if os.system(COMMAND) != 0: raise Exception
                        print '\nIDENTIFYING THE OPTIMAL PARAMETER SET FOR ...' + DBNAME
                        PROGRESS = Bar('Processing', max=2, suffix='%(percent)d%%')
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/HEATMAP_DATA_DGV.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME
                        if os.system(COMMAND) != 0: raise Exception
                        PROGRESS.next()
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/COLORING_DGV.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV +' ' + DBNAME + ' ' + RVERSION
                        if os.system(COMMAND) != 0: raise Exception
                        PROGRESS.next()
                        PROGRESS.finish()
                      except:
                        sys.exit()

                # CREATE BENIGN INTERVAL OUTPUT FILES
                # Loop through each of the CHR Folders
                # For each 'FINALREGIONS.txt' file; add to 'DD_BEX_PARA' dictionary with <PARAMETER> as Key; <FILE PATH> as Value
                print '\nYOUR OUTPUT IS SUBMITTED HERE: ', '/'.join(OUTDIR.split('/')[0:-1]).strip()
                DD_BEX_PARA = defaultdict(list)
                BEXFILES = glob.glob('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/CHR*')
                for bexf in BEXFILES:
                    BEXFDIR = glob.glob(bexf.strip() + '/*FINALREGIONS.txt')
                    for fexd in BEXFDIR:
                        DD_BEX_PARA['_'.join(fexd.strip().replace('_FINALREGIONS.txt', '').split('/')[-1].split('_')[1:]).strip()].append(fexd.strip())

                # Loop through each PARAMETER in DD_BEX_PARA Dictionary; Create PARAMETER Gain/Loss Files; Write Data to 
                PROGRESS = Bar('Processing', max=len(DD_BEX_PARA.keys())+2*len(OFILELIST), suffix='%(percent)d%%')
                for dbpara in DD_BEX_PARA.keys(): 
                    DDBEXOFILEG = defaultdict(list); DDBEXOFILEL = defaultdict(list)
                    BEXOFILEG_BED = open('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/' + '_'.join([FOLDERNAME, DBNAME, V_DGV, dbpara.strip()]).strip() + '_BENIGN-EX.bed', 'w')
                    BEXOFILEL_BED = open('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/' + '_'.join([FOLDERNAME, DBNAME, V_DGV, dbpara.strip()]).strip() + '_BENIGN-EX.bed', 'w')
######              BEXOFILEG = open('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/' + '_'.join([FOLDERNAME, DBNAME, V_DGV, dbpara.strip()]).strip() + '_BENIGN-EX.txt', 'w')
######              BEXOFILEL = open('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/' + '_'.join([FOLDERNAME, DBNAME, V_DGV, dbpara.strip()]).strip() + '_BENIGN-EX.txt', 'w')
                    
                    if 'GAIN' in dbpara.strip():
                        GAIN_SET = []
                        for dexfg in DD_BEX_PARA[dbpara.strip()]:
                            BEXG = open(dexfg.strip(), 'r')
                            for bexg in BEXG:
                                GAIN_SET.append(bexg.strip().split("\t"))
                        for gset in sorted(sorted(GAIN_SET, key=lambda x: int(x[1])), key=lambda x: x[0]):
                            BEXOFILEG_BED.write('\t'.join(gset[0:3]) + '\n')
######                      BEXOFILEG.write('\t'.join(gset) + '\n')
                        PROGRESS.next()
                        BEXG.close()

                    if 'LOSS' in dbpara.strip():
                        LOSS_SET = []
                        for dexfl in DD_BEX_PARA[dbpara.strip()]:
                            BEXL = open(dexfl.strip(), 'r')
                            for bexl in BEXL:
                                LOSS_SET.append(bexl.strip().split("\t"))
                        for lset in sorted(sorted(LOSS_SET, key=lambda x: int(x[1])), key=lambda x: x[0]):
                            BEXOFILEL_BED.write('\t'.join(lset[0:3]) + '\n')
######                      BEXOFILEL.write('\t'.join(gset) + '\n')
                        PROGRESS.next()
                        BEXL.close()
                
                BEXOFILEG_BED.close(); BEXOFILEL_BED.close()
######          BEXOFILEG.close(); BEXOFILEL.close()

                # Get List of OUTPUT Bed files for each Parameter Combination
                OFILES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV + '/' + '*.bed')
                DICT_SCORES=defaultdict(list)
                SCORES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_DGV + '/' + DBNAME + "*DATA" + "/" + "LIST_SCORES*")
                for s in SCORES:
                    SS = s.strip().split('/')[-1].split('_')[-1].split('.')[0] 
                    DICT_SCORES[SS.strip()]=defaultdict(list)
                    LIST1=[fi.strip().split('\t') for fi in set(open(s.strip(),"r").readlines())]
                    for lis1 in LIST1: DICT_SCORES[SS.strip()][lis1[0].strip().replace('"','').split(',')[0]].append(lis1[0].strip().replace('"','').split(',')[2])
                    for item in OFILELIST:
                        if item in ('DEFAULT', 'ALL'): FINAL = "AVERAGE"
                        else: FINAL = item.strip().split('.')[0]
                        A = '_'.join(DICT_SCORES[SS.strip()][FINAL][0].strip().replace('_NA','').split('_')[0:2])
                        B = '_'.join(DICT_SCORES[SS.strip()][FINAL][0].strip().replace('_NA','').split('_')[2:])
                        C = '_'.join([A,SS,B,'BENIGN-EX'])

                        # Copy the appropriate optimal Gain/Loss file to head OUTPUT folder
                        for f in OFILES:
                            F = f.strip().split('/')[-1]
                            ff = f.strip().replace('BENIGN-EX',FINAL.replace('_','').replace('_','')).replace(DBNAME + '_' + V_DGV + '/','')
                            if C in F: 
                                shutil.copy(f, os.path.join(ff))
                                PROGRESS.next()
                PROGRESS.finish()

            else:
                print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. PLEASE REVIEW README FILE.'
                sys.exit()

        # If DATABASE=CLINGEN & 1+ Correctly Formatted PARAMETERS
        elif len(PASSED_PARAMS) > 0 and DBNAME == 'CLINGEN':
            print '\n'
            print '*******************************'
            print 'PROCESSING CLINGEN DATABASE...'
            print '*******************************'
            EV_CLINGEN = [ evclingen.strip() for evclingen in list(DD_DBVER[DBNAME.strip()]) if 'default' not in evclingen.strip() ]

            # If BENIGN-EX=N; HEATMAP=Y
            if len(filter(None, EV_CLINGEN)) > 0 and APPROVALIST[0].strip() == 'N' and APPROVALIST[1].strip() == 'Y':
                for V_CLINGEN in EV_CLINGEN:

                  # Compute OVERLAP COEFFICIENT and CREATE HEATMAP
                  try:
                    print 'COMPUTING OVERLAP COEFFICIENT FOR...' + DBNAME
                    print 'THIS MIGHT TAKE TIME, PLEASE BE PATIENT...'
                    COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME + ' ' + V_CLINGEN.strip().split('_')[0].strip()
                    if os.system(COMMAND) != 0: raise Exception
                    COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/RANK_GENOME_ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME
                    if os.system(COMMAND) != 0: raise Exception

                    print '\nIDENTIFYING THE OPTIMAL PARAMETER SET FOR ...' + DBNAME
                    PROGRESS = Bar('Processing', max=2, suffix='%(percent)d%%')
                    COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/HEATMAP_DATA_CLINGEN.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME
                    if os.system(COMMAND) != 0: raise Exception
                    PROGRESS.next()
                    COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/COLORING_CLINGEN.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME + ' ' + RVERSION
                    if os.system(COMMAND) != 0: raise Exception
                    PROGRESS.next()
                  except:
                    sys.exit()
                  PROGRESS.next()

                  # Create FINAL OUTPUT BED Files(s)
                  OFILES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN + '/' + '*.bed')
                  DICT_SCORES=defaultdict(list)
                  SCORES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN + '/' + DBNAME + "*DATA" + "/" + "LIST_SCORES*")
                  for s in SCORES:
                      SS = s.strip().split('/')[-1].split('_')[-1].split('.')[0] 
                      DICT_SCORES[SS.strip()]=defaultdict(list)
                      LIST1=[fi.strip().split('\t') for fi in set(open(s.strip(),"r").readlines())]
                      for lis1 in LIST1: 
                          DICT_SCORES[SS.strip()][lis1[0].strip().replace('"','').split(',')[0]].append(lis1[0].strip().replace('"','').split(',')[2])
                      for item in OFILELIST:
                          if item in ('DEFAULT', 'ALL'): FINAL = "AVERAGE"
                          else: FINAL = item.strip()
                          A = '_'.join(DICT_SCORES[SS.strip()][FINAL][0].strip().replace('_NA','').split('_')[0:])
                          C = '_'.join([A,SS])
                          for f in OFILES:
                              F = f.strip().split('/')[-1]
                              ff = f.strip().replace('BENIGN-EX',FINAL.replace('_','').replace('.BED','')).replace(DBNAME + '_' + V_CLINGEN + '/','')
                              if C in F: 
                                  shutil.copy(f, os.path.join(ff))
                                  PROGRESS.next()
                  PROGRESS.finish()
                  print '\nYOUR OUTPUT IS SUBMITTED HERE: ', '/'.join(PATH_TO_SCRIPT.split('/')[0:-1]).strip()+ '/' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN

            # If BENIGN-EX=Y; HEATMAP=Y/N
            elif len(filter(None, EV_CLINGEN)) > 0 and APPROVALIST[0].strip() == 'Y' and APPROVALIST[1].strip() in ['N','Y']:

                # Get the list of CHR FILES
                # Get the list of OUTPUT CHR DIRECTORIES & Create FOLDERS (If Don't Already Exist)
                for V_CLINGEN in EV_CLINGEN:
                    GENOME = V_CLINGEN.split('_')[0]
                    TCLINGENSD = glob.glob(PATH_TO_SCRIPT.replace('CODE', DBNAME + '_SOURCE') + '/' + V_CLINGEN + '/chr*')
                    TCLINGENSD.sort(key=natural_sort_key)
                    for clingen_inpf in TCLINGENSD:
                        CHR = clingen_inpf.strip().split('/')[-1].split('_')[0].upper().strip()
                        GAIN_LOSS = clingen_inpf.strip().split('/')[-1].split('_')[-1].replace('.txt', '').strip()
                        INPUTFILE = clingen_inpf.strip()
                        OUTDIR = '/'.join(PATH_TO_SCRIPT.split('/')[0:-1]).strip() + '/' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN + '/' + CHR
                        try:
                            os.makedirs(OUTDIR)
                        except OSError:
                            if not os.path.isdir(OUTDIR): raise

                        # Create OUTPUT file
                        # Get parameter arguments for CLINGEN_ALGORITHM
                        # Store CLINGEN_ALGORITHM command in <COMMAND> & RUN CLINGEN_ALGORITHM
                        PROGRESS = Bar('Processing', max=len(PASSED_PARAMS), suffix='%(percent)d%%')
                        for args in PASSED_PARAMS:
                            OUTFILE = OUTDIR + '/' + CHR + '_' + '_'.join([ sgra.replace('/', '_').strip() for sgra in args ]).strip() + '_' + GAIN_LOSS + '.txt'
                            COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/' + DBNAME + '/' + DBNAME + '_ALGORITHM.py -i ' + INPUTFILE + ' -o ' + OUTFILE + ' -e ' + args[0].strip() + ' -s ' + args[1].strip()
                            try:
				if os.system(COMMAND) != 0: raise Exception
                            except:
                                sys.exit()
                            PROGRESS.next()
                        print '\t' + CHR.strip() +' '+ GAIN_LOSS
                        PROGRESS.finish()
                        
                    # if HEATMAP=N
                    if APPROVALIST[1].strip() == 'N':
                      try:
                        print 'COMPUTING OVERLAP COEFFICIENT FOR...' + DBNAME
                        print 'THIS MIGHT TAKE TIME, PLEASE BE PATIENT...'
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME + ' ' + V_CLINGEN.strip().split('_')[0].strip()
                        if os.system(COMMAND) != 0: raise Exception
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/RANK_GENOME_ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME
                        if os.system(COMMAND) != 0: raise Exception
                      except:
                        sys.exit()
  
                    # if HEATMAP=Y
                    elif APPROVALIST[1].strip() == 'Y':
                      try:
                        print 'COMPUTING OVERLAP COEFFICIENT FOR...' + DBNAME
                        print 'THIS MIGHT TAKE TIME, PLEASE BE PATIENT...'
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME + ' ' + V_CLINGEN.strip().split('_')[0].strip()
                        if os.system(COMMAND) != 0: raise Exception
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/RANK_GENOME_ONE_OVERLAP.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME
                        if os.system(COMMAND) != 0: raise Exception
                        print '\nIDENTIFYING THE OPTIMAL PARAMETER SET FOR ...' + DBNAME
                        PROGRESS = Bar('Processing', max=2, suffix='%(percent)d%%')
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/HEATMAP_DATA_CLINGEN.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME
                        if os.system(COMMAND) != 0: raise Exception
                        PROGRESS.next()
                        COMMAND = PYVERSION + ' ' + PATH_TO_SCRIPT + '/GENERIC/COLORING_CLINGEN.py ' + 'OUTPUT' + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN +' ' + DBNAME + ' ' + RVERSION
                        if os.system(COMMAND) != 0: raise Exception
                        PROGRESS.next()
                        PROGRESS.finish()
                      except:
                        sys.exit()

                # CREATE BENIGN INTERVAL OUTPUT FILES

                # Loop through each of the CHR Folders
                # For each 'FINALREGIONS.txt' file; add to 'DD_BEX_PARA' dictionary with <PARAMETER> as Key; <FILE PATH> as Value

                print '\nYOUR OUTPUT IS SUBMITTED HERE: ', '/'.join(OUTDIR.split('/')[0:-1]).strip()
                DD_BEX_PARA = defaultdict(list)
                BEXFILES = glob.glob('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/CHR*')
                for bexf in BEXFILES:
                    BEXFDIR1 = glob.glob(bexf.strip() + '/*GAIN.txt')
                    BEXFDIR2 = glob.glob(bexf.strip() + '/*LOSS.txt')
                    BEXFDIR = BEXFDIR1 + BEXFDIR2
                    for fexd in BEXFDIR:
                        DD_BEX_PARA['_'.join(fexd.strip().replace('.txt', '').split('/')[-1].split('_')[1:]).strip()].append(fexd.strip())

                # Loop through each PARAMETER in DD_BEX_PARA Dictionary; Create PARAMETER Gain/Loss Files; Write Data to File
                PROGRESS = Bar('Processing', max=len(DD_BEX_PARA.keys())+2*len(OFILELIST), suffix='%(percent)d%%')
                for dbpara in DD_BEX_PARA.keys():  
                    DDBEXOFILEG = defaultdict(list); DDBEXOFILEL = defaultdict(list)
                    BEXOFILEG_BED = open('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/' + '_'.join([FOLDERNAME, DBNAME, V_CLINGEN, dbpara.strip()]).strip() + '_BENIGN-EX.bed', 'w')
                    BEXOFILEL_BED = open('/'.join(OUTDIR.split('/')[0:-1]).strip() + '/' + '_'.join([FOLDERNAME, DBNAME, V_CLINGEN, dbpara.strip()]).strip() + '_BENIGN-EX.bed', 'w')
                    
                    if 'GAIN' in dbpara.strip():
                        GAIN_SET = []
                        for dexfg in DD_BEX_PARA[dbpara.strip()]:
                            BEXG = open(dexfg.strip(), 'r')
                            for bexg in BEXG:
                                GAIN_SET.append(bexg.strip().split("\t"))
                        BEXG.close()
                        for gset in sorted(sorted(GAIN_SET, key=lambda x: int(x[1])), key=lambda x: x[0]):
                            BEXOFILEG_BED.write('\t'.join(gset[0:3]) + '\n')
                        PROGRESS.next()

                    if 'LOSS' in dbpara.strip():
                        LOSS_SET = []
                        for dexfl in DD_BEX_PARA[dbpara.strip()]:
                            BEXL = open(dexfl.strip(), 'r')
                            for bexl in BEXL:
                                LOSS_SET.append(bexl.strip().split("\t"))
                        BEXL.close()
                        for lset in sorted(sorted(LOSS_SET, key=lambda x: int(x[1])), key=lambda x: x[0]):
                            BEXOFILEL_BED.write('\t'.join(lset[0:3]) + '\n')
                        PROGRESS.next()
                    
                    BEXOFILEG_BED.close(); BEXOFILEL_BED.close()

                # Create FINAL OUTPUT BED Files(s)
                OFILES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN + '/' + '*.bed')
                DICT_SCORES=defaultdict(list)
                SCORES = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT') + '/' + FOLDERNAME + '/' + DBNAME + '_' + V_CLINGEN + '/' + DBNAME + "*DATA" + "/" + "LIST_SCORES*")
                for s in SCORES:
                    SS = s.strip().split('/')[-1].split('_')[-1].split('.')[0] 
                    DICT_SCORES[SS.strip()]=defaultdict(list)
                    LIST1=[fi.strip().split('\t') for fi in set(open(s.strip(),"r").readlines())]
                    for lis1 in LIST1: DICT_SCORES[SS.strip()][lis1[0].strip().replace('"','').split(',')[0]].append(lis1[0].strip().replace('"','').split(',')[2])
                    for item in OFILELIST:
                        if item in ('DEFAULT', 'ALL'): FINAL = "AVERAGE"
                        else: FINAL = item.strip()
                        A = '_'.join(DICT_SCORES[SS.strip()][FINAL][0].strip().replace('_NA','').split('_')[0:])
                        C = '_'.join([A,SS])
                        for f in OFILES:
                            F = f.strip().split('/')[-1]
                            ff = f.strip().replace('BENIGN-EX',FINAL.replace('_','').replace('.BED','')).replace(DBNAME + '_' + V_CLINGEN + '/','')
                            if C in F: 
                                shutil.copy(f, os.path.join(ff))
                                PROGRESS.next()
                PROGRESS.finish()

            else:
                print '\nERROR: AUTOMATION FILE FORMAT INCORRECT. PLEASE REVIEW README FILE.'
                sys.exit()

 
 def FINALO(d,STATUS):
    # Open Final OUTPUT File for Writing  
    if d == 'AVERAGE' : NEWOFILE = open(PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/') + 'BENIGN-EX_' + STATUS + '_FINAL.txt', 'w')
    else: NEWOFILE = open(PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/') + 'BENIGN-EX_' + STATUS + '_FINAL_' + d + '.txt', 'w')

    # Loop over Intermediate Files; Merge & Sort Interval; Save to Final OUTPUT File
    DICT_CHR=defaultdict(list)
    READPFILE=open(PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/') + d + '_' + STATUS + '.txt', 'r').readlines()
    for line in READPFILE:
        DATALINE = line.strip().split('\t') # Split data into array
        DICT_CHR[DATALINE[0]].append('-'.join([DATALINE[1],DATALINE[2]]))
    for CHR in sorted(DICT_CHR.keys()): 
        tree=IntervalTree()
        for LINE in DICT_CHR[CHR]: 
            tree.add(Interval(int(LINE.strip().split('-')[0]), int(LINE.strip().split('-')[1]))) # Add Coordinates to IntervalTree
        tree.merge_overlaps(strict=False)
        for item in sorted(tree): NEWOFILE.write( CHR +'\t'+ repr(item[0]) +'\t'+ repr(item[1]) +'\n')

    # Remove Intermediate File
    os.remove(PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/') + d + '_' + STATUS + '.txt')

    # Close File
    NEWOFILE.close()

except:
  print '\n\nERROR: BENIGN-EX RUN CANCELLED.'
  sys.exit()

from intervaltree import Interval, IntervalTree
if APPROVALIST[1].strip() == 'Y':
   print '\n'
   print '********************************************'
   print 'GENERATING FINAL BENIGN-EX BENIGN REGIONS...'
   print '********************************************'
 
   # Remove Prior 'Optimal' BED Files
   REMOVE = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/' + 'BENIGN-EX*FINAL*'))
   for r in REMOVE: os.remove(r)

   # Get List of all 'Optimal' BED Files
   GAIN = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/' + '*GAIN*'))
   LOSS = glob.glob(PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/' + '*LOSS*'))
   ALL = GAIN + LOSS

   # For each 'Optimal' BED File; Extract PathList
   DOUTPUT=[] 
   for i in ALL: DOUTPUT.append(i.strip().split('/')[-1].split('_')[-1].split('.')[0])
   DOUTPUT = list(set(DOUTPUT))
   for d in DOUTPUT:
     gOUTPUT = []; lOUTPUT = []

     # Concatenate then Sort/Merge DGV/CLINGEN GAIN Files; Remove Intermediate File
     for g in GAIN:
        if d in g: gOUTPUT.append(g)
     COMMAND = 'cat ' + ' '.join(map(str, gOUTPUT)) + ' >> ' + PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/') + d + '_' + 'GAIN' + '.txt'
     os.system(COMMAND)
     FINALO(d,'GAIN')

     # Concatenate then Sort/Merge DGV/CLINGEN LOSS Files; Remove Intermediate File
     for l in LOSS:
        if d in l: lOUTPUT.append(l)
     COMMAND = 'cat ' + ' '.join(map(str, lOUTPUT)) + ' >> ' + PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME + '/') + d + '_' + 'LOSS' + '.txt'
     os.system(COMMAND)
     FINALO(d,'LOSS')
  
   # Remove Intermediate Files
   for a in ALL: os.remove(a)
   print 'YOUR OUTPUT IS SUBMITTED HERE: ', PATH_TO_SCRIPT.replace('CODE','OUTPUT' + '/' + FOLDERNAME)

