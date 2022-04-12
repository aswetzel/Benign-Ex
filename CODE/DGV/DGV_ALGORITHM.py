#change the '/' operator to mean divide by
from __future__ import division 

# load required python modules
import sys, getopt, networkx, re, os
from collections import defaultdict
from intervaltree import Interval, IntervalTree

while len(sys.argv) not in [15, 17]:
    print 'DGV_ALGORITHM.py -i <INFILE> -o <OUTFILE> -e <EVENTS> -n <SAMPLESIZE> -t <TYPE> -f <FREQ> -m <METHODS FILTER> -y <YEAR>'
    sys.exit(2)

try:

 # Create "main" function 
 para = []

 def main(argv):
    inputfile = ''
    outputfile = ''
    EVENTS = ''
    TYPE = ''
    MICROARRAY = ''
    SAMPLESIZE = ''
    FREQ = ''
    YEAR = ''
    try:
        opts, args = getopt.getopt(argv, 'hi:o:e:n:t:m:f:y:', [ # Define the parameter arguments used
         'ifile=',
         'ofile=',
         'EVENTS=',
         'SAMPLESIZE=',
         'TYPE=',
         'MICROARRAY=',
         'FREQ=',
         'YEAR='])
    except getopt.GetoptError: # If a required parameter is [blank] output usage instructions
        print 'DGV_ALGORITHM.py -i <INFILE> -o <OUTFILE> -e <EVENTS> -n <SAMPLESIZE> -t <TYPE> -f <FREQ> -m <METHODS FILTER> -y <YEAR>'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':  # If using the help flag (-h) exit & output usage instructions
            print 'DGV_ALGORITHM.py -i <INFILE> (*) -o <OUTFILE> (*) -e <EVENTS> (#) -n <SAMPLESIZE> (#) -t <TYPE> (#) -f <FREQ> (#) -m <METHODS FILTER> (#) -y <YEAR> (#)'
            print '* Required\t# Optional'
            sys.exit()
        elif opt in ('-i', '--ifile'): # Else ... print each of the arguments to the array 'para'
            inputfile = arg
            para.append(inputfile)
        elif opt in ('-o', '--ofile'):
            outputfile = arg
            para.append(outputfile)
        elif opt in ('-e', '--EVENTS'):
            EVENTS = arg
            para.append(EVENTS)
        elif opt in ('-n', '--SAMPLESIZE'):
            SAMPLESIZE = arg
            para.append(SAMPLESIZE)
        elif opt in ('-t', '--TYPE'):
            TYPE = arg
            para.append(TYPE)
        elif opt in ('-f', '--FREQ'):
            FREQ = arg
            para.append(FREQ)
        elif opt in ('-m', '--MICROARRAY'):
            MICROARRAY = arg
            para.append(MICROARRAY)
        elif opt in ('-y', '--YEAR'):
            YEAR = arg
            para.append(YEAR)
    return para

 # Allows the 'main' function to be imported into a different script.
 if __name__ == '__main__':
    main(sys.argv[1:])

 # Exception Handler; Try to execute block - if Exception thrown print/execute new
 # Institute YEAR Filter 
 iFile=[]
 cg=0
 try:
   para[0] # if input file exists
   try:
      if para[7].isdigit(): # if year filter is a year (XXXX)
       # For each line of file <FLPNE> do:
          # open for reading
          # remove spaces at the beginning and end of string
          # store each line of file as element in list
          # get the year (yr) {specifically from last section of study name}; replace 'b' with 'nothing' {e.g. for 2008b}
             # Manually adjust year for 1000 Genomes/gnomAD
             # if the year is >= {YEAR FILTER} then add whole line to the file "iFile"
        for FLPNE in open(para[0].strip(),'r').readlines():
          cg+=1
          yr=FLPNE.strip().split('\t')[6].strip().split('_')[-1].strip().replace('b','')

          if 'gnomAD' in FLPNE.strip().split('\t')[6]: yr='2019'
          elif '1000_Genomes_Consortium_Pilot_Project' in FLPNE.strip().split('\t')[6]: yr='2010'
          elif '1000_Genomes_Consortium_Phase_1' in FLPNE.strip().split('\t')[6]: yr='2012'
          elif '1000_Genomes_Consortium_Phase_3' in FLPNE.strip().split('\t')[6]: yr='2015'

          if yr.isdigit() and int(yr)>=int(para[7].strip()): iFile.append(FLPNE.strip())

    # If there was no year filter provided; just copy all lines to 'iFile'
      else: iFile=[FLPNE.strip() for FLPNE in open(para[0].strip(),'r').readlines()]

   except IndexError:
     print 'DGV_ALGORITHM.py ERROR ::: Missing YEAR argument. Please refer to the README file.'
     sys.exit()

 # If there was no input file provided.
 except IndexError:
    print 'There is no input file provided. Please refer to the README file.'
    sys.exit()

 #If the there is nothing in 'iFile' exit & output usage instructions
 if len(iFile)<1:
    print '\nERROR: INCORRECT PARAMETER FORMATS IDENTIFIED. PLEASE CHECK PARAMETER INPUT FORMAT IN THE README FILE.'
    sys.exit()
    
 # Institute DGV Sample Level Filtering: { Direction, Sample Size, Frequency, Methods }
 iFile2 = open(para[1].strip(), 'w') # open output file for writing 'iFile2'
 CHR_NAME = open(para[0].strip(),'r').readlines()[0].strip().split('\t')[1].strip(); #open input file for reading; get chrom from first column of row #1
 CHR_NAME=[CHR_NAME if 'chr' in CHR_NAME else 'chr'+CHR_NAME][0]; # add "chr" to chrom if needed
 dict_all = defaultdict(list)
 microarray = ['BAC aCGH', 'SNP array', 'Oligo aCGH', 'Sequencing'] # Methods filter
 list1 = []
 list_ori = []
 list_all = []

 for k in iFile: # for each line <k> of iFile 
    ele = filter(None, k.strip().split('\t')) # split each line put into elements of list
    v = ele[2] + '_' + ele[3]                 # start_stop
    middle = 'Type:' + ele[4] + ' & ' +\
             'Lesion:' + ele[5] + ' & ' +\
             'Study:' + ele[6] + ' & ' +\
             'Method:' + ele[7].replace(' ','_') + ' & ' +\
             'Sample_size:' + ele[8]          # Middle = Type:<> & Lesion:<> & Study:<> & Method:<> & SampleSize: <>

    # Add dictionary entry which appends the variant, type, lesion, study, method, sample size, gain freq, loss freq
    dict_all[v].append('Name:' + ele[0] + ' & ' + middle + ' & ' + 'Freq_Gain:' + repr(int(ele[9]) / int(ele[8])) + ' & ' + 'Freq_Loss:' + repr(int(ele[10]) / int(ele[8])))
    TECHNIQUES = [ ma.strip().replace('"', '').replace('_',' ') for ma in ele[7].strip().split(',') ] # Get Technique(s)
    marray = len(list(set(TECHNIQUES).intersection(set(microarray)))) # Count the number of Techniques that intersect with methods filter :: ['BAC aCGH', 'SNP array', 'Oligo aCGH']
    list_all.append(ele[2] + '-' + ele[3]) # append to FULL list of coordinates

    # If User-Supplied Direction Filter = "GAIN"
        # Calculate Gain Freq
        # If User-Supplied Methods Filter = "Y"
            # If Gain Freq >= User-Supplied Frequency Filter <AND> Sample Size >= User-Supplied Min Sample Size <AND> Methods Filter Count >= 1 THEN
                # Create list1 :: start, stop
                # Create list_ori :: start-stop
        # Else If User-Supplied Methods Filter = "N"
           # If Gain Freq >= User-Supplied Frequency Filter <AND> Sample Size >= User-Supplied Min Sample Size THEN
                # Create list1 :: start, stop
                # Create list_ori :: start-stop

    if para[4].strip() == 'GAIN':
        freq = int(ele[9]) / int(ele[8]) * 100
        if para[6].strip() == 'Y':
            if freq >= float(para[5].strip()) and int(ele[8]) >= int(para[3].strip()) and marray >= 1:
                list1.append([int(ele[2]), int(ele[3])])
                list_ori.append(ele[2].strip() + '-' + ele[3].strip())
        elif para[6].strip() == 'N':
            if freq >= float(para[5].strip()) and int(ele[8]) >= int(para[3].strip()):
                list1.append([int(ele[2]), int(ele[3])])
                list_ori.append(ele[2].strip() + '-' + ele[3].strip())

    # Else If User-Supplied Direction Filter = "LOSS"  {Same Process}
    elif para[4].strip() == 'LOSS':
        freq = int(ele[10]) / int(ele[8]) * 100
        if para[6].strip() == 'Y':
            if freq >= float(para[5].strip()) and int(ele[8]) >= int(para[3].strip()) and marray >= 1:
                list1.append([int(ele[2]), int(ele[3])])
                list_ori.append(ele[2].strip() + '-' + ele[3].strip())
        elif para[6].strip() == 'N':
            if freq >= float(para[5].strip()) and int(ele[8]) >= int(para[3].strip()):
                list1.append([int(ele[2]), int(ele[3])])
                list_ori.append(ele[2].strip() + '-' + ele[3].strip())

 # Look for overlaps with DGV Entries that meet sample-level filtering criteria
 # For each set of coordinates in list1 <h>; IF start != stop THEN add coordinates to interval tree FORMAT :: {Interval(start,stop,start-stop)}
 tree = IntervalTree()
 for h in list1:
    if repr(h[0]).strip() != repr(h[1]).strip():
        tree[int(h[0]):int(h[1])] = repr(h[0]) + '-' + repr(h[1])

 # For each set of coordinates in list1 <j>
    # put coordinates into TMPJ { format :: start-end }
    # create TMPJH and TMPH arrays

    # Loop through the Interval tree (from <h> must be >1bp in length)
        # If the interval tree entry <s> overlaps the set of coordinates <j> { But is not <j> } THEN add to TMPJH array
            # Loop through the TMPJH array <jh> and associate the overlaps with the original entry <j> STORE AS: TMPH
 list2 = []
 for j in list1:
    TMPJ = [repr(j[0]) + '-' + repr(j[1])]
    TMPJH = []
    TMPH = []
    for s in sorted(tree[int(j[0]):int(j[1])]):
        if s.data not in TMPJ: TMPJH.append(s.data.split('-'))
    if len(TMPJH) > 0: 
        for jh in TMPJH: TMPH.append('\t'.join([repr(j[0]), repr(j[1]), jh[0], jh[1]]))
    for tmph in TMPH: list2.append(tmph)

 # For each set of coordinates <q> in list2, create a dictionary entry for <dict1> with <j> as the lookup value and <TMPH> as associations
 dict1 = defaultdict(list)
 for q in list2:
    key1 = q.split('\t')[0] + '_' + q.split('\t')[1]
    dict1[key1].append(q.split('\t')[2] + '-' + q.split('\t')[3])

 # For each dictionary key <j> add the overlapping dgv entries as edges to graph { connect <j> to <TMPH> }
 g = networkx.Graph() #Create empty graph
 for j in dict1.keys():
    for ej in dict1[j.strip()]:
        g.add_edge(j.strip().replace('_', '-'), ej.strip())

 # Determine which dgv entries are connected to one another. { if A::B, A::C, B::D, E::F then {ABCD; EF} } 
 merged_lists = []
 merged_lists = networkx.connected_components(g)
 merge_merged_lists = [val for sublist in merged_lists for val in sublist]
 diff = list(set(list_ori).difference(set(merge_merged_lists))) #Gets a list of all the dgv coordinates not in the network; ie no overlaps.
 for j in diff: g.add_edge(j,j) # add non-overlapping entires to the graph
 merged_lists = networkx.connected_components(g)

 # Create function to create 500bp intervals for given start/stop range
 def range_sequence(start, stop, step):
    result = list(zip(range(start, stop, step), range(start + step, stop+1, step)))
    if (stop - start) % step != 0:
        last_fst_elem = result[-1][-1] if result else start
        result.append((last_fst_elem, stop))
    return result

 # Create function to lookup how many dgv entries overlap the 500bp interval
 def compare_benign(lk, lo, listn, counter):
    region = lk[0] + '-' + lk[1]
    DDCB = defaultdict(list)
    for a in lo:
        if len(range(max(int(lk[0]), int(a[0])), min(int(lk[1]), int(a[1])) + 1)) > 1:
            DDCB[region].append(a)
            counter += 1
    listn.append([region, counter])
    return listn

 # For each set of connected dgv entires, 
    # Split coordinates into start/end coordinate arrays 
    # Combine start/end coordinate arrays to single array <ALL>
    # Sort coordinate array <ALL>
    # Loop through the coordinate array and create a series of start/stop coordinates based on the order of coordinates in the array
        # If the distance between new start/stop is > 500 bp create 500 bp intervals
        # If the distance between new start/stop is < 500 bp use coordinates
    # Find out how many qualifying (Year, N, Method, AF) dgv entries overlap the 500bp interval

 merge_merged_lists = []
 countert = 0
 for o in merged_lists:
    S = []
    E = []
    combine = [ j.split('-') for j in o ]
    S = [ int(I[0]) for I in combine ]
    E = [ int(I[1]) for I in combine ]

    ALL = []
    ALL = S+E
    ALL.sort(key=int)

    targets_count = []
    windows=[]
    all_windows = []
    i=0; j=1

    while i < len(ALL)-1:
        while j <= len(ALL)-1:
            if ALL[j]-ALL[i] > 500: 
                windows = range_sequence(int(ALL[i]), int(ALL[j]), 500)
                all_windows = all_windows + windows

            else: all_windows.append([int(ALL[i]), int(ALL[j])])
            i = i + 1; j = j + 1

    for each in all_windows:
        compare_benign([repr(each[0]), repr(each[1])], combine, targets_count, countert)
        # lk == [repr(each[0]), repr(each[1])] == 500bp window
        # lo == combined == all intervals from connected set of overlapping intervals
        # listn == targets_count == where to save counts
        # counter == countert == reset value to zero for each new 500bp window
    merge_merged_lists.append(targets_count)


 # Create Variables and Define Functions
 tmp4_final = []
 output = []

 def key_func(s):
    return [ (int(x) if x.isdigit() else x) for x in re.findall('\\D+|\\d+', s) ]

 def merge(lol, lol_final):
    lol_final = [ [n, m] for n, m in merge_fast(lol) ]
    return lol_final

 def merge_fast(data):
    data = sorted([ [int(i), int(j)] for i, j in data ])
    it = iter(data)
    a, b = next(it)
    for c, d in it:
        if b + 500 >= c:
            b = max(b, d)
        else:
            yield (a, b)
            a, b = c, d
    yield (a, b)

 # If the number of overlaps exceeds the USER DEFINED MIN # EVENTS { para[2] } THEN
    # Add 500 bp window to dict2
    # If dict2 is NOT empty THEN for each # Events threshold (4,5,6, etc)
        # add 500bp regions to the tmp4 array and sort
        # get MIN start and MAX end as start/end of benign region
 for y in merge_merged_lists:
    tmp4 = []
    dict2 = defaultdict(list)
    for g in y:
        if int(g[1]) >= int(para[2].strip()):
            dict2[g[1]].append(g[0].strip())
    dict2.keys().sort(key=int)

    if len(dict2.keys()) > 0:
        for q in dict2.keys():
            tmp2 = dict2[q]
            tmp4 = tmp4 + tmp2

        sorted_keys = sorted(list(set(tmp4)), key=key_func)
        tmp4 = [ i.split('-') for i in sorted_keys ]
        u = merge(tmp4, tmp4_final)
        if len(u) > 1:
            for e in u: output.append([ int(i) for i in e ])
        else: output.append([ int(i) for i in u[0] ])
    tmp4 = []


 # Define Function
 def compare(lk, lo, listn):
    if len(range(max(lk[0], lo[0]), min(lk[1], lo[1]) + 1)) > 1:
        var = repr(lk[0]) + '\t' + repr(lk[1]) + '\t' + repr(lo[0]) + '\t' + repr(lo[1])
        listn.append(var)
    return listn

 # For each benign region <h>
    # Loop through all dgv entries <[x in list_all]> 
    # If dgv entry overlaps benign region THEN get benign interval; variant interval <tmp4>
    # If dgv entry overlaps benign region THEN lookup variant info in dict_all <tmp>
    # Count the number of acceptable methods for each variant; store count in <m_marray>; 'Method:Merging,SNP array' = 1; Method:Sequencing = 0; Method:Oligo aCGH,SNP array = 2 
        # ^^ For use only if Methods=Y
    # CALLED_INDEX = DGV entries which meet filtering criteria
    # UNCALLED_INDEX = ALL OTHER DGV entries
    # print 'chr' 'benign region' 'CALLED' 'UNCALLED'

 microarray1 = [] # Methods filter
 for a in microarray: microarray1.append(a.replace(' ','_'))

 tempo = list(set(list_all))
 list_all = [ i.split('-') for i in tempo ]
 for h in output:
    tmp4 = []
    for x in list_all:
        f = [ int(i) for i in x ]
        tmp = compare(h, f, tmp4)
    tmp = [ d.split(' & ') for d in sum([ dict_all['_'.join(i.split('\t')[2:])] for i in tmp4 ], []) ] # Get variant info for variants that overlap benign region
    TOTAL = range(0, len(tmp))
    m_marray = [ len(list(set(TMP[4].replace('"', '').split(':')[1].split(',')).intersection(set(microarray1)))) for TMP in tmp ] # Count number of acceptable methods

    # If User-Supplied Direction Filter = "GAIN" and User-Supplied Methods Filter = "Y"
        # Create CALLED_INDEX array
        # If Gain Freq >= User-Supplied Frequency Filter <AND> Sample Size >= User-Supplied Min Sample Size <AND> m_marray count > 1 THEN
            # Add Variant Info to CALLED_INDEX ARRAY
    # If User-Supplied Direction Filter = "GAIN" and User-Supplied Methods Filter = "N"
        # Create CALLED_INDEX array
        # If Gain Freq >= User-Supplied Frequency Filter <AND> Sample Size >= User-Supplied Min Sample Size THEN
            # Add Variant Info to CALLED_INDEX ARRAY
    # SAME PROCESS FOR "LOSS"    

    #print para[4].strip(), para[6].strip()

    if para[4].strip() == 'GAIN' and para[6].strip() == 'Y':
        CALLED_INDEX = [ tmp.index(j) for j in tmp if float(j[-2].split(':')[1]) * 100 >= float(para[5].strip()) and int(j[-3].split(':')[1]) >= int(para[3].strip()) and m_marray[tmp.index(j)] >= 1 ]
    elif para[4].strip() == 'GAIN' and para[6].strip() == 'N':
        CALLED_INDEX = [ tmp.index(j) for j in tmp if float(j[-2].split(':')[1]) * 100 >= float(para[5].strip()) and int(j[-3].split(':')[1]) >= int(para[3].strip()) ]
    elif para[4].strip() == 'LOSS' and para[6].strip() == 'Y':
        CALLED_INDEX = [ tmp.index(j) for j in tmp if float(j[-1].split(':')[1]) * 100 >= float(para[5].strip()) and int(j[-3].split(':')[1]) >= int(para[3].strip()) and m_marray[tmp.index(j)] >= 1 ]
    elif para[4].strip() == 'LOSS' and para[6].strip() == 'N':
        CALLED_INDEX = [ tmp.index(j) for j in tmp if float(j[-1].split(':')[1]) * 100 >= float(para[5].strip()) and int(j[-3].split(':')[1]) >= int(para[3].strip()) ]

    UNCALLED_INDEX = list(set(TOTAL).difference(CALLED_INDEX))   
 
    CALLED = sorted([ ' '.join(tmp[f]) for f in CALLED_INDEX ])
    UNCALLED = sorted([ ' '.join(tmp[f]) for f in UNCALLED_INDEX ])

    if len(UNCALLED) > 0:
        iFile2.write(CHR_NAME + '\t' + '\t'.join([ repr(w) for w in h ]) + '\t' + 'CALLED: ' + ';'.join(CALLED) + '\t' + 'UNCALLED: ' + ';'.join(UNCALLED) + '\n')
    else:
        iFile2.write(CHR_NAME + '\t' + '\t'.join([ repr(w) for w in h ]) + '\t' + 'CALLED: ' + ';'.join(CALLED) + '\n')

 iFile2.close()

 # If File is empty; Delete file
 if os.stat(para[1].strip()).st_size == 0:
    os.remove(para[1].strip())

except KeyboardInterrupt:
 print ''
 sys.exit(1)
