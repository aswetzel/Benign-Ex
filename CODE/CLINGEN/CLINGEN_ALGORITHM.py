# load required python modules
import sys, os, getopt
from collections import defaultdict
from intervaltree import Interval, IntervalTree

if len(sys.argv) != 9:
    print 'CLINGEN_ALGORITHM.py -i <INFILE> -o <OUTFILE> -e <EVENTS> -s <STATUS>'
    sys.exit(2)

try:

 # Create "main" function
 para = []
 def main(argv):
    inputfile = ''
    outputfile = ''
    EVENTS=''
    try:
        opts, args = getopt.getopt(argv,"hi:o:e:s:", [ # Define the parameter arguments used
         'ifile=',
         'ofile=',
         'EVENTS=',
         'STATUS='])
    except getopt.GetoptError: # If a required parameter is [blank] output usage instructions
        print 'CLINGEN_ALGORITHM.py -i <INFILE> -o <OUTFILE> -e <EVENTS> -s <STATUS>'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':  # If using the help flag (-h) exit & output usage instructions
            print 'CLINGEN_ALGORITHM.py -i <INFILE> (*) -o <OUTFILE> (*) -e <EVENTS> (#) -s <STATUS> (#)'
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
        elif opt in ('-s', '--STATUS'):
            STATUS = arg
            para.append(STATUS)
    return para

 # Allows the 'main' function to be imported into a different script.
 if __name__ == '__main__':
    main(sys.argv[1:])

 iFile=open(para[0].strip()); # open input file
 iFile2=open(para[1].strip(),"w") # open output file for writing 'iFile2'
 CHR_NAME=''.join(para[0].strip().split('/')[-1].split('_')[0]) # Get Chromosome

 dict_all=defaultdict(list); list_all=[]; list1=[];list_ori=[];
 dict_status={"B":'Benign',"US":'Uncertain_Significance',"USLB":'Uncertain significance: likely benign'};

 for k in iFile: # for each line <k> of iFile
    ele=k.strip().split('\t')           # split each line put into elements of list
    v=ele[1]+'-'+ele[2]                 # start_stop
    dict_all[v].append(ele[3].strip())  # Add dictionary entry with variant name as <key> start_stop as <value>
    list_all.append(ele[1]+'-'+ele[2])  # append to FULL list of coordinates


    # If Acceptable Classifications = "B"; Add benign intervals to list1 and list_ori 
    if(para[3].strip()=="B"):
        cnv_status_s=dict_status[para[3].strip()]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())

    # If Acceptable Classifications = "USLB"; Add benign intervals to list1 and list_ori 
    elif(para[3].strip()=="USLB"):
        cnv_status_s=dict_status["USLB"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())

    # If Acceptable Classifications = "VUS"; Add benign intervals to list1 and list_ori 
    elif(para[3].strip()=="VUS"):
        cnv_status_s=dict_status["US"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())

    # If Acceptable Classifications = "B/USLB"; Add benign/likely benign intervals to list1 and list_ori 
    elif(para[3].strip()=="B/USLB"):
        cnv_status_s=dict_status["B"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())
        cnv_status_s=dict_status["USLB"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())

    # If Acceptable Classifications = "VUS/USLB"; Add benign/likely benign intervals to list1 and list_ori 
    elif(para[3].strip()=="USLB/VUS"):
        cnv_status_s=dict_status["US"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())
        cnv_status_s=dict_status["USLB"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())

    # If Acceptable Classifications = "B/USLB/VUS"; Add benign/likely benign intervals to list1 and list_ori 
    elif(para[3].strip()=="B/USLB/VUS"):
        cnv_status_s=dict_status["B"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())
        cnv_status_s=dict_status["USLB"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())
        cnv_status_s=dict_status["US"]
        if (ele[4].strip()==cnv_status_s):
            list1.append([int(ele[1]),int(ele[2])])
            list_ori.append(ele[1].strip()+'-'+ele[2].strip())

 iFile.close()


 # For each set of coordinates in list1 <h>; Add coordinates to interval tree
 tree=IntervalTree()
 for h in list1: tree.add(Interval(int(h[0]),int(h[1]))) # Add Unique Coordinates to IntervalTree
 tree.merge_overlaps(strict=False)

 # For each clingen variant (that passes filter); compare with interval tree
 # Add to d_OVLP with intervaltree as <KEY>; variant as <VALUE>
 d_OVLP=defaultdict(list)
 for j in list_ori:
    for item in tree[ int(j.split('-')[0]) : int(j.split('-')[1]) ]: 
        COORD = repr(item[0]) +'-'+ repr(item[1])
        d_OVLP[COORD].append(j)

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
        if len(range(max(int(lk[0]), int(a[0])), min(int(lk[1]), int(a[1])) + 1)) > 1: #If the clingen entry overlaps the 500bp interval
            DDCB[region].append(a)
            counter += 1
    listn.append([region, counter])
    return listn

 # For each set of overlapping clingen variants 
    # Split coordinates into start/end coordinate arrays 
    # Combine start/end coordinate arrays to single array <ALL>
    # Sort coordinate array <ALL>
    # Loop through the coordinate array and create a series of start/stop coordinates based on the order of coordinates in the array
        # If the distance between new start/stop is > 500 bp create 500 bp intervals
        # If the distance between new start/stop is < 500 bp use coordinates
    # Find out how many variants overlap the 500bp interval

 # Get all start/end coordinates for the larger interval
 combine = [ j.split('-') for j in list_ori ]
 merge_merged_lists = []
 countert = 0
 for key in d_OVLP:
    S=[]; E=[];

    for each in d_OVLP[key]:
        S.append(each.split('-')[0])
        E.append(each.split('-')[1])
    ALL = []
    ALL = S+E
    ALL = list(set(ALL))
    ALL.sort(key=int)
    
    targets_count = []
    windows=[]
    all_windows = []
    i=0; j=1   

    while i < len(ALL)-1:
        while j <= len(ALL)-1:
            if int(ALL[j])-int(ALL[i]) > 500:
                windows = range_sequence(int(ALL[i]), int(ALL[j]), 500)
                all_windows = all_windows + windows
            elif int(ALL[j])-int(ALL[i]) > 1: all_windows.append([int(ALL[i]), int(ALL[j])])
            i = i + 1; j = j + 1

    for each in all_windows:
        compare_benign([repr(each[0]), repr(each[1])], combine, targets_count, countert)
        # lk == [repr(each[0]), repr(each[1])] == 500bp window
        # lo == combine == all intervals from connected set of overlapping intervals
        # listn == targets_count == where to save counts
        # counter == countert == reset value to zero for each new 500bp window
    merge_merged_lists.append(targets_count)
 flattened = [val for sublist in merge_merged_lists for val in sublist]

 # If the number of overlaps exceeds the USER DEFINED MIN # EVENTS { para[2] } AND the overlap >1bp THEN Add 500 bp window to intervaltree
 tree=IntervalTree()
 for y in flattened: 
    if ( (y[1] >= int(para[2].strip())) and (int(y[0].strip().split("-")[1])-int(y[0].strip().split("-")[0]) > 1) ):
            tree.add(Interval(int(y[0].strip().split("-")[0]), int(y[0].strip().split("-")[1])))
 tree.merge_overlaps(strict=False)

 # Merge benign intervals if <500bp between
 ttree = [ t for t in tree ]
 ttree.sort()
 tree=IntervalTree()
 i=0;j=1
 while i < len(ttree)-1:
    while j <= len(ttree)-1:
        if int(ttree[j][0])-int(ttree[i][1]) < 500: 
            tree.add(Interval(ttree[i][0], ttree[j][1]))
        else: 
            tree.add(Interval(ttree[i][0], ttree[i][1]))
            tree.add(Interval(ttree[j][0], ttree[j][1]))
        i=i+1; j=j+1
 tree.merge_overlaps()

 # Write benign regions to file
 for h in tree:
    iFile2.write(CHR_NAME + '\t' + repr(h[0]) + '\t' + repr(h[1]) + '\n')
 iFile2.close()

except KeyboardInterrupt:
 print ''
 sys.exit(1)
