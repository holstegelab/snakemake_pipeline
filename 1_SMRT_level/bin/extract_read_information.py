# script to plot read length distribution 

# Libraries
print('** Loading libraries')
import pysam
import pandas as pd
import numpy as np
import sys
import statistics

# Input data
input_ccs = sys.argv[1]
input_nonccs = sys.argv[2]
print('** Input CCS: %s' %(input_ccs))
print('** Input NON-CCS: %s' %(input_nonccs))
in_ccs = pysam.AlignmentFile(input_ccs, "rb", check_sq=False)
in_nonccs = pysam.AlignmentFile(input_nonccs, "rb", check_sq=False)
#in_ccs = pysam.AlignmentFile(in_ccs, "rb", check_sq=False)
#in_nonccs = pysam.AlignmentFile(in_nonccs, "rb", check_sq=False)

# Output data
reads_info = []

# Main
# extract read information for ccs and nonccs
for f in [in_ccs, in_nonccs]:
    for line in f:
        read_name = line.query_name
        read_len = line.query_sequence
        # extract tags
        info = line.tags
        # find tag position for np and rq
        np_pos = 'not_found'
        rq_pos = 'not_found'
        for i in range(len(info)):
            if info[i][0] == 'np':
                np_pos = i
            elif info[i][0] == 'rq':
                rq_pos = i
        # check if read is a ccs of non-ccs
        if (np_pos != 'not_found') and (rq_pos != 'not_found'):
            if f == in_ccs:
                reads_info.append([read_name, len(read_len), info[np_pos][1], info[rq_pos][1], 'ccs'])
            else:
                reads_info.append([read_name, len(read_len), info[np_pos][1], info[rq_pos][1], 'nonccs'])
        else:
            print('!!! NP and RQ not found in read --> %s !!!' %(line.query_name))
            break

# convert list to dataframe
df = pd.DataFrame(reads_info)
df.columns = ['read_name', 'read_length', 'read_pass', 'read_quality', 'read_type']

# write table
if 'ccs.demultiplexed.' in input_ccs:
    tmp_outname = input_ccs.split('.')[0]
    demultiplex_info = input_ccs.split('demultiplexed.')[-1].split('.primrose')[0]
    outname = '%s_%s.ccs.reads_summary.txt' %(tmp_outname, demultiplex_info)
else:
    outname = sys.argv[1].split('.')[0] + '.ccs.reads_summary.txt'
print('** Output summary: %s' %(outname))
df.to_csv(outname, sep = ',')
