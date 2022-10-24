###########################################################################################
# Try to find an efficient way to loop over BAM file read and grep only a subset of reads
###########################################################################################

# Libraries
import pysam
import sys

# Help message
if sys.argv[1] in ["-h", "help", "--help", "-help"]:
    print("**************************************")
    print("**************************************")
    print("This script should extract CCS for deepConsensus")
    print("Required inputs are:")
    print("1. Input BAM file with all reads. ")
    print("2. Output BAM file with CCS reads.")
    print("OUTPUT --> reads with NP>=3 AND RQ>=0.88")
    print("**************************************")
    print("**************************************")
else:
    print("**************************************")
    print("**************************************")
    print("This script should extract CCS reads")
    print("Your inputs are:")
    print("1. Input BAM file with all reads --> %s" %(sys.argv[1]))
    print("2. Output BAM file with CCS reads for deepConsensus --> %s" %(sys.argv[2]))
    print("**************************************")
    print("**************************************")
    print("\n\n")
    print("## Analysis started!")

    # Open BAM file
    in_bam = pysam.AlignmentFile(sys.argv[1], "rb", check_sq=False)
    out_bam_deepC = pysam.AlignmentFile(sys.argv[2], "wb", template=in_bam)

    # Create intervals for user update
    print("## Start looping on reads..")
    counter = 0
    # Obtain reads
    for line in in_bam:
        counter += 1
        #print('** %s reads processed' %(counter), end = '\r')
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
            if (info[np_pos][1] >= 3) and (info[rq_pos][1] >= 0.88):
                out_bam_deepC.write(line)
        else:
            print('!!! NP and RQ not found in read --> %s !!!' %(line.query_name))
            break

    in_bam.close()
    out_bam_deepC.close()
    print("## Done!")