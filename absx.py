from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from decimal import Decimal
import glob
import os
import shutil
import csv
import sys
import argparse

parser = argparse.ArgumentParser(description='ABSX (Alignment-Based Sequence eXtraction)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-b', metavar='table', help='alignment table file',dest="blastfilearg",required=True)
required.add_argument('-f', choices=['S','M','SM'], help='file / folder mode',dest="filefolder", required=True)
required.add_argument('-x', choices=['n','s','a','b'], help='extraction type: n (whole contig), s (only best hit region), a (extract all hit regions and join them), b (extract region between two outmost hit regions)',dest="extractiontype",required=True)

optional.add_argument('-t', metavar='assembly', help='assembly file',dest="targetf")
optional.add_argument('-q', metavar='query', help='query file(s) to which extracted results are to be appended; if not specified, sequences are extracted into blank files',dest="queryf")
optional.add_argument('-o', metavar='output', help='output folder for modified files with extracted sequences',dest="output", default="absx_out")
optional.add_argument('--om', choices=['query','target', 'combined'], help='output mode: group in files per query [query], per target [target], or combine in a single file [combined]',dest="outM", default="query")
optional.add_argument('-e', metavar='N', help='evalue cutoff',dest="evalue", type=float, default=0.01)
optional.add_argument('-i', metavar='N', help='identity cutoff',dest="identity", type=float, default=0.0)
optional.add_argument('-B', metavar='N', help='bitscore cutoff',dest="bitscore", type=float, default=0.0)
optional.add_argument('-c', metavar='N', help='number of contigs to extract, if set to 0, then extract all contigs',dest="contignum", type=int, default=0)
optional.add_argument('--fl', metavar='N', help='flanks on each side in bp',dest="flanks", type=int, default=0)
optional.add_argument('--lr', dest='local_rec', choices=['none','actual','range'], help='local reciprocator setting', default='range')
optional.add_argument('--is', dest='interstich', action='store_true', help='perform contig stiching', default=False)
optional.add_argument('--translate', dest='trans_out', action='store_true', help='translate output (for -x s or -x a)', default=False)
optional.add_argument('--hit-ovlp', metavar='N', help='allowed hit overlap on query, in bp',dest="hit_ovlp", type=int, default=5)
optional.add_argument('--ctg-ovlp', metavar='N', help='allowed contig overlap on query, in bp',dest="ctg_ovlp", type=int, default=1)
optional.add_argument('--recip-ovlp', metavar='N', help='contig overlap on query for reciprocator selection, in bp',dest="recip_ovlp", type=int, default=10)
optional.add_argument('--bt', choices=['blast','hmmer22', 'hmmer18', 'hmmer15'], help='alignment table type',dest="bt", default="blast")
optional.add_argument('--ac', choices=['dna-dna', 'tdna-aa', 'aa-tdna', 'aa-aa', 'tdna-tdna'], help='alignment coordinate type',dest="ac", default="dna-dna")
optional.add_argument('--acr', choices=['dna-dna', 'tdna-aa', 'aa-tdna', 'aa-aa', 'tdna-tdna'], help='reciprocal alignment coordinate type',dest="acr", default="dna-dna")
optional.add_argument('--acR', choices=['dna-dna', 'tdna-aa', 'aa-tdna', 'aa-aa', 'tdna-tdna'], help='reference alignment coordinate type',dest="acR", default="dna-dna")
optional.add_argument('-r', metavar='file/folder', help='reciprocal search output file or folder',dest="rec_search")
optional.add_argument('-R', metavar='file', help='target locus to reference contig correspondence file',dest="target_ref_file")
optional.add_argument('-m', choices=['e/b-i','e-b-i','b-e-i','i-b-e','i-e-b','b-i-e','e-i-b'], help='order of metrics to use to select best matches (e - evalue, b - bitscore, i - identity)',dest="metric", default="e/b-i")
optional.add_argument('--rescale-metric', dest='metricR', action='store_true', help='divide metric value by length of hit region', default=False)
optional.add_argument('--metric-merge-corr', metavar='N', help='reduce combined metric by this value',dest="metricC", type=float, default=0.75)
optional.add_argument('--no-hs', dest='no_hs', action='store_true', help='do not run hit sticher', default=False)
optional.add_argument('--ref-hs', dest='ref_hs', action='store_true', help='run hit sticher on reciprocal table (slow)', default=False)
optional.add_argument('--keep-strand', dest='keep_strand', action='store_true', help='keep original contig direction', default=False)
optional.add_argument('--rm-rec-not-found', dest='rmrecnf', action='store_true', help='remove hits without matches in reciprocal search', default=False)
optional.add_argument('--max-gap', metavar='N', help='max gap between hits on either query or target, use 0 for no filtering',dest="max_gap", type=int, default=0)
#add possibility of single blast file and multiple fasta files
#tied to that is extension in the query name (*.fas)
#modify --lr to allow literally one query per contig. check that it works correctly with -x n etc...
#bed parser for reference

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()
else:
    args = parser.parse_args()

    blastfilearg = vars(args)["blastfilearg"]
    filefolder = vars(args)["filefolder"]
    bt = vars(args)["bt"]
    if bt == "hmmer18":
        if vars(args)["extractiontype"] != "n":
            print "only whole contig extraction is supported for hmmer18 tables, option -x ignored"
        if vars(args)["ac"] != "aa-aa":
            print "hmmer18 tables are only for AA vs AA search, option --ac ignored"
        extractiontype = "n"
        ac = "aa-aa"
    elif bt == "hmmer15":
        if vars(args)["ac"] != "dna-dna":
            print "hmmer15 tables are only for DNA vs DNA search, option --ac ignored"
        ac = "dna-dna"
        extractiontype = vars(args)["extractiontype"]
    else:
        ac = vars(args)["ac"]
        extractiontype = vars(args)["extractiontype"]
    if vars(args)["no_hs"] or extractiontype == "n":
        print "without hit sticher, contig stiching is disabled, option --is ignored"
        run_hs = False
        interstich = False
    else:
        run_hs = True
        interstich = vars(args)["interstich"]
    if extractiontype == "n" and vars(args)["keep_strand"]:
        keep_strand = True
    else:
        keep_strand = False
    if vars(args)["trans_out"]:
        if extractiontype == "s" or extractiontype == "a":
            if ac == "aa-aa" or ac == "tdna-aa":
                print "should not attempt to translate AA sequence, option --translate ignored"
                trans_out = False
            else:
                trans_out = True
        else:
            print "will only translate exact matches (options -x s or -x a), option --translate ignored"
            trans_out = False
    else:
        trans_out = False

    if vars(args)["targetf"] == None:
        dry_run = True
        noq = True
        print "option -q ignored"
    else:
        targetf = vars(args)["targetf"]
        dry_run = False
        outM = vars(args)["outM"]
        if outM == "query":
            if vars(args)["queryf"] == None:
                noq = True
            else:
                queryf = vars(args)["queryf"]
                noq = False
        else:
            noq = True
    
    evalue = vars(args)["evalue"]
    bitscore = vars(args)["bitscore"]
    identity = vars(args)["identity"]
    contignum = vars(args)["contignum"]
    local_rec = vars(args)["local_rec"]
    if vars(args)["rec_search"] == None:
        rec_search = None
    else:
        rec_search = vars(args)["rec_search"]
        if vars(args)["target_ref_file"] == None:
            print "please specify -R"
            sys.exit()
        else:
            target_ref_file = vars(args)["target_ref_file"]
        acr = vars(args)["acr"]
        acR = vars(args)["acR"]
    rmrecnf = vars(args)["rmrecnf"]        

    hit_ovlp = vars(args)["hit_ovlp"]
    ctg_ovlp = vars(args)["ctg_ovlp"]
    recip_ovlp = vars(args)["recip_ovlp"]
    flanks = vars(args)["flanks"]
    output_dir = vars(args)["output"]
    if vars(args)["metric"] == "e/b-i":
        metric = "e/b-i"
    else:
        metric = [int(x) for x in vars(args)["metric"].replace("e","0").replace("b","1").replace("i","2").split("-")]
    metricR = vars(args)["metricR"]
    metricC = vars(args)["metricC"]
    ref_hs = vars(args)["ref_hs"]
    max_gap = vars(args)["max_gap"]


dashb = "#"*75
dash = "-"*50
warninglist = []

#function to display messages and write them into the debugfile
def messagefunc(msg, columns, f, fl=True):
    if fl:
        sys.stdout.write((" "*columns)+"\r")
        sys.stdout.write(msg[:columns]+"\r")
        sys.stdout.flush()
    else:
        sys.stdout.write((" "*columns)+"\r")
        print msg
    print >> f, msg


#function for creating an output folder. old stuff will be deleted
def mkdirfunc(dir1):
    if not os.path.exists (dir1):
        os.makedirs(dir1) #creating folder if necessary
    else:
        shutil.rmtree(dir1) #removing old files
        os.makedirs(dir1)

#function for copying alignment files before changing them
#all alignments are copied regardless of whether they will be modified or not
def copyfunc(dir1, cols, debugfile):
    messagefunc("copying files *.fa*", cols, debugfile, False)
    copyfunc_c = 0
    for x in glob.glob(queryf+"/*.fa*"):
        locusfname = x.split("/")[-1]
        #print locusfname
        if not os.path.exists (dir1+"/"+locusfname):
            prog = "copying "+str(locusfname)+"..."
            messagefunc(prog, cols, debugfile)
            shutil.copy2(queryf+"/"+locusfname, dir1)
            copyfunc_c += 1
    messagefunc("copied "+str(copyfunc_c)+" files", cols, debugfile, False)

#function for parsing a blast output file
def readblastfilefunc(b, evalue1, bitscore1, identity1, as_target, ac3, recstats, cols, debugfile):
    messagefunc("processing "+b, cols, debugfile, False)
    returndict = {}
    blastfile = open(b, "rU")
    reader = csv.reader(blastfile, delimiter='\t')
    linecounter = 0
    for row in reader:
        if evalue1: 
            if float(row[10]) <= evalue1 and float(row[11]) >= bitscore1 and float(row[2]) >= identity1:
                qname = row[0].split("/")[-1] #also allow *.fas
                tname = row[1]
                if recstats:
                    init_queries.add(qname)
                    init_targets.add(tname)
                if as_target:
                    #populate target table
                    if tname in returndict:
                        if qname in returndict[tname]:
                            returndict[tname][qname][linecounter] = rowfunc(row, ac3)
                        else:
                            returndict[tname][qname] = {linecounter: rowfunc(row, ac3)}
                    else:
                        returndict[tname] = {qname: {linecounter: rowfunc(row, ac3)}}
                else:
                    #populate query table
                    if qname in returndict:
                        if tname in returndict[qname]:
                            returndict[qname][tname][linecounter] = rowfunc(row, ac3)
                        else:
                            returndict[qname][tname] = {linecounter: rowfunc(row, ac3)}
                    else:
                        returndict[qname] = {tname: {linecounter: rowfunc(row, ac3)}}
            linecounter += 1
        else:
            qname = row[0].split("/")[-1]
            tname = row[1]
            if as_target:
                #populate target table
                if tname in returndict:
                    if qname in returndict[tname]:
                        returndict[tname][qname][linecounter] = rowfunc(row, ac3)
                    else:
                        returndict[tname][qname] = {linecounter: rowfunc(row, ac3)}
                else:
                    returndict[tname] = {qname: {linecounter: rowfunc(row, ac3)}}
            else:
                #populate query table
                if qname in returndict:
                    if tname in returndict[qname]:
                        returndict[qname][tname][linecounter] = rowfunc(row, ac3)
                    else:
                        returndict[qname][tname] = {linecounter: rowfunc(row, ac3)}
                else:
                    returndict[qname] = {tname: {linecounter: rowfunc(row, ac3)}}
            linecounter += 1
    blastfile.close()
    return returndict

#function to return a list for dict like this:
#dict[query] = [target_f, target_r, target_b, query_f, query_r, query_b, eval, bitscore, identity]
#query - query name, target_f - target start pos, target_r - target end, target_b - forward or reverse target direction
#query_f - query start, query_r - query end, query_b - query direction

def rowfunc(row, aligntype):
    if int(row[8]) < int(row[9]):
        target_b = True
    else:
        target_b = False
    if aligntype == "tdna-aa":
        if target_b:
            target_f = int(row[8])*3-2
            target_r = int(row[9])*3
        else:
            target_f = int(row[8])*3
            target_r = int(row[9])*3-2
    else:
        target_f = int(row[8])
        target_r = int(row[9])
    #check query
    if int(row[6]) < int(row[7]):
        query_b = True
    else:
        query_b = False
    if aligntype == "aa-tdna":
        if query_b:
            query_f = int(row[6])*3-2
            query_r = int(row[7])*3
        else:
            query_f = int(row[6])*3
            query_r = int(row[7])*3-2
    else:
        query_f = int(row[6])
        query_r = int(row[7])
    return [target_f, target_r, target_b, query_f, query_r, query_b, float(row[10]), float(row[11]),float(row[2])]


#function for parsing hmmer output tables
def readhmmerfilefunc(b, evalue1, bitscore1, bt1, ac1, cols, debugfile):
    messagefunc("processing "+b, cols, debugfile, False)
    targetdict = {}
    hmmfile = open(b, "rU")
    linecounter = 0
    recordcounter = 0
    if bt1 == "hmmer18":
        qname1 = 2
        tname1 = 0
        eval1 = 4
        bit1 = 5
    elif bt1 == "hmmer15":
        qname1 = 2
        query1 = 4
        query2 = 5
        tname1 = 0
        target1 = 6
        target2 = 7
        eval1 = 12
        bit1 = 13
    elif bt1 == "hmmer22":
        qname1 = 3
        query1 = 15
        query2 = 16
        tname1 = 0
        target1 = 17
        target2 = 18
        eval1 = 12
        bit1 = 13
    for row in hmmfile:
        if row[0] != "#":
            line = row.strip().split()
            if float(line[eval1]) <= evalue1 and float(line[bit1]) >= bitscore1:
                qname = line[qname1]+".fas"
                if ac1 == "aa-tdna":
                    tname = "_".join(line[tname1].split("_")[:-1])
                    strand = line[tname1].split("_")[-1][0]
                    frame = int(line[tname1].split("_")[-1][1])
                    ctg_length = int(line[tname1].split("_")[-1][3:])
                else:
                    tname = line[tname1]
                init_queries.add(qname)
                init_targets.add(tname)
                if bt1 != "hmmer18":
                    if int(line[target1]) < int(line[target2]):
                        target_b = True
                    else:
                        target_b = False
                    if int(line[query1]) < int(line[query2]):
                        query_b = True
                    else:
                        query_b = False
                    target_f = int(line[target1])
                    target_r = int(line[target2])
                    query_f = int(line[query1])
                    query_r = int(line[query2])
                else:
                    target_f = 0
                    target_r = 1
                    target_b = True
                    query_f = 0
                    query_r = 1
                    query_b = True
                #no identity is reported, so it is set to 0.0 and does not interfere with scoring system
                if ac1 == "aa-tdna":
                    query_f = int(line[query1])*3-2
                    query_r = int(line[query2])*3
                    if strand == "f":
                        target_b = True
                    else:
                        target_b = False
                    if target_b:
                        target_f = int(line[target1])*3-2+frame
                        target_r = int(line[target2])*3+frame
                        #test this adjustment
                        if target_r > ctg_length:
                            target_r = ctg_length
                    else:
                        target_f = ctg_length-int(line[target1])*3+frame
                        #test this adjustment
                        if target_f < 0:
                            target_f = 0
                        target_r = ctg_length-int(line[target2])*3-2+frame
                    outrow = [target_f, target_r, target_b, query_f, query_r, query_b, float(line[eval1]), float(line[bit1]), 0.0]
                    # print >> debugfile, outrow
                    # print >> debugfile, line
                    # sys.exit()
                else:
                    outrow = [target_f, target_r, target_b, query_f, query_r, query_b, float(line[eval1]), float(line[bit1]), 0.0]
                #populate target table
                if tname in targetdict:
                    if qname in targetdict[tname]:
                        targetdict[tname][qname][linecounter] = outrow
                    else:
                        targetdict[tname][qname] = {linecounter: outrow}
                else:
                    targetdict[tname] = {qname: {linecounter: outrow}}
            linecounter += 1
    hmmfile.close()
    return targetdict

#first main function
def target_processor(inpdict, local_rec, metric, metricR, hit_overlap, recip_overlap, ac1, run_hs1, max_gap1, cols, debugfile):
    messagefunc(dashb, cols, debugfile)
    messagefunc("running target processor...", cols, debugfile)
    outdict = {}
    for targetkey, targetval in inpdict.items():
        messagefunc(dash, cols, debugfile)
        messagefunc("target processor on target "+targetkey, cols, debugfile)
        #run local actual reciprocal check (check each hit of target only matches one query)
        if local_rec == "actual":
            if len(targetval) > 1:
                messagefunc("running actual (per hit) reciprocity check", cols, debugfile)
                targetval = actual_reciprocator(targetval,recip_overlap, metric, metricR)
            else:
                messagefunc("only one hit for this target, no reciprocity check", cols, debugfile)
        #run hit overlapper and sticher 
        messagefunc("running hit processor", cols, debugfile)
        if run_hs1:
            tgt_proc_out = hit_sticher(targetval, metric, metricR, hit_overlap, ac1, max_gap1, cols, debugfile) #stiched subcontigs per query
        else:
            tgt_proc_out = reformat_hits(targetval, metric, metricR, cols, debugfile)
        #run local range reciprocal check
        if local_rec == "range":
            if len(tgt_proc_out) > 1:
                messagefunc("running range reciprocity check", cols, debugfile)
                tgt_proc_out = range_reciprocator(targetkey, tgt_proc_out, metric, metricR, recip_overlap, cols, debugfile) #check that each subcontig matches only 1 Q
            else:
                messagefunc("only one query for this target, no reciprocity check", cols, debugfile)
        outdict[targetkey] = tgt_proc_out
    return outdict


#function to split hits by strand (same vs opposite)
#input is dictionary {linecounter: [target_f, target_r, target_b, query_f, query_r, query_b, float(row[10]), float(row[11]),float(row[2])]}
#return list of dictionaries
def strand_selector(inpdict):
    #split into things per direction
    clusterF = {}
    clusterR = {}
    for hitkey, hitval in inpdict.items():
        if hitval[2] == hitval[5]:
            clusterF[hitkey] = hitval
        else:
            clusterR[hitkey] = hitval
    return [clusterF, clusterR]
    
#function to check that each target hit region matches to only one query
# !!! implement direction check for translated searches
def actual_reciprocator(inpdict, recip_overlap, metric, metricR):    
    returndict = inpdict
    for refquerykey, refqueryval in inpdict.items():
        for refhit in refqueryval.keys():
            refhit_t_range = refqueryval[refhit][0:2]
            refhit_score = refqueryval[refhit][6:9]
            for testquerykey, testqueryval in inpdict.items():
                if testquerykey != refquerykey:
                    for testhit in testqueryval.keys():
                        testhit_t_range = testqueryval[testhit][0:2]
                        testhit_score = testqueryval[testhit][6:9]
                        if getOverlap(refhit_t_range, testhit_t_range) > recip_overlap:
                            comp1 = compare_scores(refhit_t_range, refhit_score, testhit_t_range, testhit_score, metric, metricR)
                            if comp1 == 0:
                                del returndict[testquerykey][testhit]
    return returndict


#helper function for hit sticher
#get hit dictionary, bucket dictionary, and bucket counter, process and return these variables
def sticher_helper(refhitkey1, hitkey1, hitdict1, bucketdict1, nbuck1, opt1, cols1, debugfile1):
    if opt1 == "add":
        #both are in hitdict, however may not be in bucketdict
        #keep both, asign different buckets
        if refhitkey1 not in bucketdict1: #otherwise it must have already gotten the bucket number
            bucketdict1[refhitkey1] = nbuck1
            nbuck1 += 1
        if hitkey1 not in bucketdict1: #same here
            bucketdict1[hitkey1] = nbuck1
            nbuck1 += 1
    elif opt1 == "samebuck":
        #both are in hitdict, if in bucketdict - assign same.
        if refhitkey1 in bucketdict1: #if one present, the other must have not been checked before
            bucketdict1[hitkey1] = bucketdict1[refhitkey1] #put in ref bucket
        elif hitkey1 in bucketdict1:
            bucketdict1[refhitkey1] = bucketdict1[hitkey1]
        else:
            bucketdict1[refhitkey1] = nbuck1
            bucketdict1[hitkey1] = nbuck1
        #no need to up the current bucket since only old buckets were worked with
    elif opt1 == "rmhit":
        if refhitkey1 not in bucketdict1: #otherwise it must have already gotten the bucket number
            bucketdict1[refhitkey1] = nbuck1
            nbuck1 += 1
        del hitdict1[hitkey1]
        if hitkey1 in bucketdict1:
            del bucketdict1[hitkey1]
    elif opt1 == "rmref":
        if hitkey1 not in bucketdict1: #otherwise it must have already gotten the bucket number
            bucketdict1[hitkey1] = nbuck1
            nbuck1 += 1
        del hitdict1[refhitkey1]
        if refhitkey1 in bucketdict1:
            del bucketdict1[refhitkey1]
    return hitdict1, bucketdict1, nbuck1

def reformat_hits(inpdict, metric, metricR, cols, debugfile):
    returnlist = {} #all queries for the target go here
    for querykey, queryval in inpdict.items():
        clusters = strand_selector(queryval)
        directions = [True, False]
        best_index1 = None
        best_dir1 = None
        best_range1 = None
        best_score1 = None
        best_data1 = None
        for cluster_index in range(2):
            direct = directions[cluster_index]
            hits1 = clusters[cluster_index].values()
            for hit_index in range(len(hits1)):
                if best_index1 == None:
                    best_index1 = hit_index
                    best_dir1 = [direct]
                    best_range1 = hits1[hit_index][0:2]+hits1[hit_index][3:5]
                    best_score1 = hits1[hit_index][6:9]
                    best_data1 = best_range1+best_score1
                else:
                    hit_scores = hits1[hit_index][6:9]
                    hit_ranges = hits1[hit_index][0:2]+hits1[hit_index][3:5]
                    hit_data = hit_ranges+hit_scores
                    comp1 = compare_scores(best_range1, best_score1, hit_ranges, hit_scores, metric, metricR)
                    if comp1 == 1:
                        best_index1 = hit_index
                        best_dir1 = [direct]
                        best_range1 = hits1[hit_index][0:2]+hits1[hit_index][3:5]
                        best_score1 = hits1[hit_index][6:9]
                        best_data1 = best_range1+best_score1
        returnlist[indexer_function(querykey,str(0))] = [best_dir1, best_score1, best_range1, best_data1]
        print >> debugfile, querykey, ", selected best hit:", returnlist[indexer_function(querykey,str(0))]
    return returnlist


#function to process hit of the same target (remove redundant, order and stich hits, make subcontigs)
def hit_sticher(inpdict, metric, metricR, hit_overlap, ac2, max_gap2, cols, debugfile):
    returnlist = {} #all queries for the target go here
    for querykey, queryval in inpdict.items():
        clusters = strand_selector(queryval)
        directions = [True, False]
        stiched_subcontigs = [] #stiched subcontigs for a query go here
        for cluster_index in range(2):
            direct = directions[cluster_index]
            #add conditions for no hits and just 1 hit
            if len(clusters[cluster_index]) > 1:
                messagefunc("running hit overlapper...", cols, debugfile)
                #this is a hitlist {hit index: [hit val]}
                hitdict = clusters[cluster_index]
                messagefunc("direction: "+str(direct)+", number of hits: "+str(len(clusters[cluster_index])), cols, debugfile)
                ovlp = True
                #this part will overlap hits of same query=target areas, as well as remove spurious hits
                #this will be a dict with {hit index as in hitlist: bucket}. -1 bucket for dumped hits
                bucketdict = {}
                nbuck = 0
                while ovlp:
                    ovlp = False
                    if len(hitdict) == 1:
                        break
                    else:
                        for refhitkey in hitdict.keys():
                            #processing a particular refhitkey
                            if ovlp:
                                break
                            else:
                                for hitkey, hitval in hitdict.items():
                                    if ovlp:
                                        break
                                    else:
                                        if hitkey != refhitkey:
                                            tovlp = getOverlap(hitval[0:2],hitdict[refhitkey][0:2])
                                            qovlp = getOverlap(hitval[3:5],hitdict[refhitkey][3:5])
                                            ##process overlaps, only one option will work
                                            if tovlp == 0:
                                                if 0 <= qovlp <= hit_overlap:
                                                    ovlp = False
                                                    hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "add", cols, debugfile)
                                                    
                                                #this is this case for subcontig splitter
                                                elif qovlp > hit_overlap:
                                                    #keep both but assign same bucket
                                                    ovlp = False
                                                    hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "samebuck", cols, debugfile)
                                                    
                                            else:
                                                #keep only one, remove the other by removing from hitdict and bucketdict
                                                if tovlp == qovlp:
                                                    # check frame when translating options are used
                                                    if ac2 == "tdna-aa" or ac2 == "tdna-tdna" or ac2 == "aa-tdna":
                                                        startshift = abs(min(hitval[0], hitval[1])-min(hitdict[refhitkey][0], hitdict[refhitkey][1]))-1
                                                        if startshift > 0 and startshift%3 == 0:
                                                            # in frame, extending range of refhit
                                                            hitdict[refhitkey] = extend_hit(hitdict[refhitkey],hitdict[hitkey], tovlp)
                                                            hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "rmhit", cols, debugfile)
                                                        else:
                                                            # not in frame, removing worst
                                                            comp0 = compare_scores(hitdict[refhitkey][0:2],hitdict[refhitkey][6:9], hitdict[hitkey][0:2],hitdict[hitkey][6:9], metric, metricR)
                                                            if comp0 == 0 or comp0 == 2:
                                                                # print >> debugfile, hitdict[hitkey], "deleted4"
                                                                hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "rmhit", cols, debugfile)
                                                            else:
                                                                # print >> debugfile, hitdict[refhitkey], "deleted5"
                                                                hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "rmref", cols, debugfile)
                                                                
                                                    else:
                                                        # frame does not matter for blastn, blastp, hmmer, extend
                                                        hitdict[refhitkey] = extend_hit(hitdict[refhitkey],hitdict[hitkey], tovlp)
                                                        hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "rmhit", cols, debugfile)
                                                else:
                                                    # print "target ovlp is larger than query, remove the worst, set bool to true"
                                                    comp0 = compare_scores(hitdict[refhitkey][0:2],hitdict[refhitkey][6:9], hitdict[hitkey][0:2],hitdict[hitkey][6:9], metric, metricR)
                                                    if comp0 == 0 or comp0 == 2:
                                                        # print >> debugfile, hitdict[hitkey], "deleted7",hitdict[refhitkey][6:9],hitdict[hitkey][6:9]
                                                        hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "rmhit", cols, debugfile)

                                                    else:
                                                        # print >> debugfile, hitdict[refhitkey], "deleted8"
                                                        hitdict, bucketdict, nbuck = sticher_helper(refhitkey, hitkey, hitdict, bucketdict, nbuck, "rmref", cols, debugfile)
                                                ovlp = True
                                                # print hitdict, bucketdict
                                                break
                        if not ovlp:
                            messagefunc("no more overlaps", cols, debugfile)
                #this part will split hits into different subcontig clusters and
                # print hitdict, bucketdict
                subcontigs = [] #this will contain multiple subcontigs in case of splitting
                while True:
                    sbctg = [] #this will contain hits of a particular subcontig, for future stiching
                    for buck2 in set(bucketdict.values()):
                        tempbuck = {}
                        for hitkey, bucket in bucketdict.items():
                            if bucket == buck2:
                                # assemble tempbuck as a subset of hitdict for a given bucket
                                tempbuck[hitkey] = hitdict[hitkey]
                        if len(tempbuck) > 1:
                            #sort hits, pick best , 6--eval float(row[10]), 7--bit float(row[11]), 8--ident float(row[2])
                            bhit = -1
                            bhitval = []
                            for htnm in tempbuck.keys():
                                if len(sbctg) == 0 or synteny_check(sbctg, tempbuck[htnm], direct, max_gap2, cols, debugfile):
                                    if bhit == -1:
                                        bhit = htnm
                                        bhitval = tempbuck[htnm]
                                    else:
                                        comp1 = compare_scores(bhitval[0:2],bhitval[6:9], tempbuck[htnm][0:2],tempbuck[htnm][6:9], metric, metricR)
                                        #if tempbuck[htnm] is better. otherwise keep best as is
                                        if comp1 == 1:
                                            bhit = htnm
                                            bhitval = tempbuck[htnm]
                                            #potentially need breaking and restarting the loop?
                            if bhit != -1:
                                sbctg.append(tempbuck[bhit]) #add best to subcontig
                                del bucketdict[bhit] #remove the best from further assessments
                                del hitdict[bhit] #same
                        elif len(tempbuck) == 1:
                            bhit = tempbuck.keys()[0]
                            if len(sbctg) == 0 or synteny_check(sbctg, tempbuck[bhit], direct, max_gap2, cols, debugfile):
                                sbctg.append(tempbuck[bhit]) #add best to subcontig
                                del bucketdict[bhit] #remove the best from further assessments
                                del hitdict[bhit] #same
                        #else do nothing
                    if len(sbctg) > 0:
                        subcontigs.append(sbctg)
                    else:
                        break #exit while loop
                #this part will stich hits of subcontigs
                stiched_subcontigs = []
                messagefunc("running hit sticher...", cols, debugfile)
                for sbctg in subcontigs:
                    stiched_subcontig = join_chunks(sbctg, direct)
                    stiched_subcontigs.append(stiched_subcontig)
            
            if len(clusters[cluster_index]) == 1:
                sbctg = clusters[cluster_index].values()
                stiched_subcontig = join_chunks(sbctg, direct)
                stiched_subcontigs.append(stiched_subcontig)

        for subcont_index in range(len(stiched_subcontigs)):
            returnlist[indexer_function(querykey,str(subcont_index))] = stiched_subcontigs[subcont_index]
    return returnlist

#checks syntheny and gap between hits
def synteny_check(inplist, item1, direct1, max_gap3, cols1, debugfile1):
    cond = True
    item1medT = median([item1[0],item1[1]])
    item1medQ = median([item1[3],item1[4]])
    min_deltaT = 0
    min_deltaQ = 0
    for item2 in inplist:
        item2medT = median([item2[0],item2[1]])
        item2medQ = median([item2[3],item2[4]])
        deltaT = item2medT - item1medT
        deltaQ = item2medQ - item1medQ
        if direct1:     
            if deltaQ > 0 and deltaT < 0:
                cond = False
                break
            elif deltaQ < 0 and deltaT > 0:
                cond = False
                break
        else:
            if deltaQ > 0 and deltaT > 0:
                cond = False
                break
            elif deltaQ < 0 and deltaT < 0:
                cond = False
                break
        if max_gap3 > 0:
            if abs(deltaT) > abs(min_deltaT):
                min_deltaT = deltaT
            if abs(deltaQ) > abs(min_deltaQ):
                min_deltaQ = deltaQ
    if not cond:
        messagefunc("synteny check: subcontig split, dQ "+str(deltaQ)+", dT "+str(deltaT), cols1, debugfile1)
    if max_gap3 > 0:
        if abs(min_deltaT) > max_gap3 or abs(min_deltaQ) > max_gap3:
            cond = False
            messagefunc("gap check: subcontig split, min dQ "+str(min_deltaQ)+", mit dT "+str(min_deltaT), cols1, debugfile1)
    return cond


#join chunks in correct order
#all coordinates are returned in forward orientation, direction maintained by [direct]
def join_chunks(sbctg1, direct1):
    # print sbctg1, direct1
    median_query = {}
    start_query = {}
    end_query = {}
    start_target = {}
    end_target = {}
    eval_hit = {}
    bit_hit = {}
    ident_hit = {}
    scores_sbctg1 = []
    for chunk in range(len(sbctg1)):
        median_query[chunk] = median([sbctg1[chunk][3],sbctg1[chunk][4]])
        start_query[chunk] = min(sbctg1[chunk][3],sbctg1[chunk][4])
        end_query[chunk] = max(sbctg1[chunk][3],sbctg1[chunk][4])
        start_target[chunk] = min(sbctg1[chunk][0],sbctg1[chunk][1])
        end_target[chunk] = max(sbctg1[chunk][0],sbctg1[chunk][1])
        eval_hit[chunk] = sbctg1[chunk][6]
        bit_hit[chunk] = sbctg1[chunk][7]
        ident_hit[chunk] = sbctg1[chunk][8]
        if chunk == 0:
            scores_sbctg1 = [eval_hit[chunk], bit_hit[chunk], ident_hit[chunk]]
        else:
            #could use max score as an alternative
            scores_sbctg1 = [min(scores_sbctg1[0], eval_hit[chunk]), max(scores_sbctg1[1], bit_hit[chunk]), max(scores_sbctg1[2], ident_hit[chunk])]
            # scores_sbctg1 = amalgamate_scores(scores_sbctg1, [eval_hit[chunk], bit_hit[chunk], ident_hit[chunk]])
        #ranges
        start_target_subctg = min(start_target.values())
        end_target_subctg = max(end_target.values())
        start_query_subctg = min(start_query.values())
        end_query_subctg = max(end_query.values())
    gapstart = 0
    # stiched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score], gap, [range, score], gap ...]]
    stiched_sbctg1 = []
    stiched_sbctg1.append([direct1]) #add direction
    stiched_sbctg1.append(scores_sbctg1) #add scores
    stiched_sbctg1.append([start_target_subctg, end_target_subctg, start_query_subctg, end_query_subctg]) #add ranges
    #add per hit information
    stiched_hits = []
    for key in sorted(median_query, key=lambda x: median_query[x]):
        if gapstart > 0:
            stiched_hits.append(start_query[key]-gapstart)
        stiched_hits.append([start_target[key], end_target[key], start_query[key], end_query[key], eval_hit[key], bit_hit[key], ident_hit[key]])
        gapstart = end_query[key]
    stiched_sbctg1.append(stiched_hits)
    # print stiched_sbctg1
    return stiched_sbctg1

#function to join contigs in correct order
def join_contigs(inplist):
    median_query = {}
    start_query = {}
    end_query = {}
    target_dict = {}
    for target in inplist:
        median_query[target[0]] = median([min(target[1][2][2],target[1][2][3]),max(target[1][2][2],target[1][2][3])])
        start_query[target[0]] = min(target[1][2][2],target[1][2][3])
        end_query[target[0]] = max(target[1][2][2],target[1][2][3])
        target_dict[target[0]] = target[1]
    gapstart = 0
    returnlist = []
    start_query_superctg = min(start_query.values())
    end_query_superctg = max(end_query.values())
    for key in sorted(median_query, key=lambda x: median_query[x]):
        if gapstart > 0:
            returnlist.append(start_query[key]-gapstart)
        returnlist.append([key,target_dict[key]])
        gapstart = end_query[key]
    return [start_query_superctg, end_query_superctg], returnlist

#function to merge overlapping hits. resulting hit gets highest score
def extend_hit(item1, item2, tovlp):
    outlist = []    
    # len1 = float(max(item1[0], item1[1]) - min(item1[0], item1[1]))
    # len2 = float(max(item2[0], item2[1]) - min(item2[0], item2[1]))
    # if len1 >= len2:
    #     #item1 is better
    #     penalty = (len2 - tovlp) / float(len2)
    #     if penalty <= 0:
    #         scores = item1[6:9]
    #     else:
    #         item_scores = item2[6:9]
    #         item_scores[0] = item_scores[0] / penalty
    #         item_scores[1] = item_scores[1] * penalty
    #         scores = amalgamate_scores(item1[6:9], item_scores)
    # else:
    #     #item2 is better
    #     penalty = (len1 - tovlp) / float(len1)
    #     if penalty <= 0:
    #         scores = item2[6:9]
    #     else:
    #         item_scores = item1[6:9]
    #         item_scores[0] = item_scores[0] / penalty
    #         item_scores[1] = item_scores[1] * penalty
    #         scores = amalgamate_scores(item_scores, item2[6:9])
    scores = (min(item1[6], item2[6]), max(item1[7], item2[7]), max(item1[8], item2[8]))
    #using item1 as benchmark for resulting direction
    direction1 = item1[2]
    direction2 = item1[5]
    #if ref item is True
    if direction1:
        outlist.append(min(item1[0], item2[0],item1[1], item2[1])) #min coord first
        outlist.append(max(item1[0], item2[0],item1[1], item2[1])) #max coord second
    else:
        outlist.append(max(item1[0], item2[0],item1[1], item2[1])) #max coord first
        outlist.append(min(item1[0], item2[0],item1[1], item2[1])) #min coord second
    outlist.append(direction1)
    #if query is true
    if direction2:
        outlist.append(min(item1[3], item1[4], item2[3], item2[4])) #min coord first
        outlist.append(max(item1[3], item1[4], item2[3], item2[4])) #max coord second
    else:
        outlist.append(max(item1[3], item1[4], item2[3], item2[4])) #max coord first
        outlist.append(min(item1[3], item1[4], item2[3], item2[4])) #min coord second
    outlist.append(direction2)
    for score in scores:
        outlist.append(score)
    return outlist

#function to combine scores (only used for contig combination for now)
def amalgamate_scores(item1, item2, metricC):
    neweval = item1[0] * item2[0] #probability product
    #test correction
    if neweval != 0.0:
        (sign, digits, exponent) = Decimal(neweval).as_tuple()
        neweval = 1*(10**int(round((exponent+len(digits))*metricC)))
    newbit = (item1[1] + item2[1])*metricC #sum of bits
    newident = (item1[2] + item2[2]) / 2 #average of identities
    return [neweval, newbit, newident]
    
#function to compare scores of two items, takes ranges and scores
#return 0 if first is better
#return 1 if second is better
#return 2 if they are equal
#metric [0,1,2] or [2,1,0] or "e/b-i"
def compare_scores(item1ranges, item1scores, item2ranges, item2scores, metric, metricR):
    len1 = float(max(item1ranges[0], item1ranges[1]) - min(item1ranges[0], item1ranges[1]))
    len2 = float(max(item2ranges[0], item2ranges[1]) - min(item2ranges[0], item2ranges[1]))
    if metricR:
        eval1 = item1scores[0] / len1
        bit1 = item1scores[1] / len1
        ident1 =  item1scores[2] / len1
        eval2 = item2scores[0] / len2
        bit2 = item2scores[1] / len2
        ident2 =  item2scores[2] / len2
    else:
        eval1 = item1scores[0]
        bit1 = item1scores[1]
        ident1 =  item1scores[2]
        eval2 = item2scores[0]
        bit2 = item2scores[1]
        ident2 =  item2scores[2]
    
    metricL = [[eval1,eval2],[bit1,bit2],[ident1,ident2]]
    if metric == "e/b-i":
        if metricL[0][0] < metricL[0][1] and metricL[1][0] > metricL[1][1]:
            return 0
        elif metricL[0][0] > metricL[0][1] and metricL[1][0] < metricL[1][1]:
            return 1
        else:
            if metricL[2][0] > metricL[2][1]:
                return 0
            elif metricL[2][0] < metricL[2][1]:
                return 1
            else:
                return 2
    else:
        #first (0)
        #if 0 higher than 1
        if metricL[metric[0]][0] > metricL[metric[0]][1]:
            #if eval, return 1 as best
            if metric[0] == 0:
                return 1
            #else return 0 as best
            else:
                return 0
        elif metricL[metric[0]][0] < metricL[metric[0]][1]:
            #if eval, return 0 as best
            if metric[0] == 0:
                return 0
            #else return 1 as best
            else:
                return 1
        #first order equal, try second order
        else:
            #if 0 higher than 1
            if metricL[metric[1]][0] > metricL[metric[1]][1]:
                #if eval, return 1 as best
                if metric[1] == 0:
                    return 1
                #else return 0 as best
                else:
                    return 0
            elif metricL[metric[1]][0] < metricL[metric[1]][1]:
                #if eval, return 0 as best
                if metric[1] == 0:
                    return 0
                #else return 1 as best
                else:
                    return 1
            #second order equal, try third order
            else:
                #if 0 higher than 1
                if metricL[metric[2]][0] > metricL[metric[2]][1]:
                    #if eval, return 1 as best
                    if metric[2] == 0:
                        return 1
                    #else return 0 as best
                    else:
                        return 0
                elif metricL[metric[2]][0] < metricL[metric[2]][1]:
                    #if eval, return 0 as best
                    if metric[2] == 0:
                        return 0
                    #else return 1 as best
                    else:
                        return 1
                #third order is equal, return 2
                else:
                    return 2
    
#function to check reciprocity based on whole ranges
#!!! implement direction check for translated searches
def range_reciprocator(targetkey1, inpdict, metric, metricR, recip_overlap, cols, debugfile):
    returnlist = {} #all queries for the target go here
    querylist = inpdict.keys() #list of all queries
    badkeys = set()
    for querykey, queryval in inpdict.items(): #queries for a given target
        if querykey not in badkeys:
            cond = True
            ref_query_range = inpdict[querykey][2]
            ref_query_scores = inpdict[querykey][1]
            for key in querylist: #running loop over other queries
                if key not in badkeys and key != querykey: #all other queries
                    current_query_range = inpdict[key][2]
                    current_query_scores = inpdict[key][1]
                    if getOverlap([current_query_range[0],current_query_range[1]],[ref_query_range[0],ref_query_range[1]]) > recip_overlap:
                        # print >> debugfile, "comparing", key, current_query_range, current_query_scores, "with ref", querykey, ref_query_range, ref_query_scores
                        comp1 = compare_scores(ref_query_range, ref_query_scores, current_query_range, current_query_scores, metric, metricR)
                        if comp1 == 1:
                            cond = False
                            # print >> debugfile, ref_query_range, ref_query_scores, current_query_range, current_query_scores
                            badkeys.add(querykey)
                            break
                        elif comp1 == 0:
                            badkeys.add(key)
                        elif comp1 == 2:
                            wrn = "warning, target "+targetkey1+" at query "+querykey+" has equal hits to query "+key+", saved for both!"
                            warninglist.append(wrn)
                            messagefunc(wrn, cols, debugfile)
                            messagefunc("match to current query: ref_query_scores[0] "+str(ref_query_scores[0])+", bitmax "+str(ref_query_scores[1]), cols, debugfile)
                            messagefunc("match to "+key+": "+",".join(map(str, current_query_scores)), cols, debugfile)
            if cond: #only return good items
                returnlist[querykey] = queryval
    messagefunc(str(len(returnlist))+" queries survived, "+str(len(badkeys))+" queries were filtered out by range reciprocity check", cols, debugfile)
    if len(badkeys) > 0:
        print >> debugfile, badkeys
    return returnlist
    
#second main function
def query_processor(inpdict, rec_dict, target_ref, metric, metricR, metricC, contignum, contig_overlap, interstich, hit_overlap, ac1, recip_overlap, ref_hs, rmrecnf, cols, debugfile):
    returndict = {}
    # 1 reformat target table as query table
    # 2 run reference reciprocation: for each query each target must match back to query transcript from reference
    messagefunc(dashb, cols, debugfile)
    query_dict = reformat_dict(inpdict, rec_dict, target_ref, metric, metricR, hit_overlap, ac1, recip_overlap, ref_hs, rmrecnf, cols, debugfile) #do steps 1 and 2
    # 3 for each query rank targets by scores
    # 4 run sticher, i.e. fill in size of query with contigs in best order, set aside, continue
    messagefunc(dashb, cols, debugfile)
    for querykey, queryval in query_dict.items():
        messagefunc(dash, cols, debugfile)
        messagefunc("query processor on query "+querykey, cols, debugfile)
        if len(queryval) > 1:
            messagefunc("rank targets", cols, debugfile)
            targetlist = rank_targets(queryval, metric, metricR) # step 3
            messagefunc("running contig sticher...", cols, debugfile)
            stiched_targets = contig_sticher(targetlist, metric, metricR, metricC, contig_overlap, interstich, cols, debugfile) # step 4
        else:
            messagefunc("only one target, no ranking and stiching", cols, debugfile)
            reformatted_queryval = [queryval.keys()[0], queryval.values()[0]]
            stiched_targets = [[0, [[True], reformatted_queryval[1][1], reformatted_queryval[1][2],[reformatted_queryval]]]]
        # 5 get top [contignum] contigs from step 4, output
        for tgt in stiched_targets:
            if len(tgt[1][3]) == 1:
                ctgn = "1"
            else:
                ctgn = str(int(len(tgt[1][3]) / 2 )+int(len(tgt[1][3]) % 2 ))
            messagefunc("supercontig "+str(tgt[0])+", number of contigs: "+ctgn+", score "+" ".join([str(x) for x in tgt[1][1]]), cols, debugfile)
        if len(stiched_targets) > contignum and contignum > 0:
            stiched_targets = stiched_targets[:contignum]
        messagefunc("subset by max number of contigs, "+str(len(stiched_targets))+" supercontig passed", cols, debugfile)
        returndict[querykey] = stiched_targets
    return returndict

#function to reformat target based dict into query based, running reference reciprocity check along the way
# returnlist[querykey+"_"+str(subcont_index)] = stiched_subcontigs[subcont_index]
# stiched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score], gap, [range, score], gap ...]]
def reformat_dict(inpdict, rec_dict, target_ref, metric, metricR, hit_overlap, ac1, recip_overlap, ref_hs, rmrecnf, cols, debugfile):
    messagefunc("reformatting dictionaries...", cols, debugfile)
    returnlist = {}
    for targetkey, targetval in inpdict.items():
        for querykey, queryval in targetval.items():
            messagefunc(dash, cols, debugfile)
            messagefunc("processing target "+targetkey+" and query "+querykey, cols, debugfile)
            querykeynew, queryindex = indexer_function(querykey, None)
            targetkeynew = indexer_function(targetkey, queryindex)
            # reciprocal search section
            if rec_dict != None:
                if reference_reciprocator(querykeynew, queryval, rec_dict, target_ref, metric, metricR, hit_overlap, ac1, recip_overlap, targetkey, ref_hs, rmrecnf, cols, debugfile):
                    if querykeynew in returnlist:
                        returnlist[querykeynew][targetkeynew] = queryval
                    else:
                        returnlist[querykeynew] = {targetkeynew: queryval}
            else:
                if querykeynew in returnlist:
                    returnlist[querykeynew][targetkeynew] = queryval
                else:
                    returnlist[querykeynew] = {targetkeynew: queryval}
    return returnlist

#function to run hit sticher on reciprocal and reference tables
def process_aux_tables(inpdict, metric, metricR, hit_overlap, ac2, max_gap1, cols, debugfile):
    returndict = {}
    for querykey, queryval in inpdict.items():
        stiched_hit_dict = hit_sticher(queryval, metric, metricR, hit_overlap, ac2, max_gap1, cols, debugfile)
        returndict[querykey] = stiched_hit_dict
    return returndict

#function to run reference based reciprocity check
#!!! implement direction for translated search
# returnlist[querykey+"_"+str(subcont_index)] = stiched_subcontigs[subcont_index]
# stiched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score], gap, [range, score], gap ...]]
def reference_reciprocator(query, queryval, rec_dict, target_ref, metric, metricR, hit_overlap, ac1, recip_overlap, targetkey1, ref_hs1, rmrecnf1, cols, debugfile):
    messagefunc("running reference based reciprocity check...", cols, debugfile)
    returnlist = {}
    cond = True
    #for a Q in the T, here are the regions
    sample_q_region = queryval[2][2:4]
    sample_t_region = queryval[2][0:2]
    # check out refence matches:
    if query in target_ref: #checking best contig for Q in reference
        messagefunc("reference check", cols, debugfile)
        best_ref_name = []
        best_ref_val = None
        #targets and values for the Q in the ref
        # print target_ref[query].items()
        for target_ref_name, target_ref_val in target_ref[query].items():
            #this will have amalgamated scores and contig ranges for this target - compare with others
            # target_ref_val_stiched = hit_sticher({target_ref_name: target_ref_val}, metric, metricR, hit_overlap, ac1, cols, debugfile)[target_ref_name+"@$0"] #stiched subcontigs per query
            #these are one or more stiched hits
            # only check same Q region in the reference
            target_ref_base_name = indexer_function(target_ref_name, None)[0]
            reference_q_region = target_ref_val[2][2:4]
            reference_t_region = target_ref_val[2][0:2]
            if getOverlap(sample_q_region,reference_q_region) > recip_overlap:
                # the first overlapping region to consider
                if len(best_ref_name) == 0 and best_ref_val == None:
                    best_ref_name = [target_ref_base_name]
                    best_ref_val = target_ref_val
                # otherwise compare with already stored best - as with range reciprocator
                else:
                    comp1 = compare_scores(best_ref_val[2], best_ref_val[1], target_ref_val[2], target_ref_val[1], metric, metricR)
                    if comp1 == 1:
                        best_ref_name = [target_ref_base_name]
                        best_ref_val = target_ref_val
                    elif comp1 == 2:
                        best_ref_name.append(target_ref_base_name)
                        best_ref_val = target_ref_val
        if best_ref_val == None:
            msg = "no same region matches to query "+query+" in ref [should not happen, possibly queries for forward and reciprocal search differ]"
            messagefunc(msg, cols, debugfile, False)
            warninglist.append(msg)
        else:
            msg = "best ref: "+",".join(best_ref_name)+"; "+" ".join([str(x) for x in best_ref_val[2]])+" "+" ".join([str(x) for x in best_ref_val[1]])
            messagefunc(msg, cols, debugfile)
            # check out sample to reference matches
            if targetkey1 in rec_dict: #checking for best contig for Q in current sample
                messagefunc("sample check", cols, debugfile)
                best_rec_name = None
                best_rec_val = None
                for rec_target, rec_hits in rec_dict[targetkey1].items():
                    if ref_hs1:
                        rec_target_base = indexer_function(rec_target, None)[0]
                        # here we checking if it's the same region on query in forward table and reverse table, and then check best match
                        reciprocal_q_region = rec_hits[2][2:4]
                        reciprocal_t_region = rec_hits[2][0:2]
                        if getOverlap(sample_t_region, reciprocal_q_region) > recip_overlap:
                            # the first overlapping region to consider
                            if best_rec_name == None and best_rec_val == None:
                                best_rec_name = rec_target_base
                                best_rec_val = rec_hits
                            else:
                                comp1 = compare_scores(best_rec_val[2], best_rec_val[1], rec_hits[2], rec_hits[1], metric, metricR)
                                if comp1 == 1:
                                    best_rec_name = rec_target_base
                                    best_rec_val = rec_hits
                            # situation when two are equal is not considered
                    else:
                        rec_target_base = rec_target
                        for hit_num, hit_val in rec_hits.items():
                            reciprocal_q_region = hit_val[3:5]
                            reciprocal_t_region = hit_val[0:2]
                            if getOverlap(sample_t_region, reciprocal_q_region) > recip_overlap:
                                # the first overlapping region to consider
                                if best_rec_name == None and best_rec_val == None:
                                    best_rec_name = rec_target_base
                                    best_rec_val = hit_val
                                else:
                                    comp1 = compare_scores(best_rec_val[0:2], best_rec_val[6:9], hit_val[0:2], hit_val[6:9], metric, metricR)
                                    if comp1 == 1:
                                        best_rec_name = rec_target_base
                                        best_rec_val = hit_val
                                # situation when two are equal is not considered
                    
                if best_rec_name == None:
                    msg = "no same region matches to target "+targetkey1+" in reciprocal table [strange, possibly queries for forward and reciprocal search differ]"
                    messagefunc(msg, cols, debugfile)
                    warninglist.append(msg)
                    # cond = False
                else:
                    msg = "best target: "+best_rec_name+", "+" ".join([str(x) for x in best_rec_val])
                    messagefunc(msg, cols, debugfile)
                    if best_rec_name not in best_ref_name:
                        cond = False
                        messagefunc("reciprocator: target "+targetkey1+" removed from query "+query+": reciprocal condition violated", cols, debugfile)
                msg = "coordinates in forward search: "+" ".join([str(x) for x in queryval[2]])+", scores "+" ".join([str(x) for x in queryval[1]])
                messagefunc(msg, cols, debugfile)
            else:
                msg = "no matches to target "+targetkey1+" in reciprocal table [strange, possibly queries for forward and reciprocal search differ]"
                messagefunc(msg, cols, debugfile)
                if rmrecnf:
                    cond = False
                warninglist.append(msg)
    else:
        msg = "no matches to query "+query+" in ref [should not happen, possibly queries for forward and reciprocal search differ]"
        messagefunc(msg, cols, debugfile)
        warninglist.append(msg)
        if rmrecnf:
            cond = False
    print >> debugfile, cond
    return cond

#function to rank targets based on their scores
# {target:[[direction],[scores],[ranges],[[hit range and score], gap, etc]]}
# {target:[ 0[direction], 1[0eval, 1bit, 2ident], 2[ranges], 3[[hit range and score], gap, etc]]}
def rank_targets(inpdict, metric, metricR):
    returnlist = []
    processed_keys = set()
    while len(processed_keys) < len(inpdict):
        best_key = None
        best_val = None
        for targetkey, targetval in inpdict.items():
            if targetkey not in processed_keys: #ignore already sorted
                if best_key == None:
                    best_key = targetkey
                    best_val = targetval
                else:
                    comp1 = compare_scores(best_val[2], best_val[1], targetval[2], targetval[1], metric, metricR)
                    if comp1 == 1:
                        best_key = targetkey
                        best_val = targetval
        returnlist.append([best_key, best_val])
        processed_keys.add(best_key)
    return returnlist
    # [target, [ 0[direction], 1[0eval, 1bit, 2ident], 2[ranges], 3[[hit range and score], gap, etc]]]

#function to check if contigs are overlapping and call functions to merge contigs, their scores, and rank supercontigs after operation is completed
def contig_sticher(inplist, metric, metricR, metricC, contig_overlap, interstich, cols, debugfile):
    ctg_counter = 0
    returndict = {}
    returnlist = []
    removed_contigs = []
    while len(removed_contigs) < len(inplist):
        current_list = []
        supercontig_scores = []
        for contig1 in inplist:
            if contig1 not in removed_contigs:
                if interstich:
                    if len(current_list) == 0:
                        current_list.append(contig1)
                        removed_contigs.append(contig1) #record processed target
                        supercontig_scores = [contig1[1][1][0], contig1[1][1][1], contig1[1][1][2]]
                    else:
                        cond = True
                        for contig2 in current_list:
                            sticher_overlap = getOverlap(contig2[1][2][2:4],contig1[1][2][2:4])
                            if sticher_overlap > contig_overlap:
                                cond = False
                                break
                        if cond:
                            current_list.append(contig1)
                            removed_contigs.append(contig1)
                            supercontig_scores = amalgamate_scores(supercontig_scores, [contig1[1][1][0], contig1[1][1][1], contig1[1][1][2]], metricC)
                else:
                    messagefunc("overlapping disabled", cols, debugfile)
                    current_list.append(contig1)
                    removed_contigs.append(contig1) #record processed target
                    supercontig_scores = [contig1[1][1][0], contig1[1][1][1], contig1[1][1][2]]
                    break
        if len(current_list) > 1:
            supercontig_Qranges, supercontig_contigs = join_contigs(current_list)
            returndict[ctg_counter] = [[True], supercontig_scores, supercontig_Qranges, supercontig_contigs]
            ctg_counter += 1
        else:
            returndict[ctg_counter] = [[True], supercontig_scores, current_list[0][1][2][2:4], [current_list[0]]]
            ctg_counter += 1
    returnlist = rank_targets(returndict, metric, metricR)
    
    return returnlist

    # stiched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score], gap, [range, score], gap ...]]

#function to subset final table and output it in target based format
def reformat_table(inpdict, cols, debugfile):
    returnlist = {}
    for querykey, queryval in inpdict.items():
        survived_queries.add(querykey)
        # total[ contig[ dir[], score[ ], Qrange[], contigs[ contig [ name, info[ ]]]]]
        #        N         N-0       N-1    N-2        N-3    N-3-M   N-3-M-0
        # print "queryval",queryval
        for elemN in queryval: #elemN = N
            for elemM in elemN[1][3]:
                if type(elemM) is not int:
                    targetname = indexer_function(elemM[0], None)[0]
                    survived_targets.add(targetname)
                    if targetname in returnlist:
                        returnlist[targetname].append(querykey)
                    else:
                        returnlist[targetname] = [querykey]
    return returnlist

#function to output final tables into files
def bltableout(output, bltableout_file, table_type):
    for key, value in sorted(output.items()):
        if table_type == "target":
            print >> bltableout_file, key+","+str(len(value))+","+",".join(value)
        elif table_type == "query":
            print >> bltableout_file, key, value


#function to compute overlap between two ranges supplied as lists with start and end
#returns overlap value
def getOverlap(a, b):
    a0=min(a)-1 #start one pos less to detect 1 position overlap
    a1=max(a)
    b0=min(b)-1 #start one pos less to detect 1 position overlap
    b1=max(b)
    return max(0, min(a1, b1) - max(a0, b0))

#function for writing actual sequence to file
def seqwritefunc(sequence, qname, tname, seqname, outM1, dir1, contignum1, lentarget):
    if outM1 == "target":
        fhandle = open(dir1+"/"+tname, "a")
        finalseq = SeqRecord(sequence)
        finalseq.id = seqname
    elif outM1 == "query":
        fhandle = open(dir1+"/"+qname, "a")
        finalseq = SeqRecord(sequence)
        if contignum1 == 1 or seqname == "none":
            finalseq.id = tname
        else:
            finalseq.id = tname+"|"+seqname
    elif outM1 == "combined":
        fhandle = open(dir1+"/combined.fas", "a")
        finalseq = SeqRecord(sequence)
        if lentarget == 1:
            finalseq.id = seqname
        else:
            finalseq.id = tname+"|"+seqname
    finalseq.name =""
    finalseq.description =""
    SeqIO.write(finalseq, fhandle, "fasta")
    fhandle.close()

#function to compute median
#taken from https://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
def median(lst):
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0


#function to extract sequence from file and modify it
def get_sequence(inplist, seq, extractiontype, fls, trans_out1, ac2, metric, metricR, keep_strand):
    finalseq = Seq("")
    seqlen = len(seq.seq)
    direct = inplist[0][0]
    scores = inplist[1]
    ranges = inplist[2]
    hits = inplist[3]
    if extractiontype == "a":
        for i in hits:
            if type(i) is not int:
                if ac2 == "tdna-aa":
                    start = (i[0]+2)/3-1
                    end = i[1]/3
                else:
                    start = i[0]-1
                    end = i[1]
                if direct:
                    tempseq = seq.seq[start:end]
                else:
                    tempseq = seq.seq[start:end].reverse_complement()
                if trans_out1:
                    finalseq += tempseq.translate()
                else:
                    finalseq += tempseq
            else:
                if trans_out1 or ac2 == "aa-aa" or ac2 == "tdna-aa":
                    tempseq2 = "X"
                else:
                    tempseq2 = "N"
                if i > 0:
                    finalseq += Seq(tempseq2*i)
                else:
                    finalseq += Seq(tempseq2)
    elif extractiontype == "b":
        if ac2 == "tdna-aa":
            start = (ranges[0]+2)/3-1
            end = ranges[1]/3
        else:
            start = ranges[0]-1
            end = ranges[1]
        if start - fls < 0:
            start = 0
        else:
            start = start - fls
        if end + fls > seqlen:
            end = seqlen
        else:
            end = end + fls
        finalseq = seq.seq[start:end]
    elif extractiontype == "n":
        finalseq = seq.seq
    elif extractiontype == "s":
        best_hit = None
        for i in hits:
            if type(i) is not int:
                if best_hit == None:
                    best_hit = i
                else:
                    comp1 = compare_scores(best_hit[0:2], best_hit[4:7], i[0:2], i[4:7], metric, metricR)
                    if comp1 == 1:
                        best_hit = i
        if ac2 == "tdna-aa":
            start = (best_hit[0]+2)/3-1
            end = best_hit[1]/3
        else:
            start = best_hit[0]-1
            end = best_hit[1]
        if start - fls < 0:
            start = 0
        else:
            start = start - fls
        if end + fls > seqlen:
            end = seqlen
        else:
            end = end + fls
        if direct:
            tempseq = seq.seq[start:end]
        else:
            tempseq = seq.seq[start:end].reverse_complement()
        if trans_out1:
            finalseq += tempseq.translate()
        else:
            finalseq += tempseq
    if not direct and extractiontype != "a" and extractiontype != "s" and keep_strand == False:
        finalseq = finalseq.reverse_complement()
    return finalseq

#function to join contig sequences and output final sequence
def dumper(inplist, extractiontype, trans_out2, ac2):
    #need to fix spacer as per alphabet
    finalseq = Seq("")
    if trans_out2 or ac2 == "aa-aa" or ac2 == "tdna-aa":
        gapsymbol = "X"
    else:
        gapsymbol = "N"
    for i in inplist:
        if type(i) is not int:
            finalseq += i[1]
        else:
            if i > 0 and extractiontype != "n":
                finalseq += Seq(gapsymbol*i)
            else:
                finalseq += Seq(gapsymbol)
    return finalseq

#function to control sequence processing
def process_fasta(target_db_name1, inputf1, final_table1, final_target_table1, extractiontype1, flanks1, trans_out1, outM1, output_dir1, contignum1, append_name1, ac1, keep_strand, cols1, debugfile1):
    c1 = len(final_target_table1)
    messagefunc("searching for contigs in: "+target_db_name1+", total number of contigs: "+str(c1), cols1, debugfile1, False)
    if outM1 == "query":
        target_set = {}
        for qkey in final_table1.keys():
            target_set[qkey] = set()
    else:
        target_set = {"main": set()}
    for seq in inputf1: #going over seqs in target file
        if seq.id in final_target_table1: #if the seq in target file
            seq_suffix_c = 1
            for qname in final_target_table1[seq.id]: #checking queries of the target
                if outM1 == "query":
                    check_name = qname
                else:
                    check_name = "main"
                for sprcontig in range(len(final_table1[qname])): #going over supercontig indexes
                    num_contigs = len(final_table1[qname][sprcontig][1][3]) #number of contigs in supercontig
                    for t in range(num_contigs): #looking for target in the query table
                        if type(final_table1[qname][sprcontig][1][3][t]) is not int: #if not gap integer
                            tgt_index_name = final_table1[qname][sprcontig][1][3][t][0]
                            targetname = indexer_function(tgt_index_name, None)[0] #process target name
                            if targetname == seq.id and tgt_index_name not in target_set[check_name]: #found target in the query table and this particular version of target wasnt used
                                #get the sequence
                                messagefunc(str(c1)+" EXTRACTING: contig "+targetname+", query "+qname, cols1, debugfile1)
                                s1 = get_sequence(final_table1[qname][sprcontig][1][3][t][1], seq, extractiontype1, flanks1, trans_out1, ac1, metric, metricR, keep_strand)
                                print >> debugfile1, "- EXTRACTING: final seq", s1[:10], "ranges", final_table1[qname][sprcontig][1][3][t][1]
                                if num_contigs == 1:
                                    #check if same target was used:
                                    use_suffix = False
                                    for ts in target_set[check_name]:
                                        if indexer_function(ts, None)[0] == targetname:
                                            use_suffix = True
                                            seq_suffix_c += 1
                                            break
                                    if extractiontype1 == "n":
                                        if not use_suffix:
                                            seqwritefunc(s1, qname, target_db_name1, targetname, outM1, output_dir1, contignum1, append_name1)
                                        #else do nothing
                                    else:
                                        if use_suffix:
                                            seqwritefunc(s1, qname, target_db_name1, indexer_function(targetname, str(seq_suffix_c)), outM1, output_dir1, contignum1, append_name1)
                                        else:
                                            seqwritefunc(s1, qname, target_db_name1, targetname, outM1, output_dir1, contignum1, append_name1)
                                    target_set[check_name].add(tgt_index_name)
                                else:
                                    #need to disable contig stiching when n is selected
                                    final_table1[qname][sprcontig][1][3][t][1] = s1
                                    messagefunc(str(c1)+" BUCKET: contig "+targetname+", query "+qname, cols1, debugfile1)
                                    target_set[check_name].add(tgt_index_name)
                                    #check the bucket
                                    dump_bucket = True
                                    for elem1 in final_table1[qname][sprcontig][1][3]:
                                        if type(elem1) is not int:
                                            if type(elem1[1]) is list:
                                                dump_bucket = False
                                                break
                                    if dump_bucket:
                                        # bucket filled, dump
                                        s1 = dumper(final_table1[qname][sprcontig][1][3], extractiontype1, trans_out1, ac1)
                                        messagefunc(str(c1)+" EXTRACTING: bucket "+qname+" dumped", cols1, debugfile1)
                                        print >> debugfile1, "- EXTRACTING: final seq", s1[:10]#, "ranges", final_table[qname][sprcontig][1]
                                        seqwritefunc(s1, qname, target_db_name1, "Merged_"+qname+"_supercontig_"+str(sprcontig), outM1, output_dir1, contignum1, append_name1)
                                        # clean up
                                        for subb in final_table1[qname][sprcontig][1][3]:
                                            if type(subb) is not int:
                                                subb[1] = ""
            #clean up after all qs are done
            del final_target_table1[seq.id]
            c1 = c1 - 1
        if len(final_target_table1) == 0:
            messagefunc(str(c1)+" search finished", cols1, debugfile1, False)
            break

#function to encode repeating keys and decode them later on
def indexer_function(base1, index1):
    if index1 != None:
        #encode
        return base1+"@"+index1
    else:
        #decode
        splitlist = base1.split("@")
        if len(splitlist) == 2:
            return splitlist[0], splitlist[1]
        else:
            #if @$ was used inside, rejoin but without last element
            return "@".join(splitlist[:-1]), splitlist[-1]

#function to estimate how many queries and targets from the original tables survived
def estimate_survival(init_queries, init_targets, survived_queries, survived_targets, cols, debugfile):
    msg = "number of queries with OK results: "+str(len(survived_queries))
    messagefunc(msg, cols, debugfile, False)
    msg = "number of queries without OK results: "+str(len(init_queries - survived_queries))
    messagefunc(msg, cols, debugfile, False)
    print >> debugfile, "empty query list:"
    print >> debugfile, " ".join(list(init_queries - survived_queries))
    msg = "number of passed targets: "+str(len(survived_targets))
    messagefunc(msg, cols, debugfile, False)
    msg = "number of filtered out targets: "+str(len(init_targets - survived_targets))
    messagefunc(msg, cols, debugfile, False)
    # print >> debugfile, "filtered out targets:"
    # print >> debugfile, " ".join(list(init_targets - survived_targets))

#---------------------------------------------------------------------#
############################ main script ##############################
debugfile_generic = open("absx.log", "w")

#get terminal window size
if not os.popen('stty size', 'r').read():
    cols = 100
else:
    rows, cols = os.popen('stty size', 'r').read().split()
    cols = int(cols)-10

messagefunc("absx run with option "+filefolder+" selected", cols, debugfile_generic, False)
messagefunc("command line parameters: "+' '.join(sys.argv), cols, debugfile_generic, False)


#multi db option
if filefolder == "M":
    #reading the blastfile
    if bt == "blast":
        blastlist = glob.glob(blastfilearg+"/*.blast")
    else:
        blastlist = glob.glob(blastfilearg+"/*.hmmer")
    if not dry_run:
        targetlist = glob.glob(targetf+"/*.fasta")
    if rec_search != None:
        rec_list = glob.glob(rec_search+"/*.blast")
elif filefolder == "S" or "SM":
    blastlist = [blastfilearg]
    if not dry_run:
        targetlist = [targetf]
    if rec_search != None:
        rec_list = [rec_search]

#dry / not dry run
if dry_run:
    messagefunc("dry run, no target files", cols, debugfile_generic, False)
else:
    messagefunc("list of target fasta files detected (mask *.fasta):", cols, debugfile_generic, False)
    for l in targetlist:
        messagefunc(l, cols, debugfile_generic)
    #make modified dir
    messagefunc("make modified dir...", cols, debugfile_generic, False)
    mkdirfunc(output_dir)
    #copy files
    if not noq:
        messagefunc("copy files...", cols, debugfile_generic, False)
        copyfunc(output_dir, cols, debugfile_generic)

#reciprocal table input
if rec_search != None:
    target_ref = readblastfilefunc(target_ref_file, evalue, bitscore, identity, False, acR, False, cols, debugfile_generic)
    target_ref = process_aux_tables(target_ref, metric, metricR, hit_ovlp, acR, max_gap, cols, debugfile_generic)
else:
    target_ref = None


#parsing blast files
messagefunc("parsing blast files...", cols, debugfile_generic, False)

b1 = 0
if len(blastlist) > 1:
    append_name = True
else:
    append_name = False
for b in blastlist:
    b1 += 1 #counter
    init_queries = set()
    init_targets = set()
    survived_queries = set()
    survived_targets = set()
    messagefunc("target table "+str(b1)+" out of "+str(len(blastlist)), cols, debugfile_generic, False)
    #set up sample debug file
    debugfile = open(b.split("/")[-1]+"_absx.log", "w")
    messagefunc("target log started: "+b, cols, debugfile_generic, False)
    #read alignment table
    if bt == "blast":
        output = readblastfilefunc(b, evalue, bitscore, identity, True, ac, True, cols, debugfile) #output 0 is query, 1 is target
    else:
        output = readhmmerfilefunc(b, evalue, bitscore, bt, ac, cols, debugfile)
    #read reciprocal alignment table
    if rec_search != None:
        if rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast" in rec_list:
            rec_out = readblastfilefunc(rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast", None, None, None, False, acr, False, cols, debugfile)
            if ref_hs:
                rec_out = process_aux_tables(rec_out, metric, metricR, hit_ovlp, acr, cols, debugfile)
        elif len(rec_list) == 1:
            rec_out = readblastfilefunc(rec_list[0], None, None, None, False, acr, False, cols, debugfile)
            if ref_hs:
                rec_out = process_aux_tables(rec_out, metric, metricR, hit_ovlp, acr, cols, debugfile)
        else:
            print "problem with finding the reciprocal search file"
            print rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast", rec_list
            sys.exit()
    else:
        rec_out = None

    final_table = {}
    final_target_table = {}
    
    #run target processor
    target_table = target_processor(output, local_rec, metric, metricR, hit_ovlp, recip_ovlp, ac, run_hs, max_gap, cols, debugfile)
    
    #run query processor
    final_table = query_processor(target_table, rec_out, target_ref, metric, metricR, metricC, contignum, ctg_ovlp, interstich, hit_ovlp, ac, recip_ovlp, ref_hs, rmrecnf, cols, debugfile)

    final_target_table = reformat_table(final_table, cols, debugfile)

    estimate_survival(init_queries, init_targets, survived_queries, survived_targets, cols, debugfile)

    qout = open(b.split("/")[-1]+"_qtable.tab", "w")
    tout = open(b.split("/")[-1]+"_ttable.tab", "w")
    bltableout(final_table, qout, "query")
    bltableout(final_target_table,tout, "target")
    qout.close()
    tout.close()

    #####-----------------------------------------------------------------------
    if dry_run:
        messagefunc("dry run, search through target file skipped", cols, debugfile_generic, False)
    else:
        messagefunc("scanning the target fasta file(s)...", cols, debugfile_generic, False)
        
        if filefolder == "SM":
            for seqname in targetlist:
                inputf = SeqIO.parse(seqname, "fasta")
                target_db_name = seqname.split("/")[-1]
                messagefunc("--------scanning the target---------", cols, debugfile)
                process_fasta(target_db_name, inputf, final_table, final_target_table, extractiontype, flanks, trans_out,outM,output_dir, contignum, append_name, ac, keep_strand, cols, debugfile)
        else:
            if filefolder == "M":
                seqname = b[:-6].split("/")[-1]
                target_db_match = [target_db_m1 for target_db_m1 in targetlist if seqname in target_db_m1]
                if len(target_db_match) == 1:
                    inputf = SeqIO.parse(target_db_match[0], "fasta")
                    target_db_name = seqname.split("/")[-1]
                elif len(target_db_match) == 0:
                    msg = "error, the target fasta file "+seqname+" is not found"
                    messagefunc(msg, cols, debugfile, False)
                else:
                    msg = "error, several matches to "+seqname+" are found in the folder"
                    messagefunc(msg, cols, debugfile, False)
            else:
                seqname = targetlist[0]
                inputf = SeqIO.parse(seqname, "fasta")
                target_db_name = seqname.split("/")[-1]
            messagefunc("--------scanning the target---------", cols, debugfile)
            process_fasta(target_db_name, inputf, final_table, final_target_table, extractiontype, flanks, trans_out,outM,output_dir, contignum, append_name, ac, keep_strand, cols, debugfile)

    debugfile.close()

print len(warninglist)
# for w in warninglist:
#     print w

print >> debugfile_generic, "done"
debugfile_generic.close()

print "done"