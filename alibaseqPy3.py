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
import math

parser = argparse.ArgumentParser(description='ALiBaSeq (Alignment-Based Sequence extraction)',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-b', metavar='table', help='alignment table file',
                        dest="blastfilearg",required=True)
required.add_argument('-f', choices=['S','M','SM'], help='file / folder mode',
                        dest="filefolder", required=True)
required.add_argument('-x', choices=['n','s','a','b'],
                        help='extraction type: n (whole contig), '+
                        's (only best hit region), '+
                        'a (extract all hit regions and join them), '+
                        'b (extract region between two outmost hit regions)',
                        dest="extractiontype",required=True)

optional.add_argument('-t', metavar='assembly', help='assembly file',dest="targetf")
optional.add_argument('-q', metavar='query',
                        help='query file(s) to which extracted results are to be appended; '+
                        'if not specified, sequences are extracted into blank files',
                        dest="queryf")
optional.add_argument('-o', metavar='output',
                        help='output folder for modified files with extracted sequences',
                        dest="output", default="alibaseq_out")
optional.add_argument('-s', metavar='logsuffix', help='output log suffix',
                        dest="logsuffix", default="default")
optional.add_argument('--om', choices=['query','target', 'combined'],
                        help='output mode: group in files per bait [query], '+
                        'per sample [target], or combine in a single file [combined]',
                        dest="outM", default="query")
optional.add_argument('-e', metavar='N', help='evalue cutoff',dest="evalue",
                        type=float, default=0.01)
optional.add_argument('-i', metavar='N', help='identity cutoff',dest="identity",
                        type=float, default=0.0)
optional.add_argument('-B', metavar='N', help='bitscore cutoff',dest="bitscore",
                        type=float, default=0.0)
optional.add_argument('-c', metavar='N', help='number of contigs to extract; '+
                        'if set to 0, then extract all contigs; '+
                        'if set to -1, then extract the best and all close matches',
                        dest="contignum", type=int, default=1)
optional.add_argument('--fl', metavar='N', help='flanks on each side in bp',
                        dest="flanks", type=int, default=0)
optional.add_argument('--lr', dest='local_rec', choices=['none','actual','range'],
                        help='local reciprocator setting', default='range')
optional.add_argument('--is', dest='interstitch', action='store_true',
                        help='perform contig stitching', default=False)
# optional.add_argument('--translate', dest='trans_out', action='store_true',
#                         help='translate output (for -x s or -x a)', default=False)
optional.add_argument('--hit-ovlp', metavar='N', help='allowed hit overlap on query, '+
                        '>= 1 in bp, or relative 0 < N < 1',dest="hit_ovlp",
                        type=float, default=0.1)
optional.add_argument('--ctg-ovlp', metavar='N', help='allowed contig overlap on query, '+
                        '>= 1 in bp, or relative 0 < N < 1',dest="ctg_ovlp",
                        type=float, default=0.2)
optional.add_argument('--recip-ovlp', metavar='N', help='contig overlap on query for '+
                        'reciprocator selection, >= 1 in bp, or relative 0 < N < 1',
                        dest="recip_ovlp", type=float, default=10)
optional.add_argument('--bt',
                        choices=['blast','hmmer22', 'hmmer18', 'hmmer15', "lastz", "sam", "bam"],
                        help='alignment table type',dest="bt", default="blast")
optional.add_argument('--btR', choices=['blast','bed'], help='reference alignment table type',
                        dest="btR", default="blast")
optional.add_argument('--ac', choices=['dna-dna', 'tdna-aa', 'aa-tdna', 'aa-aa', 'tdna-tdna'],
                        help='alignment coordinate type',dest="ac", default="dna-dna")
optional.add_argument('--acr', choices=['dna-dna', 'tdna-aa', 'aa-tdna', 'aa-aa', 'tdna-tdna'],
                        help='reciprocal alignment coordinate type',dest="acr", default="dna-dna")
optional.add_argument('--acR', choices=['dna-dna', 'tdna-aa', 'aa-tdna', 'aa-aa', 'tdna-tdna'],
                        help='reference alignment coordinate type',dest="acR", default="dna-dna")
optional.add_argument('-r', metavar='file/folder', help='reciprocal search output file or folder',
                        dest="rec_search")
optional.add_argument('-R', metavar='file', help='bait to reference contig correspondence file',
                        dest="target_ref_file")
optional.add_argument('-m', choices=['e/b-i','e-b-i','b-e-i','i-b-e','i-e-b','b-i-e','e-i-b'],
                        help='order of metrics to use to select best matches '+
                        '(e - evalue, b - bitscore, i - identity)',
                        dest="metric", default="e/b-i")
optional.add_argument('--rescale-metric', dest='metricR', action='store_true',
                        help='divide metric value by length of hit region', default=False)
optional.add_argument('--metric-merge-corr', metavar='N',
                        help='modify combined metric by this value',dest="metricC",
                        type=float, default=1.0)
optional.add_argument('--no-hs', dest='no_hs', action='store_true',
                        help='do not run hit stitcher', default=False)
optional.add_argument('--ref-hs', dest='ref_hs', action='store_true',
                        help='run hit stitcher on reciprocal table (slow)', default=False)
optional.add_argument('--keep-strand', dest='keep_strand', action='store_true',
                        help='keep original contig direction', default=False)
optional.add_argument('--rm-rec-not-found', dest='rmrecnf', action='store_true',
                        help='remove hits without matches in reciprocal search', default=False)
optional.add_argument('--hmmer-global', dest='hmmerg', action='store_true',
                        help='use HMMER contig score instead of domain score', default=False)
optional.add_argument('--amalgamate-hits', dest='amlghitscore', action='store_true',
                        help='combine score for different hits of the same contig',
                        default=False)
optional.add_argument('--max-gap', metavar='N',
                        help='max gap between HSP regions in either query or hit, '+
                        'use 0 for no filtering',dest="max_gap", type=int, default=0)
optional.add_argument('--cname', dest='cname', action='store_true',
                        help='append original contig name to output sequence name',
                        default=False)
optional.add_argument('--both-strands', dest='bstrands', choices=['1','0'],
                        help='allow both strands of the same contig region to be considered',
                        default='1')
optional.add_argument('--srt', metavar='N',
                        help='score ratio threshold, '+
                        'greater which the hits considered be close matches', dest="srt",
                        type=float, default=0.9)
optional.add_argument('--samScore', metavar='metric',
                        help='metric to use for scoring matches', dest="samscore",
                        default="MAPQ")
optional.add_argument('--dd', dest='dd', choices=['all','none','random'],
                        help='in case hit matches several query with exactly equal score, '+
                        'assign such hit to [all queries / '+
                        'none of the queries / at random to only one]', default='none')
optional.add_argument('--log-header', dest='loghead', action='store_true',
                        help='add a header to the table-like log files', default=False)
optional.add_argument('--synteny', dest='syntcheck', choices=['1','0'],
                        help='only stitch hits that are in synteny to query', default='1')
optional.add_argument('-d', metavar='namedelim',
                        help='sample/contigID delimiter, '+
                        'should be a character not encountered in contig or sample names',
                        dest="namedelim", default="|")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()
else:
    args = parser.parse_args()

    blastfilearg = vars(args)["blastfilearg"]
    filefolder = vars(args)["filefolder"]
    bt = vars(args)["bt"]
    hmmerg = vars(args)["hmmerg"]
    if bt == "hmmer18":
        if vars(args)["extractiontype"] != "n":
            print ("only whole contig extraction is supported "+
                    "for hmmer18 tables, option -x ignored")
        if vars(args)["ac"] != "aa-aa":
            print ("hmmer18 tables are only for AA vs AA search, "+
                    "option --ac ignored")
        extractiontype = "n"
        ac = "aa-aa"
    elif bt == "hmmer15":
        if vars(args)["ac"] != "dna-dna":
            print ("hmmer15 tables are only for DNA vs DNA search, "+
                    "option --ac ignored")
        ac = "dna-dna"
        extractiontype = vars(args)["extractiontype"]
    elif bt == "sam" or bt == "bam":
        print ("SAM / BAM input, importing Pysam")
        import pysam
        samscore = vars(args)["samscore"]
        if vars(args)["ac"] != "dna-dna":
            print ("SAM / BAM files are only for DNA vs DNA search, "+
                    "option --ac ignored")
        ac = "dna-dna"
        extractiontype = vars(args)["extractiontype"]
    else:
        ac = vars(args)["ac"]
        extractiontype = vars(args)["extractiontype"]
    if vars(args)["no_hs"] or extractiontype == "n" or extractiontype == "s":
        print ("without hit stitcher, contig stitching is disabled, "+
                "option --is ignored")
        run_hs = False
        interstitch = False
    else:
        run_hs = True
        interstitch = vars(args)["interstitch"]
    if extractiontype == "n" and vars(args)["keep_strand"]:
        keep_strand = True
    else:
        keep_strand = False
    # if vars(args)["trans_out"]:
    #     if extractiontype == "s" or extractiontype == "a":
    #         if ac == "aa-aa" or ac == "tdna-aa":
    #             print ("should not attempt to translate AA sequence, "+
    #                     "option --translate ignored")
    #             trans_out = False
    #         else:
    #             trans_out = True
    #     else:
    #         print ("will only translate exact matches (options -x s or -x a), "+
    #                 "option --translate ignored")
    #         trans_out = False
    # else:
    #     trans_out = False
    trans_out = False

    outM = vars(args)["outM"]
    if vars(args)["targetf"] == None:
        dry_run = True
        noq = True
        print ("option -q ignored")
    else:
        targetf = vars(args)["targetf"]
        dry_run = False
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
    if contignum == 1 and outM == "query":
        cname = vars(args)["cname"]
    else:
        cname = True
    local_rec = vars(args)["local_rec"]
    if vars(args)["rec_search"] == None:
        rec_search = None
    else:
        rec_search = vars(args)["rec_search"]
        if vars(args)["target_ref_file"] == None:
            print ("please specify -R")
            sys.exit()
        else:
            target_ref_file = vars(args)["target_ref_file"]
        acr = vars(args)["acr"]
        acR = vars(args)["acR"]
        btR = vars(args)["btR"]
    rmrecnf = vars(args)["rmrecnf"]

    hit_ovlp = vars(args)["hit_ovlp"]
    if hit_ovlp < 0:
        print ("overlap value must be 0 or positive")
        sys.exit()
    ctg_ovlp = vars(args)["ctg_ovlp"]
    if ctg_ovlp < 0:
        print ("overlap value must be 0 or positive")
        sys.exit()
    recip_ovlp = vars(args)["recip_ovlp"]
    if recip_ovlp < 0:
        print ("overlap value must be 0 or positive")
        sys.exit()
    flanks = vars(args)["flanks"]
    output_dir = vars(args)["output"]
    logsuffix = vars(args)["logsuffix"]
    if bt == "sam" or bt == "bam":
        metric = [1,0,2]
    else:
        if vars(args)["metric"] == "e/b-i":
            metric = "e/b-i"
        else:
            metric = [int(x) for x in vars(args)["metric"].
                        replace("e","0").replace("b","1").
                        replace("i","2").split("-")]
    metricR = vars(args)["metricR"]
    metricC = vars(args)["metricC"]
    ref_hs = vars(args)["ref_hs"]
    max_gap = vars(args)["max_gap"]
    amlghitscore = vars(args)["amlghitscore"]
    if vars(args)["bstrands"] == '1':
        if trans_out or contignum == 1:
            bstrands = True
        else:
            bstrands = False
    else:
        bstrands = False
    srt = vars(args)["srt"]
    dd = vars(args)["dd"]
    loghead = vars(args)["loghead"]
    syntcheck = vars(args)["syntcheck"]
    namedelim = vars(args)["namedelim"]
    if namedelim == "@":
        print ("@ is reserved, please use other delimiter")
        sys.exit()
    elif len(namedelim) != 1:
        print ("delimiter should be a single character")
        sys.exit()

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
        print (msg)
    print (msg, file=f)


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
    ignorecount = 0
    for row in reader:
        if evalue1:
            if float(row[10]) <= evalue1 and \
            float(row[11]) >= bitscore1 and \
            float(row[2]) >= identity1:
                qname = row[0].split("/")[-1]
                if qname[-4::] == ".fas":
                    qname = qname[:-4:]
                tname = row[1]
                if recstats:
                    init_queries.add(qname)
                    init_targets.add(tname)
                if as_target:
                    #populate target table
                    if tname in returndict:
                        if qname in returndict[tname]:
                            returndict[tname][qname][linecounter] = \
                            rowfunc(row, ac3)
                        else:
                            returndict[tname][qname] = \
                            {linecounter: rowfunc(row, ac3)}
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
            else:
                ignorecount += 1
            linecounter += 1
        else:
            qname = row[0].split("/")[-1]
            if qname[-4::] == ".fas":
                qname = qname[:-4:]
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
    return returndict, ignorecount

#function for the blast parser to return a list for dict like this:
#dict[query] = [target_f, target_r, target_b,
#   query_f, query_r, query_b, eval, bitscore, identity]
#query - query name, target_f - target start pos, target_r - target end,
#   target_b - forward or reverse target direction
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
    return [target_f, target_r, target_b,
            query_f, query_r, query_b,
            float(row[10]), float(row[11]),float(row[2])]


#function for parsing hmmer output tables
def readhmmerfilefunc(b, evalue1, bitscore1, bt1, ac1, hmmerg1, cols, debugfile):
    messagefunc("processing "+b, cols, debugfile, False)
    targetdict = {}
    hmmfile = open(b, "rU")
    linecounter = 0
    recordcounter = 0
    ignorecount = 0
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
                qname = line[qname1]
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
                #no identity is reported, so it is set to 0.0 and
                #   does not interfere with scoring system
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
                        if target_r > ctg_length:
                            target_r = ctg_length
                    else:
                        target_f = ctg_length-int(line[target1])*3+frame
                        if target_f < 0:
                            target_f = 0
                        target_r = ctg_length-int(line[target2])*3-2+frame
                if bt1 == "hmmer22" and hmmerg1:
                    outrow = [target_f, target_r, target_b,
                                query_f, query_r, query_b,
                                float(line[6]), float(line[7]), 0.0]
                else:
                    outrow = [target_f, target_r, target_b,
                                query_f, query_r, query_b,
                                float(line[eval1]), float(line[bit1]), 0.0]
                #populate target table
                if tname in targetdict:
                    if qname in targetdict[tname]:
                        targetdict[tname][qname][linecounter] = outrow
                    else:
                        targetdict[tname][qname] = {linecounter: outrow}
                else:
                    targetdict[tname] = {qname: {linecounter: outrow}}
            else:
                ignorecount += 1
            linecounter += 1
    hmmfile.close()
    return targetdict, ignorecount


#function for parsing a bed file (bait regions)
def readbedfilefunc(b, cols, debugfile):
    messagefunc("processing "+b, cols, debugfile, False)
    returndict = {}
    bedfile = open(b, "rU")
    reader = csv.reader(bedfile, delimiter='\t')
    linecounter = 0
    for row in reader:
        qname = row[3] # name field of BED
        tname = row[0] # chrom field of BED
        target_f = int(row[1]) # chromStart field of BED
        target_r = int(row[2]) # chromEnd field of BED
        query_f = 1
        query_r = (target_r - target_f)+1
        query_b = True
        #score field is ignored
        if row[5] == "+":
            target_b = True
        else:
            target_b = False
        returndict[qname] = {tname : {linecounter :
            [target_f, target_r, target_b, query_f, query_r, query_b, 0, 1, 100.0]}}
    print (returndict, file = debugfile)
    bedfile.close()
    return returndict

#function for parsing a lastz output file
def readlastzfilefunc(b, bitscore1, identity1, as_target, recstats, cols, debugfile):
    messagefunc("processing "+b, cols, debugfile, False)
    returndict = {}
    lastzfile = open(b, "rU")
    reader = csv.reader(lastzfile, delimiter='\t')
    linecounter = 0
    ignorecount = 0
    for row in reader:
        if float(row[0]) >= bitscore1 and float(row[14]) >= identity1:
            qname = row[6].split("/")[-1]
            if qname[-4::] == ".fas":
                qname = qname[:-4:]
            tname = row[1]
            if recstats:
                init_queries.add(qname)
                init_targets.add(tname)
            if as_target:
                #populate target table
                if tname in returndict:
                    if qname in returndict[tname]:
                        returndict[tname][qname][linecounter] = rowfunclastz(row)
                    else:
                        returndict[tname][qname] = {linecounter: rowfunclastz(row)}
                else:
                    returndict[tname] = {qname: {linecounter: rowfunclastz(row)}}
            else:
                #populate query table
                if qname in returndict:
                    if tname in returndict[qname]:
                        returndict[qname][tname][linecounter] = rowfunclastz(row)
                    else:
                        returndict[qname][tname] = {linecounter: rowfunclastz(row)}
                else:
                    returndict[qname] = {tname: {linecounter: rowfunclastz(row)}}
        else:
            ignorecount += 1
        linecounter += 1
    lastzfile.close()
    return returndict, ignorecount

#row function for lastz parser
def rowfunclastz(row):
    target_f = int(row[3])+1
    target_r = int(row[4])
    if row[2] == "+":
        target_b = True
    else:
        target_b = False
    query_f = int(row[8])+1
    query_r = int(row[9])
    if row[7] == "+":
        query_b = True
    else:
        query_b = False
    return [target_f, target_r, target_b,
            query_f, query_r, query_b,
            0, float(row[0]), float(row[14])]

#function for parsing SAM / BAM formats
def readsamformat(b, isbinary, bitscore1, samscore1, cols, debugfile):
    if isbinary:
        rmode = 'rb'
    else:
        rmode = 'r'
    messagefunc("processing "+b, cols, debugfile, False)
    returndict = {}
    linecounter = 0
    ignorecount = 0
    samfile = pysam.AlignmentFile(b, rmode)
    for read in samfile.fetch():
        if read.is_unmapped is False:
            qname = read.query_name.split(" ")[0]
            tname = read.reference_name.split(" ")[0]
            query_b = not read.is_reverse
            query_f = read.query_alignment_start+1
            query_r = read.query_alignment_end
            target_f = read.reference_start+1
            target_r = read.reference_end
            cigar_stats = read.get_cigar_stats()
            if cigar_stats[0][7] > 0:
                ident = float(cigar_stats[0][7]) / \
                        (cigar_stats[0][7]+cigar_stats[0][8])*100
            else:
                ident = 100.0
            if samscore1 == "MAPQ":
                quality = read.mapping_quality
            else:
                tags = dict(read.get_tags())
                if samscore1 in tags:
                    quality = tags[samscore1]
                else:
                    messagefunc("SAM / BAM quality metric specified is not found",
                                cols, debugfile, False)
                    sys.exit()
            rowval = [target_f, target_r, True,
                        query_f, query_r, query_b,
                        0, float(quality), ident]
            if tname in returndict:
                if qname in returndict[tname]:
                    returndict[tname][qname][linecounter] = rowval
                else:
                    returndict[tname][qname] = {linecounter: rowval}
            else:
                returndict[tname] = {qname: {linecounter: rowval}}
        else:
            ignorecount += 1
        linecounter += 1
    samfile.close()
    return returndict, ignorecount

#first main function
def target_processor(inpdict, local_rec, metric, metricR, hit_overlap,
                    recip_overlap, ac1, run_hs1, max_gap1, amlghitscore,
                    metricC,bstrands1, filtration_table1, cols, debugfile):
    messagefunc(dashb, cols, debugfile)
    messagefunc("running target processor...", cols, debugfile)
    outdict = {}
    filtration_table1["synteny"] = 0
    filtration_table1["actual"] = 0
    filtration_table1["range"] = 0
    for targetkey, targetval in inpdict.items():
        messagefunc(dash, cols, debugfile)
        messagefunc("target processor on hit "+targetkey, cols, debugfile)
        #run local actual reciprocal check (check each hit of target only matches one query)
        if local_rec == "actual":
            if len(targetval) > 1 or len(list(targetval[x] for x in targetval)[0]) > 1:
                messagefunc("running actual (per HSP) reciprocity check", cols, debugfile)
                targetval = actual_reciprocator(targetval,recip_overlap, bstrands1,
                                                filtration_table1, metric, metricR)
            else:
                messagefunc("only one HSP for this hit, no reciprocity check", cols, debugfile)
        #run hit overlapper and stitcher
        messagefunc("running hit processor", cols, debugfile)
        if run_hs1:
            tgt_proc_out = hit_stitcher(targetval, metric, metricR, hit_overlap, ac1,
                                        max_gap1, amlghitscore, metricC, filtration_table1,
                                        cols, debugfile) #stitched subcontigs per query
        else:
            tgt_proc_out = reformat_hits(targetval, metric, metricR, cols, debugfile)
        #run local range reciprocal check
        if local_rec == "range":
            if len(tgt_proc_out) > 1:
                messagefunc("running range reciprocity check", cols, debugfile)
                #check that each subcontig matches only 1 Q
                tgt_proc_out = range_reciprocator(targetkey, tgt_proc_out, metric,
                                                    metricR, recip_overlap, bstrands1,
                                                    filtration_table1,cols, debugfile)
            else:
                messagefunc("only one query for this hit, no reciprocity check", cols, debugfile)
        outdict[targetkey] = tgt_proc_out
    return outdict


#function to split hits by 'relative' strand (same vs opposite)
#input is dictionary {linecounter: [target_f, target_r, target_b,
#                                   query_f, query_r, query_b,
#                                   float(row[10]), float(row[11]),float(row[2])]}
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
def actual_reciprocator(inpdict, recip_overlap, bstrands2, filtration_table2, metric, metricR):
    returndict = inpdict
    blacklisted = set()
    for refquerykey, refqueryval in inpdict.items():
        #hits with same target and query - check per strand
        clusters = strand_selector(refqueryval)
        cluster_scoring = []
        for cluster in clusters:
            for refhit in list(x for x in refqueryval):
                refhit_t_range = refqueryval[refhit][0:2]
                refhit_score = refqueryval[refhit][6:9]
                for testhit in list(x for x in cluster):
                    if testhit != refhit and (refquerykey, testhit) not in blacklisted:
                        testhit_t_range = refqueryval[testhit][0:2]
                        testhit_score = refqueryval[testhit][6:9]
                        if getOverlap(refhit_t_range, testhit_t_range, recip_overlap) \
                        > recip_overlap:
                            comp1 = compare_scores(refhit_t_range, refhit_score,
                                                    testhit_t_range, testhit_score,
                                                    metric, metricR)
                            if comp1 == 0:
                                blacklisted.add((refquerykey, testhit))
            strand_best_hit = []
            for hitkey, hitval in cluster.items():
                if (refquerykey, hitkey) not in blacklisted:
                    if len(strand_best_hit) == 0:
                        strand_best_hit = hitval
                    else:
                        comp1 = compare_scores(strand_best_hit[0:2], strand_best_hit[6:9],
                                                hitval[0:2], hitval[6:9], metric, metricR)
                        if comp1 == 1:
                            strand_best_hit = hitval
            cluster_scoring.append(strand_best_hit)

        if not bstrands2 and len(cluster_scoring[0]) > 0 and len(cluster_scoring[1]) > 0:
            comp2 = compare_scores(cluster_scoring[0][0:2], cluster_scoring[0][6:9],
                                    cluster_scoring[1][0:2], cluster_scoring[1][6:9],
                                    metric, metricR)
            if comp2 == 1:
                for hitkey in clusters[0]:
                    blacklisted.add((refquerykey, hitkey))
            else:
                for hitkey in clusters[1]:
                    blacklisted.add((refquerykey, hitkey))

        #hits with same target but different queries - check against all strands
        for refhit in list(x for x in refqueryval):
            if (refquerykey, refhit) not in blacklisted:
                for testquerykey, testqueryval in inpdict.items():
                    if testquerykey != refquerykey:
                        for testhit in list(x for x in testqueryval):
                            if testhit != refhit and (testquerykey, testhit) not in blacklisted:
                                testhit_t_range = testqueryval[testhit][0:2]
                                testhit_score = testqueryval[testhit][6:9]
                                if getOverlap(refhit_t_range, testhit_t_range, recip_overlap) > \
                                recip_overlap:
                                    comp1 = compare_scores(refhit_t_range, refhit_score,
                                                            testhit_t_range, testhit_score,
                                                            metric, metricR)
                                    if comp1 == 0:
                                        blacklisted.add((testquerykey, testhit))
    messagefunc("filtered out HSPs: "+str(len(blacklisted))+", listed below:", cols, debugfile)
    for bad_element in blacklisted:
        print (bad_element[0],returndict[bad_element[0]][bad_element[1]], file = debugfile)
        del returndict[bad_element[0]][bad_element[1]]
    filtration_table2["actual"] += len(blacklisted)
    return returndict

#function to reformat the hit tables in case hit stitcher is not run
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
            hits1 = list(clusters[cluster_index][x] for x in clusters[cluster_index])
            for hit_index in range(len(hits1)):
                if best_index1 == None:
                    best_index1 = hit_index
                    best_dir1 = [direct]
                    best_range1 = [min(hits1[hit_index][0:2]),max(hits1[hit_index][0:2]),
                                    min(hits1[hit_index][3:5]),max(hits1[hit_index][3:5])]
                    best_score1 = hits1[hit_index][6:9]
                    best_data1 = best_range1+best_score1
                else:
                    hit_scores = hits1[hit_index][6:9]
                    hit_ranges = [min(hits1[hit_index][0:2]),max(hits1[hit_index][0:2]),
                                    min(hits1[hit_index][3:5]), max(hits1[hit_index][3:5])]
                    hit_data = hit_ranges+hit_scores
                    comp1 = compare_scores(best_range1, best_score1, hit_ranges,
                                            hit_scores, metric, metricR)
                    if comp1 == 1:
                        best_index1 = hit_index
                        best_dir1 = [direct]
                        best_range1 = [min(hits1[hit_index][0:2]),max(hits1[hit_index][0:2]),
                                        min(hits1[hit_index][3:5]), max(hits1[hit_index][3:5])]
                        best_score1 = hits1[hit_index][6:9]
                        best_data1 = best_range1+best_score1
        returnlist[indexer_function(querykey,str(0))] = [best_dir1, best_score1,
                                                        best_range1, best_data1]
        print (querykey, ", selected best hit:", returnlist[indexer_function(querykey,str(0))],
                file = debugfile)
    return returnlist

#function to split hits by 'absolute' strand
def trans_selector(inpdict):
    #split into things per direction
    cluster1 = {}
    cluster2 = {}
    for hitkey, hitval in inpdict.items():
        if hitval[2] == True:
            cluster1[hitkey] = hitval
        else:
            cluster2[hitkey] = hitval
    return [cluster1, cluster2]

#function to process hit of the same target
#   (remove redundant, order and stitch hits, make subcontigs)
def hit_stitcher(inpdict, metric, metricR, hit_overlap, ac2,
                max_gap2, amlghitscore, metricC, filtration_table2,
                cols, debugfile):
    returnlist = {} #all queries for the target go here
    for querykey, queryval in inpdict.items():
        clusters = strand_selector(queryval)
        directions = [True, False]
        stitched_subcontigs = [] #stitched subcontigs for a query go here
        for cluster_index in range(2):
            direct = directions[cluster_index]
            #add condition for no hits?
            if len(clusters[cluster_index]) > 1:
                messagefunc("running hit overlapper...", cols, debugfile)
                #this is a hitlist {hit index: [hit val]}
                hitdict = clusters[cluster_index]
                messagefunc("direction: "+str(direct)+", number of HSPs: "+
                            str(len(clusters[cluster_index])), cols, debugfile)
                subcontigs = []
                if ac2 == "tdna-aa" or ac2 == "tdna-tdna" or ac2 == "aa-tdna":
                    hitdicts = trans_selector(hitdict)
                else:
                    hitdicts = [hitdict]
                for trans_dict in hitdicts:
                    startpoints = {}
                    endpoints = {}
                    for hitkey, hitval in trans_dict.items():
                        startpoints[hitkey] = min(hitval[3], hitval[4])
                        endpoints[hitkey] = max(hitval[3], hitval[4])
                    sorted_startpoints = sorted(startpoints, key=lambda x: startpoints[x])
                    sorted_endpoints = sorted(endpoints, key=lambda x: endpoints[x])
                    currently_processing = []
                    ovlp_processed = []
                    cur_number = 0
                    for i in range(len(sorted_startpoints)):
                        if i == 0:
                            currently_processing.append(trans_dict[sorted_startpoints[i]])
                            ovlp_processed.append([])
                        else:
                            maxtovlp = []
                            maxqovlp = []
                            maxqovlpD = []
                            combovlp = []
                            for cp in range(len(currently_processing)):
                                tovlp = getOverlap(trans_dict[sorted_startpoints[i]][0:2],
                                                    currently_processing[cp][0:2], 1)
                                qovlp = getOverlap(trans_dict[sorted_startpoints[i]][3:5],
                                                    currently_processing[cp][3:5], 1)
                                qovlpD = getOverlap(trans_dict[sorted_startpoints[i]][3:5],
                                                    currently_processing[cp][3:5], hit_overlap)
                                maxtovlp.append(tovlp)
                                maxqovlp.append(qovlp)
                                maxqovlpD.append(qovlpD)
                                combovlp.append(tovlp+qovlp)
                            best_ind = combovlp.index(max(combovlp))
                            #if closest hit overlaps on query
                            if maxqovlpD[best_ind] > hit_overlap:
                                #if also overlaps on target - extend proper layer
                                if maxtovlp[best_ind] > 0:
                                    if ac2 == "tdna-aa" or ac2 == "tdna-tdna" or ac2 == "aa-tdna":
                                        if maxqovlp[best_ind] % 3 == 0:
                                            currently_processing[best_ind] = \
                                            extend_hit(currently_processing[best_ind],
                                                        trans_dict[sorted_startpoints[i]])
                                        else:
                                            comp0 = compare_scores(
                                                    currently_processing[best_ind][0:2],
                                                    currently_processing[best_ind][6:9],
                                                    trans_dict[sorted_startpoints[i]][0:2],
                                                    trans_dict[sorted_startpoints[i]][6:9],
                                                    metric, metricR)
                                            #if current is better, do nothing,
                                            #   else replace current with i
                                            if comp0 == 1:
                                                currently_processing[best_ind] = \
                                                trans_dict[sorted_startpoints[i]]
                                    else:
                                        currently_processing[best_ind] = \
                                        extend_hit(currently_processing[best_ind],
                                                    trans_dict[sorted_startpoints[i]])
                                #else overlaps on query but not target -
                                #   add separate layer, keep current in the current layer
                                else:
                                    currently_processing.append(trans_dict[sorted_startpoints[i]])
                                    ovlp_processed.append([])
                            #if non overlapping on query - add current to the layer,
                            #   make it the new current
                            else:
                                #first check synteny and gap
                                if synteny_check(currently_processing[best_ind],
                                                    trans_dict[sorted_startpoints[i]],
                                                    direct, max_gap2, cols, debugfile):
                                    ovlp_processed[best_ind].append(currently_processing[best_ind])
                                    currently_processing[best_ind] = \
                                    trans_dict[sorted_startpoints[i]]
                                #otherwise, add as separate layer
                                else:
                                    currently_processing.append(trans_dict[sorted_startpoints[i]])
                                    ovlp_processed.append([])
                                    filtration_table2["synteny"] += 1
                    #finalize current
                    for cp in range(len(currently_processing)):
                        ovlp_processed[cp].append(currently_processing[cp])
                    for sbctg in ovlp_processed:
                        subcontigs.append(sbctg)

                messagefunc(str(len(subcontigs))+" pseudocontigs survived", cols, debugfile)
                #this part will stitch hits of subcontigs
                messagefunc("running hit stitcher...", cols, debugfile)
                for sbctg in subcontigs:
                    stitched_subcontig = join_chunks(sbctg, direct, amlghitscore, metricC)
                    stitched_subcontigs.append(stitched_subcontig)
            #if only one hit
            if len(clusters[cluster_index]) == 1:
                sbctg = list(clusters[cluster_index][x] for x in clusters[cluster_index])
                stitched_subcontig = join_chunks(sbctg, direct, amlghitscore, metricC)
                stitched_subcontigs.append(stitched_subcontig)

        for subcont_index in range(len(stitched_subcontigs)):
            returnlist[indexer_function(querykey,str(subcont_index))] = \
            stitched_subcontigs[subcont_index]
    return returnlist

#checks syntheny and gap between hits
def synteny_check(item1, item2, direct1, max_gap3, cols1, debugfile1):
    cond = True
    item1medT = median([item1[0],item1[1]])
    item1medQ = median([item1[3],item1[4]])
    item2medT = median([item2[0],item2[1]])
    item2medQ = median([item2[3],item2[4]])
    deltaT = item2medT - item1medT
    deltaQ = item2medQ - item1medQ
    if syntcheck:
        if direct1:
            if deltaQ > 0 and deltaT < 0:
                cond = False
            elif deltaQ < 0 and deltaT > 0:
                cond = False
        else:
            if deltaQ > 0 and deltaT > 0:
                cond = False
            elif deltaQ < 0 and deltaT < 0:
                cond = False
    if not cond:
        messagefunc("synteny check: direction "+str(direct1)+", dQ "+
                    str(deltaQ)+", dT "+str(deltaT), cols1, debugfile1)
    if max_gap3 > 0:
        if abs(deltaT) > max_gap3 or abs(deltaQ) > max_gap3:
            cond = False
            messagefunc("gap check: direction "+str(direct1)+", min dQ "+
                        str(deltaQ)+", min dT "+str(deltaT), cols1, debugfile1)
    return cond


#join chunks in correct order
#all coordinates are returned in forward orientation, direction maintained by [direct]
def join_chunks(sbctg1, direct1, amlghitscore1, metricC1):
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
            if amlghitscore1:
                scores_sbctg1 = amalgamate_scores(scores_sbctg1, [eval_hit[chunk],
                                                    bit_hit[chunk], ident_hit[chunk]],
                                                    metricC1)
            else:
                scores_sbctg1 = [min(scores_sbctg1[0], eval_hit[chunk]),
                                max(scores_sbctg1[1], bit_hit[chunk]),
                                max(scores_sbctg1[2], ident_hit[chunk])]
        #ranges
        start_target_subctg = min(list(start_target[x] for x in start_target))
        end_target_subctg = max(list(end_target[x] for x in end_target))
        start_query_subctg = min(list(start_query[x] for x in start_query))
        end_query_subctg = max(list(end_query[x] for x in end_query))
    gapstart = 0
    # stitched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score],
    #                                                           gap,
    #                                                          [range, score],
    #                                                           gap ...]]
    stitched_sbctg1 = []
    stitched_sbctg1.append([direct1]) #add direction
    stitched_sbctg1.append(scores_sbctg1) #add scores
    stitched_sbctg1.append([start_target_subctg, end_target_subctg,
                            start_query_subctg, end_query_subctg]) #add ranges
    #add per hit information
    stitched_hits = []
    for key in sorted(median_query, key=lambda x: median_query[x]):
        if gapstart > 0:
            stitched_hits.append(start_query[key]-1-gapstart)
        stitched_hits.append([start_target[key], end_target[key],
                            start_query[key], end_query[key],
                            eval_hit[key], bit_hit[key], ident_hit[key]])
        gapstart = end_query[key]
    stitched_sbctg1.append(stitched_hits)
    return stitched_sbctg1

#function to join contigs in correct order
def join_contigs(inplist):
    median_query = {}
    start_query = {}
    end_query = {}
    target_dict = {}
    for target in inplist:
        median_query[target[0]] = median([min(target[1][2][2],target[1][2][3]),
                                        max(target[1][2][2],target[1][2][3])])
        start_query[target[0]] = min(target[1][2][2],target[1][2][3])
        end_query[target[0]] = max(target[1][2][2],target[1][2][3])
        target_dict[target[0]] = target[1]
    gapstart = 0
    returnlist = []
    start_query_superctg = min(list(start_query[x] for x in start_query))
    end_query_superctg = max(list(end_query[x] for x in end_query))
    for key in sorted(median_query, key=lambda x: median_query[x]):
        if gapstart > 0:
            returnlist.append(start_query[key]-1-gapstart)
        returnlist.append([key,target_dict[key]])
        gapstart = end_query[key]
    return [start_query_superctg, end_query_superctg], returnlist

#function to merge overlapping hits. resulting hit gets highest score
def extend_hit(item1, item2):
    outlist = []
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

#function to combine scores (used for hit (only if --amalgamate-hits) and
#   contig (always) stitching)
def amalgamate_scores(item1, item2, metricC):
    neweval = item1[0] * item2[0] #probability product
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

#function to check reciprocity based on stitched match ranges
#   (rather than individual HSP ranges)
def range_reciprocator(targetkey1, inpdict, metric, metricR, recip_overlap,
                        bstrands2, filtration_table2, cols, debugfile):
    returnlist = {} #all queries for the target go here
    querylist = list(x for x in inpdict) #list of all queries
    badkeys = set()
    for querykey, queryval in inpdict.items(): #queries for a given target
        if querykey not in badkeys:
            cond = True
            ref_query_range = inpdict[querykey][2]
            ref_query_scores = inpdict[querykey][1]
            ref_query_direct = inpdict[querykey][0]
            for key in querylist: #running loop over other queries
                if key not in badkeys:
                    #if both-strands to be considered,
                    #   do not compare strands of same contig
                    #if both-strands are off,
                    #   compare anything including diff strands of same contig,
                    #   and discard worse
                    #do not issue a warning of strands
                    if (bstrands2 and key.split("@")[0] != querykey.split("@")[0]) or \
                    (not bstrands2 and key != querykey): #all other queries
                        current_query_range = inpdict[key][2]
                        current_query_scores = inpdict[key][1]
                        current_query_direct = inpdict[key][0]
                        if getOverlap([current_query_range[0],current_query_range[1]],
                                        [ref_query_range[0],ref_query_range[1]],
                                        recip_overlap) > recip_overlap:
                            comp1 = compare_scores(ref_query_range, ref_query_scores,
                                                    current_query_range, current_query_scores,
                                                    metric, metricR)
                            if comp1 == 1:
                                cond = False
                                badkeys.add(querykey)
                                if querykey in returnlist:
                                    del returnlist[querykey]
                                break
                            elif comp1 == 0:
                                badkeys.add(key)
                            elif comp1 == 2:
                                #if different hits
                                if key.split("@")[0] != querykey.split("@")[0]:
                                    wrn = "warning, hit "+targetkey1+" at query "+ \
                                        querykey+" has equal matches to query "+key
                                    warninglist.append(wrn)
                                    messagefunc(wrn, cols, debugfile)
                                    messagefunc("match to current query "+querykey+
                                                ": ref_query_scores[0] "+
                                                str(ref_query_scores[0])+
                                                ", bitmax "+str(ref_query_scores[1]),
                                                cols, debugfile)
                                    messagefunc("match to "+key+": "+
                                                ",".join(map(str, current_query_scores)),
                                                cols, debugfile)
                                    if dd == "none":
                                        cond = False
                                        badkeys.add(key)
                                        badkeys.add(querykey)
                                        if querykey in returnlist:
                                            del returnlist[querykey]
                                        break
                                    elif dd == "random":
                                        badkeys.add(key)
                                else:
                                    #different strands of the same hit
                                    messagefunc("two strands with same score: "+querykey+
                                                ", "+key+", removing the latter...",
                                                cols, debugfile)
                                    badkeys.add(key)

            if cond: #only return good items
                returnlist[querykey] = queryval
    messagefunc(str(len(returnlist))+" queries survived, "+str(len(badkeys))+
                " queries were filtered out by range reciprocity check",
                cols, debugfile)
    if len(badkeys) > 0:
        print (badkeys, file = debugfile)
    filtration_table2["range"] += len(badkeys)
    print ("survived:", returnlist, file = debugfile)
    return returnlist

#second main function
def query_processor(inpdict, rec_dict, target_ref, metric, metricR, metricC,
                    contignum, contig_overlap, interstitch, hit_overlap, ac1,
                    recip_overlap, ref_hs, rmrecnf, srt, filtration_table1,
                    cols, debugfile):
    filtration_table1["contig"] = 0
    filtration_table1["rbh"] = 0
    returndict = {}
    # 1 reformat target table as query table
    # 2 run reference reciprocation:
    #   for each query each target must match back to query contig from reference
    messagefunc(dashb, cols, debugfile)
    query_dict = reformat_dict(inpdict, rec_dict, target_ref, metric, metricR,
                                hit_overlap, ac1, recip_overlap, ref_hs, rmrecnf,
                                filtration_table1, cols, debugfile) #do steps 1 and 2
    # 3 for each query rank targets by scores
    # 4 run stitcher, i.e. fill in size of query with contigs in best order,
    #   set aside, continue
    messagefunc(dashb, cols, debugfile)
    for querykey, queryval in query_dict.items():
        messagefunc(dash, cols, debugfile)
        messagefunc("query processor on query "+querykey, cols, debugfile)
        if len(queryval) > 1:
            if contignum == 0 and interstitch == False:
                messagefunc("no ranking (extract all and no contig stitching)",
                            cols, debugfile)
                targetlist = []
                for targetkey, targetval in queryval.items():
                    targetlist.append([targetkey, targetval])
                stitched_targets = []
                ctg_counter = 0
                for contig1 in targetlist:
                    supercontig_scores = [contig1[1][1][0],
                                            contig1[1][1][1],
                                            contig1[1][1][2]]
                    stitched_targets.append([ctg_counter, [[True], supercontig_scores,
                                            contig1[1][2][2:4], [contig1]]])
                    ctg_counter += 1
            else:
                messagefunc("rank hits", cols, debugfile)
                targetlist = rank_targets(queryval, metric, metricR) # step 3
                if interstitch:
                    messagefunc("running contig stitcher...", cols, debugfile)
                    stitched_targets = contig_stitcher(targetlist, metric, metricR,
                                                        metricC, contig_overlap,
                                                        cols, debugfile) # step 4
                else:
                    messagefunc("stitching disabled", cols, debugfile)
                    stitched_targets = []
                    ctg_counter = 0
                    for contig1 in targetlist:
                        supercontig_scores = [contig1[1][1][0],
                                                contig1[1][1][1],
                                                contig1[1][1][2]]
                        stitched_targets.append([ctg_counter, [[True], supercontig_scores,
                                                contig1[1][2][2:4], [contig1]]])
                        ctg_counter += 1
        else:
            messagefunc("only one hit, no ranking and stitching", cols, debugfile)
            reformatted_queryval = [list(x for x in queryval)[0],
                                    list(queryval[x] for x in queryval)[0]]
            stitched_targets = [[0, [[True], reformatted_queryval[1][1],
                                reformatted_queryval[1][2],[reformatted_queryval]]]]
        # 5 get top [contignum] contigs from step 4, output
        for tgt in stitched_targets:
            if len(tgt[1][3]) == 1:
                ctgn = "1"
            else:
                ctgn = str(int(len(tgt[1][3]) / 2 )+int(len(tgt[1][3]) % 2 ))
            messagefunc("supercontig "+str(tgt[0])+", number of contigs: "+ctgn+
                        ", score "+" ".join([str(x) for x in tgt[1][1]]), cols, debugfile)
        if contignum == -1:
            if len(stitched_targets) > 1:
                filtration_table1["contig"] += len(stitched_targets)
                stitched_targets = \
                stitched_targets[:subset_stiched_targets(stitched_targets, srt)]
                filtration_table1["contig"] -= len(stitched_targets)
                if len(stitched_targets) > 1:
                    msg = "bait "+querykey+" has close suboptimal hits"
                    messagefunc(msg, cols, debugfile)
                    warninglist.append(msg)
                messagefunc("subset by suboptimal scores, "+str(len(stitched_targets))+
                            " supercontig passed", cols, debugfile)
        elif contignum > 0:
            if len(stitched_targets) > contignum:
                dltE, dltB = score_ratio(stitched_targets[contignum-1],
                                            stitched_targets[contignum])
                if dltE > srt or dltB > srt:
                    msg = "bait "+querykey+" has close suboptimal hits "+ \
                            str(dltE)+", "+str(dltB)
                    messagefunc(msg, cols, debugfile)
                    warninglist.append(msg)
                filtration_table1["contig"] += len(stitched_targets)
                stitched_targets = stitched_targets[:contignum]
                filtration_table1["contig"] -= len(stitched_targets)
                messagefunc("subset by max number of contigs, "+str(len(stitched_targets))+
                            " supercontig passed", cols, debugfile)
        returndict[querykey] = stitched_targets
    return returndict

#function to select N suboptimal matches (close to the best match)
def subset_stiched_targets(targets1, srt1):
    subopt = True
    for tgt1 in range(1,len(targets1)):
        ref = targets1[tgt1-1]
        target1 = targets1[tgt1]
        dltE, dltB = score_ratio(ref,target1)
        if dltE <= srt1 and dltB <= srt1:
            subopt = False
            break
    if subopt:
        return tgt1+1
    else:
        return tgt1

#function to compute score differential
def score_ratio(scores1, scores2):
    refE = scores1[1][1][0]
    testE = scores2[1][1][0]
    if refE != 0.0 and testE != 0.0:
        deltaE = max(math.log(refE),math.log(testE))/min(math.log(refE),math.log(testE))
    else:
        deltaE = 0
    refB = scores1[1][1][1]
    testB = scores2[1][1][1]
    deltaB = min(refB, testB)/max(refB, testB)
    return (deltaE, deltaB)

#function to reformat target based dict into query based,
#   running reference reciprocity check along the way
# returnlist[querykey+"_"+str(subcont_index)] = stitched_subcontigs[subcont_index]
# stitched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score],
#                                                           gap,
#                                                           [range, score],
#                                                           gap ...]]
def reformat_dict(inpdict, rec_dict, target_ref, metric, metricR, hit_overlap,
                ac1, recip_overlap, ref_hs, rmrecnf, filtration_table2,
                cols, debugfile):
    messagefunc("reformatting dictionaries...", cols, debugfile)
    returnlist = {}
    for targetkey, targetval in inpdict.items():
        for querykey, queryval in targetval.items():
            messagefunc(dash, cols, debugfile)
            messagefunc("processing hit "+targetkey+" and query "+querykey,
                        cols, debugfile)
            querykeynew, queryindex = indexer_function(querykey, None)
            targetkeynew = indexer_function(targetkey, queryindex)
            # reciprocal search section
            if rec_dict != None:
                if reference_reciprocator(querykeynew, queryval, rec_dict,
                                            target_ref, metric, metricR, hit_overlap,
                                            ac1, recip_overlap, targetkey, ref_hs,
                                            rmrecnf, cols, debugfile):
                    if querykeynew in returnlist:
                        returnlist[querykeynew][targetkeynew] = queryval
                    else:
                        returnlist[querykeynew] = {targetkeynew: queryval}
                else:
                    filtration_table1["rbh"] += 1
            else:
                if querykeynew in returnlist:
                    returnlist[querykeynew][targetkeynew] = queryval
                else:
                    returnlist[querykeynew] = {targetkeynew: queryval}
    return returnlist

#function to run hit stitcher on reciprocal and reference tables
def process_aux_tables(inpdict, metric, metricR, hit_overlap, ac2, max_gap1,
                        amlghitscore, metricC, filtration_table1, cols, debugfile):
    returndict = {}
    for querykey, queryval in inpdict.items():
        stitched_hit_dict = hit_stitcher(queryval, metric, metricR, hit_overlap, ac2,
                                        max_gap1, amlghitscore, metricC, filtration_table1,
                                        cols, debugfile)
        returndict[querykey] = stitched_hit_dict
    return returndict

#function to run reference based reciprocity check
# returnlist[querykey+"_"+str(subcont_index)] = stitched_subcontigs[subcont_index]
# stitched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score],
#                                                           gap,
#                                                           [range, score],
#                                                           gap ...]]
def reference_reciprocator(query, queryval, rec_dict, target_ref, metric,
                            metricR, hit_overlap, ac1, recip_overlap,
                            targetkey1, ref_hs1, rmrecnf1, cols, debugfile):
    messagefunc("running reference based reciprocity check...", cols, debugfile)
    returnlist = {}
    cond = True
    #for a Q in the T, here are the regions
    sample_q_region = queryval[2][2:4]
    sample_t_region = queryval[2][0:2]
    # check out refence matches:
    if query in target_ref: #checking best contig for Q in reference
        messagefunc("reference check for "+query, cols, debugfile)
        best_ref_name = []
        best_ref_val = None
        #targets and values for the Q in the ref
        for target_ref_name, target_ref_val in target_ref[query].items():
            #this will have amalgamated scores and
            #   contig ranges for this target - compare with others
            #stitched subcontigs per query
            # target_ref_val_stitched = hit_stitcher({target_ref_name: target_ref_val},
            #                                         metric, metricR, hit_overlap,
            #                                         ac1, cols, debugfile)
            #                                         [target_ref_name+"@$0"]
            #these are one or more stitched hits
            # only check same Q region in the reference
            target_ref_base_name = indexer_function(target_ref_name, None)[0]
            reference_q_region = target_ref_val[2][2:4]
            reference_t_region = target_ref_val[2][0:2]
            if getOverlap(sample_q_region,reference_q_region, recip_overlap) > recip_overlap:
                # the first overlapping region to consider
                if len(best_ref_name) == 0 and best_ref_val == None:
                    best_ref_name = [target_ref_base_name]
                    best_ref_val = target_ref_val
                # otherwise compare with already stored best - as with range reciprocator
                else:
                    comp1 = compare_scores(best_ref_val[2], best_ref_val[1],
                                            target_ref_val[2], target_ref_val[1],
                                            metric, metricR)
                    if comp1 == 1:
                        best_ref_name = [target_ref_base_name]
                        best_ref_val = target_ref_val
                    elif comp1 == 2:
                        best_ref_name.append(target_ref_base_name)
                        best_ref_val = target_ref_val
        if best_ref_val == None:
            msg = "no same region matches to query "+query+" in ref"
            messagefunc(msg, cols, debugfile, False)
            warninglist.append(msg)
        else:
            msg = "best ref: "+",".join(best_ref_name)+"; "+ \
                " ".join([str(x) for x in best_ref_val[2]])+" "+ \
                " ".join([str(x) for x in best_ref_val[1]])
            messagefunc(msg, cols, debugfile)
            # check out sample to reference matches
            if targetkey1 in rec_dict: #checking for best contig for Q in current sample
                messagefunc("best match check on "+targetkey1, cols, debugfile)
                best_rec_name = None
                best_rec_val = None
                for rec_target, rec_hits in rec_dict[targetkey1].items():
                    if ref_hs1:
                        rec_target_base = indexer_function(rec_target, None)[0]
                        # here we checking if it's the same region on query in
                        #   forward table and reverse table, and then check best match
                        reciprocal_q_region = rec_hits[2][2:4]
                        reciprocal_t_region = rec_hits[2][0:2]
                        if getOverlap(sample_t_region, reciprocal_q_region,
                                        recip_overlap) > recip_overlap:
                            # the first overlapping region to consider
                            if best_rec_name == None and best_rec_val == None:
                                best_rec_name = rec_target_base
                                best_rec_val = rec_hits
                            else:
                                comp1 = compare_scores(best_rec_val[2], best_rec_val[1],
                                                        rec_hits[2], rec_hits[1],
                                                        metric, metricR)
                                if comp1 == 1:
                                    best_rec_name = rec_target_base
                                    best_rec_val = rec_hits
                            # situation when two are equal is not considered
                    else:
                        rec_target_base = rec_target
                        for hit_num, hit_val in rec_hits.items():
                            reciprocal_q_region = hit_val[3:5]
                            reciprocal_t_region = hit_val[0:2]
                            if getOverlap(sample_t_region, reciprocal_q_region,
                                            recip_overlap) > recip_overlap:
                                # the first overlapping region to consider
                                if best_rec_name == None and best_rec_val == None:
                                    best_rec_name = rec_target_base
                                    best_rec_val = hit_val
                                else:
                                    comp1 = compare_scores(best_rec_val[0:2],
                                                            best_rec_val[6:9],
                                                            hit_val[0:2],
                                                            hit_val[6:9],
                                                            metric, metricR)
                                    if comp1 == 1:
                                        best_rec_name = rec_target_base
                                        best_rec_val = hit_val
                                # situation when two are equal is not considered

                if best_rec_name == None:
                    msg = "no same region matches to hit "+targetkey1+" in reciprocal table"
                    messagefunc(msg, cols, debugfile)
                    if rmrecnf:
                        cond = False
                    warninglist.append(msg)
                else:
                    msg = "best hit: "+best_rec_name+", "+ \
                        " ".join([str(x) for x in best_rec_val])
                    messagefunc(msg, cols, debugfile)
                    if best_rec_name not in best_ref_name:
                        cond = False
                        messagefunc("reciprocator: hit "+targetkey1+
                                    " removed from query "+query+
                                    ": reciprocal condition violated",
                                    cols, debugfile)
                msg = "coordinates in forward search: "+ \
                        " ".join([str(x) for x in queryval[2]])+ \
                        ", scores "+" ".join([str(x) for x in queryval[1]])
                messagefunc(msg, cols, debugfile)
            else:
                msg = "no matches to hit "+targetkey1+" in reciprocal table"
                messagefunc(msg, cols, debugfile)
                if rmrecnf:
                    cond = False
                warninglist.append(msg)
    else:
        msg = "no matches to query "+query+" in ref"
        messagefunc(msg, cols, debugfile)
        warninglist.append(msg)
        if rmrecnf:
            cond = False
    return cond

#function to rank targets based on their scores
# {target:[[direction],[scores],[ranges],[[hit range and score], gap, etc]]}
# {target:[ 0[direction], 1[0eval, 1bit, 2ident],
#                           2[ranges], 3[[hit range and score], gap, etc]]}
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
                    comp1 = compare_scores(best_val[2], best_val[1],
                                            targetval[2], targetval[1],
                                            metric, metricR)
                    if comp1 == 1:
                        best_key = targetkey
                        best_val = targetval
        returnlist.append([best_key, best_val])
        processed_keys.add(best_key)
    return returnlist
    # [target, [ 0[direction], 1[0eval, 1bit, 2ident],
    #                           2[ranges], 3[[hit range and score], gap, etc]]]

#function to check if contigs are overlapping and call functions to merge contigs,
#   their scores, and rank supercontigs after operation is completed
def contig_stitcher(inplist, metric, metricR, metricC, contig_overlap, cols, debugfile):
    ctg_counter = 0
    returndict = {}
    returnlist = []
    removed_contigs = []
    while len(removed_contigs) < len(inplist):
        current_list = []
        supercontig_scores = []
        for contig1 in inplist:
            if contig1 not in removed_contigs:
                if len(current_list) == 0:
                    current_list.append(contig1)
                    removed_contigs.append(contig1) #record processed target
                    supercontig_scores = [contig1[1][1][0],
                                            contig1[1][1][1],
                                            contig1[1][1][2]]
                else:
                    cond = True
                    for contig2 in current_list:
                        #check names - do not allow same contig stitching -
                        #   reserved for hit stitcher only now:
                        if contig1[0].split("@")[0] == contig2[0].split("@")[0]:
                            cond = False
                            break
                        stitcher_overlap = getOverlap(contig2[1][2][2:4],
                                                        contig1[1][2][2:4],
                                                        contig_overlap)
                        if stitcher_overlap > contig_overlap:
                            cond = False
                            break
                    if cond:
                        current_list.append(contig1)
                        removed_contigs.append(contig1)
                        supercontig_scores = amalgamate_scores(supercontig_scores,
                                                                [contig1[1][1][0],
                                                                contig1[1][1][1],
                                                                contig1[1][1][2]],
                                                                metricC)
        # print >> debugfile, "current_list",current_list
        if len(current_list) > 1:
            supercontig_Qranges, supercontig_contigs = join_contigs(current_list)
            returndict[ctg_counter] = [[True], supercontig_scores, supercontig_Qranges,
                                        supercontig_contigs]
            ctg_counter += 1
        else:
            returndict[ctg_counter] = [[True], supercontig_scores,
                                        current_list[0][1][2][2:4],
                                        [current_list[0]]]
            ctg_counter += 1
    returnlist = rank_targets(returndict, metric, metricR)
    return returnlist

    # stitched_subcontigs = [ctg_counter, [[direct],[scores],[ranges],
    #                       [hits: [range, score], gap, [range, score], gap ...]]]

#function to subset final table and output it in target based format
def reformat_table(inpdict, cols, debugfile):
    returnlist = {}
    for querykey, queryval in inpdict.items():
        survived_queries.add(querykey)
        # total[ contig[ dir[], score[ ], Qrange[], contigs[ contig [ name, info[ ]]]]]
        #        N         N-0       N-1    N-2        N-3    N-3-M   N-3-M-0
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
    if table_type == "target":
        for key, value in sorted(output.items()):
            print (key+","+str(len(value))+","+",".join(value), file = bltableout_file)
    elif table_type == "query":
        for key, value in sorted(output.items()):
            print (key, value, file = bltableout_file)
    elif table_type == "queryH":
        print("query,supercontig,hit,PCindex,direction,PCstart,PCend,Qstart,Qend,Eval,Bit,Ident",
                file= bltableout_file)
        for key, value in sorted(output.items()):
            for superCTG in value:
                for pseudoCTG in superCTG[1][3]:
                    if type(pseudoCTG) is not int:
                        pseudoCTGfname = pseudoCTG[0].split("@")
                        CTGname = pseudoCTGfname[0]
                        pseudoCTGindex = pseudoCTGfname[1]
                        pseudoCTGdir = str(pseudoCTG[1][0][0])
                        pseudoCTGranges = ",".join(str(x) for x in pseudoCTG[1][2])
                        pseudoCTGscores = ",".join(str(x) for x in pseudoCTG[1][1])
                        print (key+","+str(superCTG[0])+","+CTGname+","+pseudoCTGindex+
                                ","+pseudoCTGdir+","+pseudoCTGranges+","+pseudoCTGscores,
                                file = bltableout_file)


#function to compute overlap between two ranges supplied as lists with start and end
#returns absolute or relative overlap value
def getOverlap(a, b, o):
    a0=min(a)-1 #start one pos less to detect 1 position overlap
    a1=max(a)
    b0=min(b)-1 #start one pos less to detect 1 position overlap
    b1=max(b)
    if 0 < o < 1:
        scaling_len = float(min(a1-a0, b1-b0))
        return (min(a1, b1) - max(a0, b0))/scaling_len
    else:
        if min(a1, b1) - max(a0, b0) < min(abs(a1-a0),abs(b1-b0)):
            return min(a1, b1) - max(a0, b0)
        else:
            return o+1



#function for writing actual sequence to file
def seqwritefunc(sequence, qname, tname, seqname, outM1, dir1, cname1, lentarget):
    if outM1 == "target":
        fhandle = open(dir1+"/"+tname, "a")
        finalseq = SeqRecord(sequence)
        finalseq.id = seqname
    elif outM1 == "query":
        if qname[-4::] != ".fas":
            fhandle = open(dir1+"/"+qname+".fas", "a")
        else:
            fhandle = open(dir1+"/"+qname, "a")
        finalseq = SeqRecord(sequence)
        if cname1 == False or seqname == "none":
            finalseq.id = tname
        else:
            finalseq.id = tname+namedelim+seqname
    elif outM1 == "combined":
        fhandle = open(dir1+"/combined.fas", "a")
        finalseq = SeqRecord(sequence)
        if lentarget == 1:
            finalseq.id = seqname
        else:
            finalseq.id = tname+namedelim+seqname
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
def get_sequence(inplist, seq, extractiontype, fls, trans_out1,
                    ac2, metric, metricR, keep_strand):
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
        if ac2 == "tdna-aa":
            start = (hits[0]+2)/3-1
            end = hits[1]/3
        else:
            start = hits[0]-1
            end = hits[1]
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
def process_fasta(target_db_name1, inputf1, final_table1, final_target_table1,
                    extractiontype1, flanks1, trans_out1, outM1, output_dir1,
                    cname1, append_name1, ac1, keep_strand, cols1, debugfile1):
    c1 = len(final_target_table1)
    messagefunc("searching for contigs in: "+target_db_name1+
                ", total number of contigs to extract: "+str(c1),
                cols1, debugfile1, False)
    if outM1 == "query":
        target_set = {}
        for qkey in list(x for x in final_table1):
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
                #going over supercontig indexes
                for sprcontig in range(len(final_table1[qname])):
                    #number of contigs in supercontig
                    num_contigs = len(final_table1[qname][sprcontig][1][3])
                    #looking for target in the query table
                    for t in range(num_contigs):
                        #if not gap integer
                        if type(final_table1[qname][sprcontig][1][3][t]) is not int:
                            tgt_index_name = final_table1[qname][sprcontig][1][3][t][0]
                            #process target name
                            targetname = indexer_function(tgt_index_name, None)[0]
                            #found target in the query table and
                            #   this particular version of target wasnt used
                            if targetname == seq.id and \
                            tgt_index_name not in target_set[check_name]:
                                #get the sequence
                                messagefunc(str(c1)+" EXTRACTING: contig "+
                                            targetname+", query "+qname,
                                            cols1, debugfile1)
                                s1 = get_sequence(final_table1[qname][sprcontig][1][3][t][1],
                                                    seq, extractiontype1, flanks1,
                                                    trans_out1, ac1, metric,
                                                    metricR, keep_strand)
                                print ("- EXTRACTING: final seq", s1[:10], "ranges",
                                        final_table1[qname][sprcontig][1][3][t][1],
                                        file = debugfile1)
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
                                            seqwritefunc(s1, qname, target_db_name1,
                                                            targetname, outM1, output_dir1,
                                                            cname1, append_name1)
                                        #else do nothing
                                    else:
                                        if use_suffix:
                                            seqwritefunc(s1, qname, target_db_name1,
                                                            indexer_function(targetname,
                                                                            str(seq_suffix_c)),
                                                                            outM1, output_dir1,
                                                                            cname1, append_name1)
                                        else:
                                            seqwritefunc(s1, qname, target_db_name1,
                                                        targetname, outM1, output_dir1,
                                                        cname1, append_name1)
                                    target_set[check_name].add(tgt_index_name)
                                else:
                                    #need to disable contig stitching when n is selected
                                    final_table1[qname][sprcontig][1][3][t][1] = s1
                                    messagefunc(str(c1)+" BUCKET: contig "+targetname+
                                                ", query "+qname, cols1, debugfile1)
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
                                        s1 = dumper(final_table1[qname][sprcontig][1][3],
                                                    extractiontype1, trans_out1, ac1)
                                        messagefunc(str(c1)+" EXTRACTING: bucket "+qname+
                                                    " dumped", cols1, debugfile1)
                                        print ("- EXTRACTING: final seq", s1[:10],
                                                file = debugfile1)
                                        seqwritefunc(s1, qname, target_db_name1,
                                                    "Merged_"+qname+"_supercontig_"+
                                                    str(sprcontig), outM1, output_dir1,
                                                    cname1, append_name1)
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
def estimate_survival(init_queries, init_targets, survived_queries, survived_targets,
                        ignored, filtration_table, cols, debugfile_generic1, debugfile):
    msg = "number of ignored records based on scoring threshold: "+str(ignored)
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "number of queries with OK results: "+str(len(survived_queries))
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "number of queries without OK results: "+str(len(init_queries - survived_queries))
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    print ("empty query list:", file = debugfile)
    print (" ".join(list(init_queries - survived_queries)), file = debugfile)
    msg = "number of passed hits: "+str(len(survived_targets))
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "number of filtered out hits: "+str(len(init_targets - survived_targets))
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "synteny / max gap check failed: "+str(filtration_table["synteny"])
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "per hit reciprocal match check failed (hit matches several queries): "+\
        str(filtration_table["actual"])
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "range reciprocal match check failed (stitched hit matches several queries): "+\
        str(filtration_table["range"])
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "reference reciprocal match check failed (RBH check fail): "+\
        str(filtration_table["rbh"])
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    msg = "filtered away by number of allowed supercontigs: "+\
        str(filtration_table["contig"])
    messagefunc(msg, cols, debugfile, False)
    print (msg, file = debugfile_generic1)
    # print >> debugfile, "filtered out targets:"
    # print >> debugfile, " ".join(list(init_targets - survived_targets))

#---------------------------------------------------------------------#
############################ main script ##############################
debugfile_generic = open("alibaseq_"+logsuffix+".log", "w")

#get terminal window size
if not os.popen('stty size', 'r').read():
    cols = 100
else:
    rows, cols = os.popen('stty size', 'r').read().split()
    cols = int(cols)-10

messagefunc("alibaseq run with option "+filefolder+" selected", cols,
            debugfile_generic, False)
messagefunc("command line parameters: "+' '.join(sys.argv), cols,
            debugfile_generic, False)
messagefunc("contignum parameter: "+str(contignum), cols,
            debugfile_generic, False)

#multi db option
if filefolder == "M":
    #reading the blastfile
    if bt == "blast":
        blastlist = []
        preblastlist = glob.glob(blastfilearg+"/*.blast")
        for blastfile in preblastlist:
            if "_reciprocal.blast" not in blastfile:
                blastlist.append(blastfile)
    elif bt == "lastz":
        blastlist = glob.glob(blastfilearg+"/*.lastz")
    elif bt == "sam":
        blastlist = glob.glob(blastfilearg+"/*.sam")
    elif bt == "bam":
        blastlist = glob.glob(blastfilearg+"/*.bam")
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
    messagefunc("dry run, no sample files", cols, debugfile_generic, False)
else:
    messagefunc("list of target fasta files detected (mask *.fasta):",
                cols, debugfile_generic, False)
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
    if btR == "blast":
        target_ref, ignoredRef = readblastfilefunc(target_ref_file, evalue,
                                                    bitscore, identity, False,
                                                    acR, False, cols, debugfile_generic)
        messagefunc("number of ignored records based on scoring threshold: "+
                    str(ignoredRef), cols, debugfile_generic, False)
        ref_filtration_table = {"synteny" : 0}
        target_ref = process_aux_tables(target_ref, metric, metricR, hit_ovlp, acR,
                                        max_gap, amlghitscore, metricC, ref_filtration_table,
                                        cols, debugfile_generic)
        messagefunc("synteny / max gap check failed: "+str(ref_filtration_table["synteny"]),
                    cols, debugfile_generic, False)
    elif btR == "bed":
        target_ref = readbedfilefunc(target_ref_file, cols, debugfile_generic)
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
    messagefunc("sample table "+str(b1)+" out of "+str(len(blastlist)),
                cols, debugfile_generic, False)
    #set up sample debug file
    samplelogname = b.split("/")[-1]+"_"+logsuffix+".log"
    debugfile = open(samplelogname, "w")
    messagefunc("sample log started: "+samplelogname, cols, debugfile_generic, False)
    #read alignment table
    if bt == "blast":
        #output 0 is query, 1 is target
        output, ignored = readblastfilefunc(b, evalue, bitscore, identity,
                                            True, ac, True, cols, debugfile)
    elif bt == "lastz":
        #output 0 is query, 1 is target
        output, ignored = readlastzfilefunc(b, bitscore, identity,
                                            True, True, cols, debugfile)
    elif bt == "sam":
        output, ignored = readsamformat(b, False, bitscore, samscore, cols, debugfile)
    elif bt == "bam":
        output, ignored = readsamformat(b, True, bitscore, samscore, cols, debugfile)
    else:
        output, ignored = readhmmerfilefunc(b, evalue, bitscore, bt,
                                            ac, hmmerg, cols, debugfile)
    #read reciprocal alignment table
    if rec_search != None:
        if rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast" in rec_list:
            rec_out, ignoredDummy = readblastfilefunc(rec_search.rstrip("/")+"/"+
                                                        b.split("/")[-1]+
                                                        "_reciprocal.blast",
                                                        None, None, None, False,
                                                        acr, False, cols, debugfile)
            if ref_hs:
                rec_filtration_table = {"synteny" : 0}
                rec_out = process_aux_tables(rec_out, metric, metricR, hit_ovlp, acr,
                                                max_gap, amlghitscore, metricC,
                                                rec_filtration_table, cols, debugfile)
                messagefunc("synteny / max gap check failed: "+
                            str(rec_filtration_table["synteny"]),
                            cols, debugfile_generic, False)
        elif len(rec_list) == 1:
            rec_out, ignoredDummy = readblastfilefunc(rec_list[0], None, None, None,
                                                        False, acr, False, cols, debugfile)
            if ref_hs:
                rec_filtration_table = {"synteny" : 0}
                rec_out = process_aux_tables(rec_out, metric, metricR,
                                                hit_ovlp, acr, max_gap,
                                                amlghitscore, metricC,
                                                rec_filtration_table,
                                                cols, debugfile)
                messagefunc("synteny / max gap check failed: "+
                            str(rec_filtration_table["synteny"]),
                            cols, debugfile_generic, False)
        else:
            print ("problem with finding the reciprocal search file")
            print (rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast", rec_list)
            sys.exit()
    else:
        rec_out = None

    final_table = {}
    final_target_table = {}
    filtration_table = {}

    #run target processor
    target_table = target_processor(output, local_rec, metric, metricR, hit_ovlp,
                                    recip_ovlp, ac, run_hs, max_gap, amlghitscore,
                                    metricC, bstrands, filtration_table, cols, debugfile)

    #run query processor
    final_table = query_processor(target_table, rec_out, target_ref, metric, metricR,
                                    metricC, contignum, ctg_ovlp, interstitch, hit_ovlp,
                                    ac, recip_ovlp, ref_hs, rmrecnf, srt,
                                    filtration_table, cols, debugfile)

    final_target_table = reformat_table(final_table, cols, debugfile)

    estimate_survival(init_queries, init_targets, survived_queries, survived_targets,
                        ignored, filtration_table, cols, debugfile_generic, debugfile)

    qout = open(b.split("/")[-1]+"_"+logsuffix+"_qtable.tab", "w")
    bltableout(final_table, qout, "query")
    qout.close()
    qoutH = open(b.split("/")[-1]+"_"+logsuffix+"_qtableH.tab", "w")
    bltableout(final_table, qoutH, "queryH")
    qoutH.close()
    tout = open(b.split("/")[-1]+"_"+logsuffix+"_ttable.tab", "w")
    if loghead:
        print("hit,num_queries,names_of_queries", file=tout)
    bltableout(final_target_table,tout, "target")
    tout.close()

    #####-----------------------------------------------------------------------
    if dry_run:
        messagefunc("dry run, search through sample file skipped",
                    cols, debugfile_generic, False)
    else:
        messagefunc("scanning the sample fasta file(s)...",
                    cols, debugfile_generic, False)

        if filefolder == "SM":
            for seqname in targetlist:
                inputf = SeqIO.parse(seqname, "fasta")
                target_db_name = seqname.split("/")[-1]
                messagefunc("--------scanning the sample---------",
                            cols, debugfile)
                process_fasta(target_db_name, inputf, final_table,
                                final_target_table, extractiontype, flanks,
                                trans_out,outM,output_dir, cname, append_name,
                                ac, keep_strand, cols, debugfile)
        else:
            if filefolder == "M":
                seqname = b[:-6].split("/")[-1]
                target_db_match = [target_db_m1 for target_db_m1 in targetlist
                                    if seqname in target_db_m1]
                if len(target_db_match) == 1:
                    inputf = SeqIO.parse(target_db_match[0], "fasta")
                    target_db_name = seqname.split("/")[-1]
                elif len(target_db_match) == 0:
                    msg = "error, the sample fasta file "+seqname+" is not found"
                    messagefunc(msg, cols, debugfile, False)
                else:
                    msg = "error, several matches to "+seqname+" are found in the folder"
                    messagefunc(msg, cols, debugfile, False)
            else:
                seqname = targetlist[0]
                inputf = SeqIO.parse(seqname, "fasta")
                target_db_name = seqname.split("/")[-1]
            messagefunc("--------scanning the sample---------",
                        cols, debugfile)
            process_fasta(target_db_name, inputf, final_table, final_target_table,
                            extractiontype, flanks, trans_out,outM,output_dir, cname,
                            append_name, ac, keep_strand, cols, debugfile)

    debugfile.close()

print (len(warninglist), file = debugfile_generic)
for w in warninglist:
    print (w, file = debugfile_generic)

print ("done", file = debugfile_generic)
debugfile_generic.close()

print ("done")
