from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import os
import shutil
import csv
import sys
import itertools
import argparse


parser = argparse.ArgumentParser(description='ABSX (Alignment-Based Sequence eXtraction)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-b', metavar='table', help='alignment table file',dest="blastfilearg",required=True)
required.add_argument('-ff', choices=['S','M'], help='file / folder mode',dest="filefolder", required=True)
required.add_argument('-et', choices=['n','s','a','b'], help='extraction type: n (normal), s (only best hit region), a (extract all hit regions and join them), b (extract region between two outmost hit regions)',dest="extractiontype",required=True)

optional.add_argument('-t', metavar='assembly', help='assembly file',dest="targetf")
optional.add_argument('-q', metavar='query', help='query file(s) to which extracted results are to be appended; if not specified, sequences are extracted into blank files',dest="queryf")
optional.add_argument('-o', metavar='output', help='output folder for modified files with extracted sequences',dest="output", default="absx_out")
optional.add_argument('-eval', metavar='N', help='evalue cutoff',dest="evalue", type=float, default=0.01)
optional.add_argument('-c', metavar='N', help='number of contigs to extract, if set to 0, then extract all contigs',dest="contignum", type=int, default=0)
optional.add_argument('-fl', metavar='N', help='flanks on each side in bp',dest="flanks", type=int, default=0)
optional.add_argument('--no-recip', dest='reciprocate', action='store_false', help='do not run reciprocator', default=True)
optional.add_argument('--is', dest='interstich', action='store_true', help='perform intercontig stiching', default=False)
optional.add_argument('--allow-ovlp', dest='allow_ovlp', action='store_true', help='allow hit overlap', default=False)
optional.add_argument('-bt', choices=['blast','hmmer'], help='alignment table type',dest="bt", default="blast")
#blastn, tblastx, blastp = 1 to 1; tblastn = aa in query, dna in db; blastx = dna in query, aa in db. hmmer is always 1 to 1.
optional.add_argument('-ac', choices=['normal','tblastn', 'blastx'], help='alignment coordinate type',dest="ac", default="normal")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()
else:
    args = parser.parse_args()

    blastfilearg = vars(args)["blastfilearg"]
    filefolder = vars(args)["filefolder"]
    extractiontype = vars(args)["extractiontype"]

    if vars(args)["targetf"] == None:
        dry_run = True
        noq = True
        print "option -q ignored"
    else:
        targetf = vars(args)["targetf"]
        dry_run = False
        if vars(args)["queryf"] == None:
            noq = True
        else:
            queryf = vars(args)["queryf"]
            noq = False
    
    evalue = vars(args)["evalue"]
    contignum = vars(args)["contignum"]
    reciprocate = vars(args)["reciprocate"]
    interstich = vars(args)["interstich"]
    allow_ovlp = vars(args)["allow_ovlp"]
    flanks = vars(args)["flanks"]
    output_dir = vars(args)["output"]
    bt = vars(args)["bt"]
    ac = vars(args)["ac"]

    print blastfilearg, evalue, filefolder,extractiontype,contignum,reciprocate, interstich, allow_ovlp, flanks, dry_run, noq

dash = "--------------------------------------------------------"
warninglist = []
recip_olvp = 10 #max overlap on query, after which contigs start compete
hit_query_ovl = 10 #max overlap on query, after which target hits are merged
hit_target_ovl = 0 #max overlap on target, after which target hits are merged
contig_ovlp = 0 #max overlap on query, after which contigs are considered overlapping and are not stiched
bitAVG = False #using max bitscore or average bitscore to compare blast hits


def messagefunc(msg, f, fl=True):
    if fl:
        sys.stdout.write(msg+"\r")
        sys.stdout.flush()
    else:
        print ""
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
def copyfunc(dir1):
    messagefunc("copying files *.fa*", debugfile, False)
    copyfunc_c = 0
    for x in glob.glob(queryf+"/*.fa*"):
        locusfname = x.split("/")[-1]
        #print locusfname
        if not os.path.exists (dir1+"/"+locusfname):
            prog = "copying "+str(locusfname)+"..."
            messagefunc(prog, debugfile)
            shutil.copy2(queryf+"/"+locusfname, dir1)
            copyfunc_c += 1
    messagefunc("copied "+str(copyfunc_c)+" files", debugfile, False)

#function for parsing a blast output file
def readblastfilefunc(b, debugfile):
    messagefunc("processing "+b, debugfile, False)
    querydict = {}
    targetdict = {}
    blastfile = open(b, "rU")
    reader = csv.reader(blastfile, delimiter='\t')
    linecounter = 0
    for row in reader:
        if float(row[10]) <= evalue:
            qname = row[0].split("/")[-1]
            tname = row[1]
            #populate query table
            if qname in querydict:
                if tname in querydict[qname]:
                    querydict[qname][tname][linecounter] = rowfunc(row, ac)
                else:
                    querydict[qname][tname] = {linecounter: rowfunc(row, ac)}
            else:
                querydict[qname] = {tname: {linecounter: rowfunc(row, ac)}}
            #populate target table
            if tname in targetdict:
                if qname in targetdict[tname]:
                    targetdict[tname][qname][linecounter] = rowfunc(row, ac)
                else:
                    targetdict[tname][qname] = {linecounter: rowfunc(row, ac)}
            else:
                targetdict[tname] = {qname: {linecounter: rowfunc(row, ac)}}
        linecounter += 1
    blastfile.close()
    return querydict, targetdict

#WIP
#nhmmer --tblout
# def readhmmerfilefunc(b, debugfile):
#     messagefunc("processing "+b, debugfile, False)
#     querydict = {}
#     targetdict = {}
#     hmmfile = open(b, "rU")
#     #reader = csv.reader(blastfile, delimiter='\t')
#     linecounter = 0
#     recordcounter = 0
#     for row in hmmfile:
#         if row[0] != "#":
#             line = row.strip().split()
#             #print line, line[12]
#             if float(line[12]) <= evalue:
#                 qname = line[2]+".fas"
#                 tname = line[0]
#                 #populate query table
#                 if qname in querydict:
#                     if tname in querydict[qname]:
#                         querydict[qname][tname][linecounter] = hmmrowfunc(line)
#                     else:
#                         querydict[qname][tname] = {linecounter: hmmrowfunc(line)}
#                 else:
#                     querydict[qname] = {tname: {linecounter: hmmrowfunc(line)}}
#                 #populate target table
#                 if tname in targetdict:
#                     if qname in targetdict[tname]:
#                         targetdict[tname][qname][linecounter] = hmmrowfunc(line)
#                     else:
#                         targetdict[tname][qname] = {linecounter: hmmrowfunc(line)}
#                 else:
#                     targetdict[tname] = {qname: {linecounter: hmmrowfunc(line)}}
#             linecounter += 1
#     hmmfile.close()
#     return querydict, targetdict

#hmmsearch --domtblout
def readhmmerfilefunc(b, debugfile):
    messagefunc("processing "+b, debugfile, False)
    querydict = {}
    targetdict = {}
    hmmfile = open(b, "rU")
    #reader = csv.reader(blastfile, delimiter='\t')
    linecounter = 0
    recordcounter = 0
    for row in hmmfile:
        if row[0] != "#":
            line = row.strip().split()
            #print line, line[12]
            if float(line[12]) <= evalue:
                qname = line[3]+".fas"
                tname = line[0]
                #populate query table
                if qname in querydict:
                    if tname in querydict[qname]:
                        querydict[qname][tname][linecounter] = hmmrowfunc(line)
                    else:
                        querydict[qname][tname] = {linecounter: hmmrowfunc(line)}
                else:
                    querydict[qname] = {tname: {linecounter: hmmrowfunc(line)}}
                #populate target table
                if tname in targetdict:
                    if qname in targetdict[tname]:
                        targetdict[tname][qname][linecounter] = hmmrowfunc(line)
                    else:
                        targetdict[tname][qname] = {linecounter: hmmrowfunc(line)}
                else:
                    targetdict[tname] = {qname: {linecounter: hmmrowfunc(line)}}
            linecounter += 1
    hmmfile.close()
    return querydict, targetdict


#implement reciprocator break value, default 10. DO percet, e.g. 0.1 of the shortes range
def reciprocator(inpdict, query, range1, range2, emax, bitscore, target):
    cond = True
    for key, val in inpdict.items():
        if key != query: #all other queries
            target_ranks_temp = compute_ranks(val) #get all pieces, compute values for this query
            if getOverlap([target_ranks_temp[2],target_ranks_temp[3]],[range1,range2]) > recip_olvp:
                if target_ranks_temp[0] < emax:
                    cond = False
                    #messagefunc("coords compared: "+str(target_ranks_temp[2])+":"+str(target_ranks_temp[3])+", "+str(range1)+":"+str(range2), debugfile)
                    messagefunc("ranks temp: "+",".join(map(str, target_ranks_temp)), debugfile)
                    break
                elif target_ranks_temp[0] == emax:
                    if target_ranks_temp[1] > bitscore:
                        cond = False
                        #messagefunc("coords compared: "+str(target_ranks_temp[2])+":"+str(target_ranks_temp[3])+", "+str(range1)+":"+str(range2), debugfile)
                        messagefunc("ranks temp: "+",".join(map(str, target_ranks_temp)), debugfile)
                        break
                    else:
                        wrn = "warning, target "+target+" has equal hits to several queries, saved for both!"
                        warninglist.append(wrn)
                        messagefunc(wrn, debugfile)
    return cond


def bltableout(output, bltableout_file, table_type):
    for key, value in sorted(output.items()):
        if table_type == "target":
            print >> bltableout_file, key+","+str(len(value))+","+",".join(value)
        elif table_type == "query":
            print >> bltableout_file, key, value

#function to return a list for dict like this:
#dict[query] = [target_f, target_r, target_b, query_f, query_r, query_b, eval, bitscore]
#query - query name, target_f - target start pos, target_r - target end, target_b - forward or reverse target direction
#query_f - query start, query_r - query end, query_b - query direction
#blastn, tblastx, blastp = 1 to 1; tblastn = aa in query, dna in db; blastx = dna in query, aa in db. hmmer is always 1 to 1.
#optional.add_argument('-ac', choices=['normal','tblastn', 'blastx'], help='alignment coordinate type',dest="ac", default="normal")
##### redo to divide by three in blastx case
def rowfunc(row, aligntype):
    if int(row[8]) < int(row[9]):
        target_b = True
    else:
        target_b = False
    if aligntype == "normal" or "tblastn":
        target_f = int(row[8])
        target_r = int(row[9])
    else: #blastx
        if target_b:
            target_f = int(row[8])*3-2
            target_r = int(row[9])*3
        else:
            target_f = int(row[8])*3
            target_r = int(row[9])*3-2
    #check query
    if int(row[6]) < int(row[7]):
        query_b = True
    else:
        query_b = False
    if aligntype == "normal" or "blastx":
        query_f = int(row[6])
        query_r = int(row[7])
    else: #tblastn
        if query_b:
            query_f = int(row[6])*3-2
            query_r = int(row[7])*3
        else:
            query_f = int(row[6])*3
            query_r = int(row[7])*3-2
    return [target_f, target_r, target_b, query_f, query_r, query_b, float(row[10]), float(row[11])]

# #nhmmer --tblout
# def hmmrowfunc(row):
#     if row[11] == "+":
#         target_b = True
#     else:
#         target_b = False
#     target_f = int(row[8])
#     target_r = int(row[9])
#     #check query
#     query_b = True
#     query_f = int(row[4])
#     query_r = int(row[5])
#     return [target_f, target_r, target_b, query_f, query_r, query_b, float(row[12]), float(row[13])]

#hmmsearch --domtblout   
def hmmrowfunc(row):
    target_b = True
    target_f = int(row[17])
    target_r = int(row[18])
    #check query
    query_b = True
    query_f = int(row[15])
    query_r = int(row[16])
    return [target_f, target_r, target_b, query_f, query_r, query_b, float(row[12]), float(row[13])]


#function to compute overlap between two ranges supplied as lists with start and end
#returns overlap value
def getOverlap(a, b):
    a0=min(a)
    a1=max(a)
    b0=min(b)
    b1=max(b)
    return max(0, min(a1, b1) - max(a0, b0))

#function for appending the seqeunce to the alignemnt
#returns 1 if success, 0 otherwise
def seqwritefunc(sequence, qname, tname, seqname, noq,dir1):
    if noq:
        fhandle = open(dir1+"/"+tname, "a")
        finalseq = SeqRecord(sequence)
        finalseq.id = seqname
    else:
        fhandle = open(dir1+"/"+qname, "a")
        finalseq = SeqRecord(sequence)
        if contignum == 1 or seqname == "none":
            finalseq.id = tname
        else:
            finalseq.id = tname+"|"+seqname
    finalseq.name =""
    finalseq.description =""
    SeqIO.write(finalseq, fhandle, "fasta")
    return 1
    fhandle.close()

def compute_ranks(hits):
    eval_max = []
    bitscore = []
    coord = []
    for key, hit in hits.items():
        eval_max.append(hit[6])
        bitscore.append(hit[7])
        coord.append(hit[0])
        coord.append(hit[1])
    if bitAVG:
        return [min(eval_max), float(sum(bitscore)) / max(len(bitscore), 1), min(coord), max(coord)]
    else:
        return [min(eval_max), max(bitscore), min(coord), max(coord)]

def hit_sticher(inpdict, extractiontype, ovlpB):
    outlist = []
    #get best item and its direction
    bhit = -1
    best = 1
    for key, val in inpdict.items():
        if bhit == -1:
            bhit = val[7]
        if val[6] < best or (val[6] == best and val[7] > bhit):
            best = val[6]
            bhit = val[7]
            if val[2] == val[5]:
                direct = True
            else:
                direct = False
            if extractiontype == "n" or extractiontype == "s" or len(inpdict.keys()) == 1:
                if extractiontype == "n":
                    outlist = [direct, [-1, -1, val[3], val[4], best, bhit]]
                else:
                    #Option -e is taken care of at the moment of seq extraction
                    outlist = [direct, [val[0],val[1], val[3], val[4], best, bhit]]
                    if extractiontype == "a" or extractiontype == "b":
                        messagefunc("only one hit region, abort", debugfile)
    if extractiontype == "a" and len(inpdict.keys()) > 1 or extractiontype == "b" and len(inpdict.keys()) > 1:
        messagefunc("running hit overlapper...", debugfile)
        stichlist = inpdict.values()
        ovlp = True
        while ovlp:
            if len(stichlist) == 1:
                break
            else:
                combos = list(itertools.combinations(range(len(stichlist)), 2))
                messagefunc("combinations: "+str(len(combos))+", number of regions: "+str(len(stichlist)), debugfile)
                for comb in range(len(combos)):
                    tovlp = getOverlap(stichlist[combos[comb][0]][0:2],stichlist[combos[comb][1]][0:2])
                    qovlp = getOverlap(stichlist[combos[comb][0]][3:5],stichlist[combos[comb][1]][3:5])
                    if tovlp > hit_target_ovl and abs(qovlp-tovlp) < hit_query_ovl: #hits overlap on target and roughly same way overlap on query
                        #stichlist[combos[comb][0]] # this is whole record with 8 elements
                        stichlist[combos[comb][0]] = [min(stichlist[combos[comb][0]][0], stichlist[combos[comb][1]][0],stichlist[combos[comb][0]][1], stichlist[combos[comb][1]][1]), max(stichlist[combos[comb][0]][0], stichlist[combos[comb][1]][0],stichlist[combos[comb][0]][1], stichlist[combos[comb][1]][1]), direct, min(stichlist[combos[comb][0]][3], stichlist[combos[comb][1]][3],stichlist[combos[comb][0]][4], stichlist[combos[comb][1]][4]), max(stichlist[combos[comb][0]][3], stichlist[combos[comb][1]][3],stichlist[combos[comb][0]][4], stichlist[combos[comb][1]][4]),direct,min(stichlist[combos[comb][0]][6], stichlist[combos[comb][1]][6]),max(stichlist[combos[comb][0]][7], stichlist[combos[comb][1]][7])]
                        del stichlist[combos[comb][1]]
                        messagefunc("overlapped, merging...", debugfile)
                        
                        ovlp = True
                        break
                    elif tovlp <= hit_target_ovl and abs(qovlp-tovlp) >= hit_query_ovl: #hits barely overlap on target but a lot on query -- possibly repeated target region? Select best based on scores
                        messagefunc("warning: repeated hit region?", debugfile)
                        print >> debugfile, stichlist[combos[comb][0]], stichlist[combos[comb][1]]
                        if ovlpB:
                            ovlp = False
                        else:
                            if stichlist[combos[comb][0]][6] < stichlist[combos[comb][1]][6] or stichlist[combos[comb][0]][7] > stichlist[combos[comb][1]][7] or abs(stichlist[combos[comb][0]][0]-stichlist[combos[comb][0]][1]) >= abs(stichlist[combos[comb][1]][0]-stichlist[combos[comb][1]][1]):
                                del stichlist[combos[comb][1]]
                            else:
                                del stichlist[combos[comb][0]]
                            ovlp = True
                            break
                    elif (float(tovlp) / abs(stichlist[combos[comb][0]][0]-stichlist[combos[comb][0]][1]))*100 > 10 and qovlp == 0: #hits overlapping on target but not on query - remove as well as before
                        messagefunc("warning: disjunct query hits, deleting...", debugfile)
                        print >> debugfile, stichlist[combos[comb][0]], stichlist[combos[comb][1]], (float(tovlp) / abs(stichlist[combos[comb][0]][0]-stichlist[combos[comb][0]][1]))*100
                        if stichlist[combos[comb][0]][6] < stichlist[combos[comb][1]][6] or stichlist[combos[comb][0]][7] > stichlist[combos[comb][1]][7] or abs(stichlist[combos[comb][0]][0]-stichlist[combos[comb][0]][1]) >= abs(stichlist[combos[comb][1]][0]-stichlist[combos[comb][1]][1]):
                            del stichlist[combos[comb][1]]
                        else:
                            del stichlist[combos[comb][0]]
                        ovlp = True
                        break
                    else:
                        ovlp = False
                if not ovlp:
                    messagefunc("no more overlaps", debugfile)
        #RUN STICHER and return margins and also sequence of regions and gaps
        messagefunc("running hit sticher...", debugfile)
        median_coords = {}
        start_coords = {}
        end_coords = {}
        start_target = {}
        end_target = {}
        for chunk in range(len(stichlist)):
            median_coords[chunk] = median([stichlist[chunk][3],stichlist[chunk][4]])
            start_coords[chunk] = min(stichlist[chunk][3],stichlist[chunk][4])
            end_coords[chunk] = max(stichlist[chunk][3],stichlist[chunk][4])
            start_target[chunk] = min(stichlist[chunk][0],stichlist[chunk][1])
            end_target[chunk] = max(stichlist[chunk][0],stichlist[chunk][1])
        gapstart = 0
        outlist = [direct]
        for key in sorted(median_coords, key=lambda x: median_coords[x]):
            if gapstart > 0:
                outlist.append(start_coords[key]-gapstart)
            outlist.append([start_target[key], end_target[key], start_coords[key], end_coords[key]])
            gapstart = end_coords[key]
    return outlist

def median(lst): #taken from https://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

def contig_overlap(inplist):
    tab = {}
    for target in inplist:
        tab[target[0]] = [min(target[1][1][2],target[1][1][3]),max(target[1][-1][2],target[1][-1][3])]
    flatlist = tab.values()
    messagefunc("checking contig overlap", debugfile)
    
    combos = list(itertools.combinations(range(len(flatlist)), 2))
    ovlp = False
    for comb in range(len(combos)):
        if getOverlap(flatlist[combos[comb][0]][:2],flatlist[combos[comb][1]][:2]) > contig_ovlp:
            messagefunc("overlapping: "+str(getOverlap(flatlist[combos[comb][0]][:2],flatlist[combos[comb][1]][:2])), debugfile)
            ovlp = True
            break
    if ovlp:
        return True
    else:
        return False

def contig_sticher(inplist):
    messagefunc("running contig sticher...", debugfile)
    median_coords = {}
    start_coords = {}
    end_coords = {}
    for target in inplist:
        median_coords[target[0]] = median([min(target[1][1][2],target[1][1][3]),max(target[1][-1][2],target[1][-1][3])])
        start_coords[target[0]] = min(target[1][1][2],target[1][1][3])
        end_coords[target[0]] = max(target[1][-1][2],target[1][-1][3])
    gapstart = 0
    output = []
    for key in sorted(median_coords, key=lambda x: median_coords[x]):
        if gapstart > 0:
            output.append(start_coords[key]-gapstart)
        output.append(key)
        gapstart = end_coords[key]
    return output

def get_sequence(inplist, seq, extractiontype, fls):
    finalseq = Seq("")
    seqlen = len(seq.seq)
    if extractiontype == "a":
        for i in inplist[1:]:
            if type(i) is not int:
                start = min(i[0],i[1])-1
                end = max(i[0],i[1])
                if inplist[0]:
                    finalseq += seq.seq[start:end]
                else:
                    finalseq += seq.seq[start:end].reverse_complement()
            else:
                if i > 0:
                    finalseq += Seq("N"*i)
                else:
                    finalseq += Seq("N")
    elif extractiontype == "b":
        f1 = True
        for i in inplist[1:]:
            if type(i) is not int:
                if f1:
                    start = min(i[:2])
                    end = max(i[:2])
                    f1 = False
                else:
                    start = min(start, i[0], i[1])
                    end = max(end, i[0], i[1])
        start = start - 1
        #fls implementation
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
        start = min(inplist[1][0],inplist[1][1])-1
        end = max(inplist[1][0],inplist[1][1])
        #fls implementation
        if start - fls < 0:
            start = 0
        else:
            start = start - fls
        if end + fls > seqlen:
            end = seqlen
        else:
            end = end + fls
        finalseq = seq.seq[start:end]
    # elif extractiontype[:2] == "-e":
    #     flank = int(extractiontype[2:])
    #     start = min(inplist[1][0],inplist[1][1])-1
    #     end = max(inplist[1][0],inplist[1][1])
    #     if start - flank < 0:
    #         start = 0
    #     else:
    #         start = start - flank
    #     if end + flank > len(seq.seq)-1:
    #         end = len(seq.seq)-1
    #     else:
    #         end = end + flank
    #     finalseq = seq.seq[start:end]
    if not inplist[0] and extractiontype != "a":
        finalseq = finalseq.reverse_complement()
    return finalseq

def dumper(inplist, extractiontype):
    finalseq = Seq("")
    for i in inplist:
        if type(i) is not int:
            finalseq += i
        else:
            if i > 0 and extractiontype != "n":
                finalseq += Seq("N"*i)
            else:
                finalseq += Seq("N")
    return finalseq
#---------------------------------------------------------------------

debugfile = open("absx.log", "w")

messagefunc("absx run with option "+filefolder+" selected", debugfile, False)
messagefunc("command line parameters: "+' '.join(sys.argv), debugfile, False)


#make modified dir
messagefunc("make modified dir...", debugfile, False)
mkdirfunc(output_dir)

#copy files
if not noq:
    messagefunc("copy files...", debugfile, False)
    copyfunc(output_dir)

#multi db option
if filefolder == "M":
    #reading the blastfile
    if bt == "blast":
        blastlist = glob.glob(blastfilearg+"/*.blast")
    else:
        blastlist = glob.glob(blastfilearg+"/*.hmmer")
    translist = glob.glob(targetf+"/*.fasta")
elif filefolder == "S":
    blastlist = [blastfilearg]
    translist = [targetf]

messagefunc("list of target fasta files detected (mask *.fasta):", debugfile, False)
for l in translist:
    messagefunc(l, debugfile)

#debug vars
number = 0
numberset = set()
totalloci = 0
#parsing blast files
messagefunc("parsing blast files...", debugfile, False)

b1 = 0
for b in blastlist:
    b1 += 1
    messagefunc("target "+str(b1)+" out of "+str(len(blastlist)), debugfile, False)
    if bt == "blast":
        output = readblastfilefunc(b, debugfile) #output 0 is query, 1 is target
    else:
        output = readhmmerfilefunc(b, debugfile)
    final_table = {}
    final_target_table = {}
    for query in output[0].keys():
        #QUERY PROCESSING: first, rank targets by highest eval, also get average bitscore
        print >> debugfile, dash
        messagefunc("Q: "+query, debugfile)
        ranks = [{},{}] #evail is first, bitscore is second
        messagefunc("computing initial target contig stats...", debugfile)
        for target, hits in output[0][query].items():
            ranks_temp = compute_ranks(hits)
            #check reciprocy
            if reciprocate == False or reciprocate == True and reciprocator(output[1][target], query, ranks_temp[2], ranks_temp[3], ranks_temp[0],ranks_temp[1], target):
                ranks[0][target] = ranks_temp[0]
                ranks[1][target] = ranks_temp[1]
            else:
                messagefunc("reciprocator: target "+target+" removed from query "+query,  debugfile)
                
        if len(ranks[0]) == 0:
            messagefunc("EMPTY "+query,  debugfile)
        else:
            sorted_evals = sorted(ranks[0], key=lambda x: ranks[0][x])
            sorted_bits = sorted(ranks[1], key=lambda x: ranks[1][x], reverse=True)
            #GET TARGETS ORDERED AND STICHED
            targets = []
            messagefunc("sorting target contigs...", debugfile)
            wrn1 = True
            for x in range(len(ranks[0])): #using length of ranks, since some contigs are removed due to better hit elsewhere
                if sorted_evals[0] != sorted_bits[0] and wrn1:
                    messagefunc("eval and bit disagree at rank "+str(x+1)+" out of "+str(len(ranks[0])), debugfile)
                    wrn1 = False
                    if x == 0:
                        wrn = "warning, eval and bit disagree for query "+query+", sorted based on evals, "+str(ranks[0][sorted_evals[0]])+" "+str(ranks[1][sorted_bits[0]])
                        warninglist.append(wrn)
                        messagefunc(wrn, debugfile)
                tname1 = sorted_evals.pop(0)
                del sorted_bits[0]
                #SELECT OPTION:
                targets.append([tname1, hit_sticher(output[0][query][tname1], extractiontype, allow_ovlp)])
            messagefunc("best matching contig: "+targets[0][0]+", total contigs: "+str(len(targets)), debugfile)
            #print >> debugfile, targets
            #CHECK TARGETS FOR OVERLAP
            if interstich:
                if len(targets) > 1:
                    if contig_overlap(targets):
                        messagefunc("contigs overlapping, no contig stiching", debugfile)        
                        stiching_schedule = "none"
                        #CUTTING OFF EXCESS CONTIGS
                        if len(targets) > contignum and contignum > 0:
                            targets = targets[:contignum]
                    else:
                        stiching_schedule = contig_sticher(targets)
                        #ALL will be stiched to just one
                else:
                    messagefunc("only 1 target, no contig stiching", debugfile)
                    stiching_schedule = "none"
            else:
                messagefunc("IS deactivated, no contig stiching", debugfile)
                stiching_schedule = "none"
                if len(targets) > contignum and contignum > 0:
                    targets = targets[:contignum]

        final_table[query] = [targets, stiching_schedule]
        for t in targets:
            if t[0] in final_target_table:
                final_target_table[t[0]].append(query)
            else:
                final_target_table[t[0]] = [query]
    qout = open(b.split("/")[-1]+"_qtable.tab", "w")
    tout = open(b.split("/")[-1]+"_ttable.tab", "w")
    bltableout(final_table,qout, "query")
    bltableout(final_target_table,tout, "target")
    qout.close()
    tout.close()

#####-----------------------------------------------------------------------

    messagefunc("scanning the target fasta file...", debugfile, False)
    
    #get the transcriptome filename, matching blast filename
    for t_file in translist:
        if b[:-6].split("/")[-1] in t_file:
            
            inputf = SeqIO.parse(t_file, "fasta")
            
            seqname = b[:-6].split("/")[-1]
            target_db_name = t_file.split("/")[-1]
            break
    #print >> debugfile, "target:", target_db_name, "; target name:", seqname
    if not inputf:
        messagefunc("error, the target fasta file is not found", debugfile, False)
        break
    c1 = len(final_target_table)
    messagefunc("searching for contigs in: "+target_db_name+", total number of contigs: "+str(c1), debugfile, False)
    for seq in inputf: #going over seqs in target file
        if seq.id in final_target_table: #looking up same seq in target file
            for qname in final_target_table[seq.id]: #checking it's queries
                for t in range(len(final_table[qname][0])): #looking for target in the query table
                    if final_table[qname][0][t][0] == seq.id: #found target in the query table
                        if final_table[qname][1] == "none":
                            #extraction
                            messagefunc(str(c1)+" EXTRACTING: contig "+final_table[qname][0][t][0]+", query "+qname, debugfile)
                            s1 = get_sequence(final_table[qname][0][t][1], seq, extractiontype, flanks)
                            print >> debugfile, "- EXTRACTING: final seq", s1[:10], "ranges", final_table[qname][0][t][1]
                            seqwritefunc(s1, qname,target_db_name, seq.id, noq, output_dir)
                        else:
                            s1 = get_sequence(final_table[qname][0][t][1], seq, extractiontype, flanks)
                            final_table[qname][1][final_table[qname][1].index(final_table[qname][0][t][0])] = s1
                            messagefunc(str(c1)+" BUCKET: contig "+final_table[qname][0][t][0]+", query "+qname, debugfile)
                            dump_bucket = True
                            for buck1 in final_table[qname][1]:
                                if type(buck1) is str:
                                    dump_bucket = False
                                    break
                            if dump_bucket:
                                s1 = dumper(final_table[qname][1], extractiontype)
                                messagefunc(str(c1)+" EXTRACTING: bucket "+qname+" dumped", debugfile)
                                print >> debugfile, "- EXTRACTING: final seq", s1[:10]#, "ranges", final_table[qname][1]
                                seqwritefunc(s1, qname,target_db_name, "Merged_"+qname, noq,output_dir)
                        #cleanup
                        del final_table[qname][0][t]
                        break # breaking from target table
            #clean up after all qs are done
            del final_target_table[seq.id]
            c1 = c1 - 1
        if len(final_target_table) == 0:
            messagefunc(str(c1)+" search finished", debugfile, False)
            break

print len(warninglist)
for w in warninglist:
    print w

print >> debugfile, "done"
debugfile.close()

print "done"