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
# optional.add_argument('--allow-ovlp', dest='allow_ovlp', action='store_true', help='allow hit overlap', default=False)
optional.add_argument('--translate', dest='trans_out', action='store_true', help='translate output (for -et s or -et a)', default=False)
optional.add_argument('--hit-ovlp', metavar='N', help='allowed hit overlap on query, in bp',dest="hit_ovlp", type=int, default=5)
optional.add_argument('--ctg-ovlp', metavar='N', help='allowed contig overlap on query, in bp',dest="ctg_ovlp", type=int, default=1)
optional.add_argument('--recip-ovlp', metavar='N', help='contig overlap on query for reciprocator selection, in bp',dest="recip_ovlp", type=int, default=10)
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
    if (extractiontype == "s" or extractiontype == "a") and vars(args)["trans_out"]:
        trans_out = True
    else:
        trans_out = False

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
    # allow_ovlp = vars(args)["allow_ovlp"]
    hit_ovlp = vars(args)["hit_ovlp"]
    ctg_ovlp = vars(args)["ctg_ovlp"]
    recip_ovlp = vars(args)["recip_ovlp"]
    flanks = vars(args)["flanks"]
    output_dir = vars(args)["output"]
    bt = vars(args)["bt"]
    ac = vars(args)["ac"]

    print blastfilearg, evalue, filefolder,extractiontype,contignum,reciprocate, interstich, flanks, dry_run, noq

dash = "--------------------------------------------------------"
warninglist = []

def messagefunc(msg, columns, f, fl=True):
    if fl:
        sys.stdout.write((" "*columns)+"\r")
        sys.stdout.write(msg[:columns]+"\r")
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
def readblastfilefunc(b, cols, debugfile):
    messagefunc("processing "+b, cols, debugfile, False)
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
def readhmmerfilefunc(b, cols, debugfile):
    messagefunc("processing "+b, cols, debugfile, False)
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
def reciprocator(inpdict, query, range1, range2, emax, bitscore, target, cols, debugfile, recip_overlap):
    messagefunc("running reciprocator on contig "+target, cols, debugfile)
    cond = True
    for key, val in inpdict.items():
        if key != query: #all other queries
            target_ranks_temp = compute_ranks(val) #get all pieces, compute values for this query
            if getOverlap([target_ranks_temp[2],target_ranks_temp[3]],[range1,range2]) > recip_overlap:
                if target_ranks_temp[0] < emax:
                    cond = False
                    #messagefunc("coords compared: "+str(target_ranks_temp[2])+":"+str(target_ranks_temp[3])+", "+str(range1)+":"+str(range2), debugfile)
                    # messagefunc("ranks temp: "+",".join(map(str, target_ranks_temp)), cols, debugfile)
                    messagefunc(target+" has better eval hit to "+key, cols, debugfile)
                    break
                elif target_ranks_temp[0] == emax:
                    if target_ranks_temp[1] > bitscore:
                        cond = False
                        #messagefunc("coords compared: "+str(target_ranks_temp[2])+":"+str(target_ranks_temp[3])+", "+str(range1)+":"+str(range2), debugfile)
                        # messagefunc("ranks temp: "+",".join(map(str, target_ranks_temp)), cols, debugfile)
                        messagefunc(target+" has better bitscore hit to "+key, cols, debugfile)
                        break
                    elif target_ranks_temp[1] == bitscore:
                        wrn = "warning, target "+target+" at query "+query+" has equal hits to query "+key+", saved for both!"
                        warninglist.append(wrn)
                        messagefunc(wrn, cols, debugfile)
                        messagefunc("match to current query: emax "+str(emax)+", bitmax "+str(bitscore), cols, debugfile)
                        messagefunc("match to "+key+": "+",".join(map(str, target_ranks_temp)), cols, debugfile)
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
    return [target_f, target_r, target_b, query_f, query_r, query_b, float(row[10]), float(row[11]),float(row[2])]

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
    ident = []
    coord = []
    for key, hit in hits.items():
        eval_max.append(hit[6])
        bitscore.append(hit[7])
        ident.append(hit[8])
        coord.append(hit[0])
        coord.append(hit[1])
    return [min(eval_max), max(bitscore), min(coord), max(coord), max(ident)]

def hit_sticher(inpdict, extractiontype, hit_overlap, cols, debugfile):
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
                    outlist = [direct, [val[0],val[1], val[3], val[4], best, bhit]]
                    if extractiontype == "a" or extractiontype == "b":
                        messagefunc("only one hit region, abort", cols, debugfile)
    if extractiontype == "a" and len(inpdict.keys()) > 1 or extractiontype == "b" and len(inpdict.keys()) > 1:
        messagefunc("running hit overlapper...", cols, debugfile)
        stichlist = inpdict.values()
        messagefunc("best direction: "+str(direct), cols, debugfile)
        for x in stichlist:
            #print x
            if x[2] == direct:
                messagefunc(" ".join([str(a) for a in x]), cols, debugfile)
            else:
                messagefunc(" ".join([str(a) for a in x])+" exclude", cols, debugfile)
        stichlist = [x for x in stichlist if x[2] == direct]
        messagefunc("hits survived: "+str(len(stichlist)), cols, debugfile)
        ovlp = True
        while ovlp:
            if len(stichlist) == 1:
                break
            else:
                combos = list(itertools.combinations(range(len(stichlist)), 2))
                messagefunc("combinations: "+str(len(combos))+", number of regions: "+str(len(stichlist)), cols, debugfile)
                for comb in range(len(combos)):
                    tovlp = getOverlap(stichlist[combos[comb][0]][0:2],stichlist[combos[comb][1]][0:2])
                    qovlp = getOverlap(stichlist[combos[comb][0]][3:5],stichlist[combos[comb][1]][3:5])
                    if tovlp == 0:
                        if qovlp == 0:
                            # print "both non ovlp, stich, set bool to false"
                            ovlp = False
                        elif 0 < qovlp <= hit_overlap:
                            # print "mild ovlp on q, stich, set bool to false"
                            ovlp = False
                        elif qovlp > hit_overlap:
                            # print "overlap on q is too large, remove the worst, set bool to true"
                            stichlist = remove_worst(stichlist, combos, comb, debugfile)
                            ovlp = True
                            break
                    else:
                        if tovlp == qovlp:
                            #print "both same overlap (compare with expansion of q and target), overlap, set bool to true"
                            # print "check frames of two hits, if in frame, then keep, else remove the worst"
                            startshift = abs(min(stichlist[combos[comb][0]][0], stichlist[combos[comb][0]][1])-min(stichlist[combos[comb][1]][0], stichlist[combos[comb][1]][1]))-1
                            if startshift > 0 and startshift%3 == 0:
                                # print "in frame, keep"
                                stichlist[combos[comb][0]] = [min(stichlist[combos[comb][0]][0], stichlist[combos[comb][1]][0],stichlist[combos[comb][0]][1], stichlist[combos[comb][1]][1]), max(stichlist[combos[comb][0]][0], stichlist[combos[comb][1]][0],stichlist[combos[comb][0]][1], stichlist[combos[comb][1]][1]), direct, min(stichlist[combos[comb][0]][3], stichlist[combos[comb][1]][3],stichlist[combos[comb][0]][4], stichlist[combos[comb][1]][4]), max(stichlist[combos[comb][0]][3], stichlist[combos[comb][1]][3],stichlist[combos[comb][0]][4], stichlist[combos[comb][1]][4]),(not direct),min(stichlist[combos[comb][0]][6], stichlist[combos[comb][1]][6]),max(stichlist[combos[comb][0]][7], stichlist[combos[comb][1]][7])]
                                del stichlist[combos[comb][1]]
                            else:
                                # print "not in frame, pick best"
                                stichlist = remove_worst(stichlist, combos, comb, debugfile)
                            ovlp = True
                            break
                        else:
                            # print "target ovlp is larger than query, remove the worst, set bool to true"
                            stichlist = remove_worst(stichlist, combos, comb, debugfile)
                            ovlp = True
                            break

                if not ovlp:
                    messagefunc("no more overlaps", cols, debugfile)

        #RUN STICHER and return margins and also sequence of regions and gaps
        messagefunc("running hit sticher...", cols, debugfile)
        for x in stichlist:
            print >> debugfile, x
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

def remove_worst(inplist1, combos1, comb1, debugfile):
    if inplist1[combos1[comb1][0]][6] < inplist1[combos1[comb1][1]][6] or inplist1[combos1[comb1][0]][7] > inplist1[combos1[comb1][1]][7] or abs(inplist1[combos1[comb1][0]][0]-inplist1[combos1[comb1][0]][1]) >= abs(inplist1[combos1[comb1][1]][0]-inplist1[combos1[comb1][1]][1]):
        print >> debugfile, inplist1[combos1[comb1][1]], "deleted"
        del inplist1[combos1[comb1][1]]
    else:
        print >> debugfile, inplist1[combos1[comb1][0]], "deleted"
        del inplist1[combos1[comb1][0]]
    return inplist1

def median(lst): #taken from https://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

def contig_overlap(inplist, ranksd, ctnum, cols, debugfile, contig_overlap):
    tab = {}
    messagefunc("running contig overlapper...", cols, debugfile)
    for target in inplist:
        tab[target[0]] = [min(target[1][1][2],target[1][1][3]),max(target[1][-1][2],target[1][-1][3])]
        # messagefunc(target[0]+" "+",".join([str(a) for a in tab[target[0]]]), cols, debugfile)
    
    if ctnum != 1:
        flatlist = tab.values()
        contnames = tab.keys()
        combos = list(itertools.combinations(range(len(flatlist)), 2))
        ovlp = False
        for comb in range(len(combos)):
            if getOverlap(flatlist[combos[comb][0]][:2],flatlist[combos[comb][1]][:2]) > contig_overlap:
                messagefunc(contnames[combos[comb][0]]+" and "+contnames[combos[comb][1]]+" overlapping: "+str(getOverlap(flatlist[combos[comb][0]][:2],flatlist[combos[comb][1]][:2])), cols, debugfile)
                ovlp = True
                break
        if ovlp:
            return (True, inplist)
        else:
            return (False, inplist)
    else:
        sorted_tabkeys = sorted(tab, key=lambda x: abs(tab[x][0]-tab[x][1]), reverse=True)
        bad_tabkeys = set()
        for tkn in range(len(sorted_tabkeys)):
            for tkb in range(tkn+1, len(sorted_tabkeys)):
                if sorted_tabkeys[tkn] not in bad_tabkeys and sorted_tabkeys[tkb] not in bad_tabkeys:
                    ovlpcalc = getOverlap(tab[sorted_tabkeys[tkn]][:2],tab[sorted_tabkeys[tkb]][:2])
                    if ovlpcalc > contig_overlap:
                        if ranksd[0][sorted_tabkeys[tkn]] <= ranksd[0][sorted_tabkeys[tkb]] and ranksd[1][sorted_tabkeys[tkn]] > ranksd[1][sorted_tabkeys[tkb]]:
                            bad_tabkeys.add(sorted_tabkeys[tkb])
                            messagefunc(sorted_tabkeys[tkb]+" removed, overlapping with "+sorted_tabkeys[tkn]+", worse at both", cols, debugfile)
                        elif ranksd[0][sorted_tabkeys[tkn]] >= ranksd[0][sorted_tabkeys[tkb]] and ranksd[1][sorted_tabkeys[tkn]] < ranksd[1][sorted_tabkeys[tkb]]:
                            bad_tabkeys.add(sorted_tabkeys[tkn])
                            messagefunc(sorted_tabkeys[tkn]+" removed, overlapping with "+sorted_tabkeys[tkb]+", worse at both", cols, debugfile)
                        else:
                            if ranksd[2][sorted_tabkeys[tkn]] >= ranksd[2][sorted_tabkeys[tkb]]:
                                bad_tabkeys.add(sorted_tabkeys[tkb])
                                messagefunc(sorted_tabkeys[tkb]+" removed, overlapping with "+sorted_tabkeys[tkn]+", based on ident", cols, debugfile)
                            elif ranksd[2][sorted_tabkeys[tkn]] < ranksd[2][sorted_tabkeys[tkb]]:
                                bad_tabkeys.add(sorted_tabkeys[tkn])
                                messagefunc(sorted_tabkeys[tkn]+" removed, overlapping with "+sorted_tabkeys[tkb]+", based on ident", cols, debugfile)
        tab_out = []
        for target in inplist:
            if target[0] not in bad_tabkeys:
                tab_out.append(target)
        return (False, tab_out)


def contig_sticher(inplist, cols, debugfile):
    messagefunc("running contig sticher...", cols, debugfile)
    messagefunc("number of contigs: "+str(len(inplist)), cols, debugfile)
    median_coords = {}
    start_coords = {}
    end_coords = {}
    for target in inplist:
        messagefunc(target[0], cols, debugfile)
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

def get_sequence(inplist, seq, extractiontype, fls, trans_out1):
    finalseq = Seq("")
    seqlen = len(seq.seq)
    if extractiontype == "a":
        for i in inplist[1:]:
            if trans_out1:
                #translation
                if type(i) is not int:
                    start = min(i[0],i[1])-1
                    end = max(i[0],i[1])
                    # messagefunc("len to translate: "+str(len(seq.seq[start:end])), debugfile)
                    if inplist[0]:
                        # if len(seq.seq[start:end]) % 3 != 0:
                        #     ## check which frame is best
                        #     stops = []
                        #     for frame in range(3):
                        #         stops.append(seq.seq[(start+frame):end].translate().count("*"))
                        #     messagefunc("min stop frame: "+str(min(stops)),debugfile)
                        finalseq += seq.seq[start:end].translate()
                    else:
                        finalseq += seq.seq[start:end].reverse_complement().translate()
                else:
                    if i > 0:
                        finalseq += Seq("X"*i)
                    else:
                        finalseq += Seq("X")
            else:
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

        if trans_out1:
            #translation
            if inplist[0]:
                finalseq += seq.seq[start:end].translate()
            else:
                finalseq += seq.seq[start:end].reverse_complement().translate()
        else:
            if inplist[0]:
                finalseq += seq.seq[start:end]
            else:
                finalseq += seq.seq[start:end].reverse_complement()
    if not inplist[0] and extractiontype != "a" and extractiontype != "s":
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

debugfile_generic = open("absx.log", "w")
#get terminal window size
rows, cols = os.popen('stty size', 'r').read().split()
cols = int(cols)-10

messagefunc("absx run with option "+filefolder+" selected", cols, debugfile_generic, False)
messagefunc("command line parameters: "+' '.join(sys.argv), cols, debugfile_generic, False)



#make modified dir
if not dry_run:
    messagefunc("make modified dir...", cols, debugfile_generic, False)
    mkdirfunc(output_dir)

#copy files
if not noq:
    messagefunc("copy files...", cols, debugfile_generic, False)
    copyfunc(output_dir, cols, debugfile_generic)

#multi db option
if filefolder == "M":
    #reading the blastfile
    if bt == "blast":
        blastlist = glob.glob(blastfilearg+"/*.blast")
    else:
        blastlist = glob.glob(blastfilearg+"/*.hmmer")
    if not dry_run:
        translist = glob.glob(targetf+"/*.fasta")
elif filefolder == "S":
    blastlist = [blastfilearg]
    if not dry_run:
        translist = [targetf]

if dry_run:
    messagefunc("dry run, no target files", cols, debugfile_generic, False)
else:
    messagefunc("list of target fasta files detected (mask *.fasta):", cols, debugfile_generic, False)
    for l in translist:
        messagefunc(l, cols, debugfile_generic)

#debug vars
number = 0
numberset = set()
totalloci = 0
#parsing blast files
messagefunc("parsing blast files...", cols, debugfile_generic, False)

b1 = 0
for b in blastlist:
    b1 += 1
    messagefunc("target "+str(b1)+" out of "+str(len(blastlist)), cols, debugfile_generic, False)
    debugfile = open(b.split("/")[-1]+"_absx.log", "w")
    messagefunc("target log started: "+b, cols, debugfile_generic, False)
    if bt == "blast":
        output = readblastfilefunc(b, cols, debugfile) #output 0 is query, 1 is target
    else:
        output = readhmmerfilefunc(b, cols, debugfile)
    final_table = {}
    final_target_table = {}
    for query in output[0].keys():
        #QUERY PROCESSING: first, rank targets by highest eval, also get average bitscore
        print >> debugfile, dash
        messagefunc("Q: "+query, cols, debugfile)
        ranks = [{},{},{}] #evail is first, bitscore is second, ident is third
        messagefunc("computing initial target contig stats...", cols, debugfile)
        for target, hits in output[0][query].items():
            ranks_temp = compute_ranks(hits)
            #check reciprocy
            if reciprocate == False or reciprocate == True and reciprocator(output[1][target], query, ranks_temp[2], ranks_temp[3], ranks_temp[0],ranks_temp[1], target, cols, debugfile, recip_ovlp):
                ranks[0][target] = ranks_temp[0]
                ranks[1][target] = ranks_temp[1]
                ranks[2][target] = ranks_temp[4]
            else:
                messagefunc("reciprocator: target "+target+" removed from query "+query, cols, debugfile)
                
        if len(ranks[0]) == 0:
            messagefunc("EMPTY "+query, cols, debugfile)
        else:
            sorted_evals = sorted(ranks[0], key=lambda x: ranks[0][x])
            sorted_bits = sorted(ranks[1], key=lambda x: ranks[1][x], reverse=True)
            #GET TARGETS ORDERED AND STICHED
            targets = []
            messagefunc("sorting target contigs...", cols, debugfile)
            wrn1 = True
            for x in range(len(ranks[0])): #using length of ranks, since some contigs are removed due to better hit elsewhere
                if sorted_evals[0] != sorted_bits[0]:
                    if wrn1:
                        messagefunc("eval and bit disagree at rank "+str(x+1)+" out of "+str(len(ranks[0])), cols, debugfile)
                        wrn1 = False
                    if ranks[2][sorted_evals[0]] >= ranks[2][sorted_bits[0]]:
                        tname1 = sorted_evals.pop(0)
                        sorted_bits.remove(tname1)
                    else:
                        tname1 = sorted_bits.pop(0)
                        sorted_evals.remove(tname1)
                else:
                    tname1 = sorted_evals.pop(0)
                    sorted_bits.remove(tname1)
                #SELECT OPTION:
                messagefunc("run hitsticher function on contig "+tname1, cols, debugfile)
                targets.append([tname1, hit_sticher(output[0][query][tname1], extractiontype, hit_ovlp, cols, debugfile)])
            messagefunc("best matching contig: "+targets[0][0]+", total contigs: "+str(len(targets)), cols, debugfile)
            #print >> debugfile, targets
            #CHECK TARGETS FOR OVERLAP
            if interstich:
                if len(targets) > 1:
                    ovlp_bin, targets = contig_overlap(targets, ranks, contignum, cols, debugfile, ctg_ovlp)
                    if ovlp_bin:
                        #CONTIGS OVERLAP
                        messagefunc("contigs overlapping, no contig stiching", cols, debugfile)
                        stiching_schedule = "none"
                        #CUTTING OFF EXCESS CONTIGS
                        if len(targets) > contignum and contignum > 0:
                            targets = targets[:contignum]
                    else:
                        #CONTIGS DON'T OVERLAP
                        if len(targets) == 1:
                            #They don't because only one survived overlap, as a result of contignum == 1
                            messagefunc("single contig, no contig stiching", cols, debugfile)        
                            stiching_schedule = "none"
                        else:
                            #many survived
                            stiching_schedule = contig_sticher(targets, cols, debugfile)
                            #ALL will be stiched to just one
                else:
                    messagefunc("only 1 target, no contig stiching", cols, debugfile)
                    stiching_schedule = "none"
            else:
                messagefunc("IS deactivated, no contig stiching", cols, debugfile)
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
    if dry_run:
        messagefunc("dry run, search through target file skipped", cols, debugfile_generic, False)
    else:
        messagefunc("scanning the target fasta file...", cols, debugfile_generic, False)
        messagefunc("--------scanning the target---------", cols, debugfile)
        
        #get the transcriptome filename, matching blast filename
        for t_file in translist:
            if b[:-6].split("/")[-1] in t_file:
                
                inputf = SeqIO.parse(t_file, "fasta")
                
                seqname = b[:-6].split("/")[-1]
                target_db_name = t_file.split("/")[-1]
                break
        #print >> debugfile, "target:", target_db_name, "; target name:", seqname
        if not inputf:
            messagefunc("error, the target fasta file is not found", cols, debugfile, False)
            break
        c1 = len(final_target_table)
        messagefunc("searching for contigs in: "+target_db_name+", total number of contigs: "+str(c1), cols, debugfile, False)
        if noq:
            target_set = set()
        for seq in inputf: #going over seqs in target file
            if seq.id in final_target_table: #looking up same seq in target file
                for qname in final_target_table[seq.id]: #checking it's queries
                    for t in range(len(final_table[qname][0])): #looking for target in the query table
                        if final_table[qname][0][t][0] == seq.id: #found target in the query table
                            if final_table[qname][1] == "none":
                                #extraction
                                messagefunc(str(c1)+" EXTRACTING: contig "+final_table[qname][0][t][0]+", query "+qname, cols, debugfile)
                                s1 = get_sequence(final_table[qname][0][t][1], seq, extractiontype, flanks, trans_out)
                                print >> debugfile, "- EXTRACTING: final seq", s1[:10], "ranges", final_table[qname][0][t][1]
                                if noq:
                                    if seq.id not in target_set:
                                        seqwritefunc(s1, qname,target_db_name, seq.id, noq, output_dir)
                                        target_set.add(seq.id)
                                else:
                                    seqwritefunc(s1, qname,target_db_name, seq.id, noq, output_dir)
                            else:
                                s1 = get_sequence(final_table[qname][0][t][1], seq, extractiontype, flanks, trans_out)
                                final_table[qname][1][final_table[qname][1].index(final_table[qname][0][t][0])] = s1
                                messagefunc(str(c1)+" BUCKET: contig "+final_table[qname][0][t][0]+", query "+qname, cols, debugfile)
                                dump_bucket = True
                                for buck1 in final_table[qname][1]:
                                    if type(buck1) is str:
                                        dump_bucket = False
                                        break
                                if dump_bucket:
                                    s1 = dumper(final_table[qname][1], extractiontype)
                                    messagefunc(str(c1)+" EXTRACTING: bucket "+qname+" dumped", cols, debugfile)
                                    print >> debugfile, "- EXTRACTING: final seq", s1[:10]#, "ranges", final_table[qname][1]
                                    seqwritefunc(s1, qname,target_db_name, "Merged_"+qname, noq,output_dir)
                            #cleanup
                            del final_table[qname][0][t]
                            break # breaking from target table
                #clean up after all qs are done
                del final_target_table[seq.id]
                c1 = c1 - 1
            if len(final_target_table) == 0:
                messagefunc(str(c1)+" search finished", cols, debugfile, False)
                break
    debugfile.close()

print len(warninglist)
for w in warninglist:
    print w

print >> debugfile_generic, "done"
debugfile_generic.close()

print "done"