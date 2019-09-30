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
optional.add_argument('-lr', dest='local_rec', choices=['none','actual','range'], help='local reciprocator setting', default='range')
optional.add_argument('--is', dest='interstich', action='store_true', help='perform intercontig stiching', default=False)
optional.add_argument('--translate', dest='trans_out', action='store_true', help='translate output (for -et s or -et a)', default=False)
optional.add_argument('--hit-ovlp', metavar='N', help='allowed hit overlap on query, in bp',dest="hit_ovlp", type=int, default=5)
optional.add_argument('--ctg-ovlp', metavar='N', help='allowed contig overlap on query, in bp',dest="ctg_ovlp", type=int, default=1)
optional.add_argument('--recip-ovlp', metavar='N', help='contig overlap on query for reciprocator selection, in bp',dest="recip_ovlp", type=int, default=10)
optional.add_argument('-bt', choices=['blast','hmmer'], help='alignment table type',dest="bt", default="blast")
# blastn, tblastx, blastp = 1 to 1; tblastn = aa in query, dna in db; blastx = dna in query, aa in db. hmmer is always 1 to 1.
optional.add_argument('-ac', choices=['normal','tblastn', 'blastx'], help='alignment coordinate type',dest="ac", default="normal")
optional.add_argument('-r', metavar='file/folder', help='reciprocal search output file or folder',dest="rec_search")
optional.add_argument('-R', metavar='file', help='target locus to reference contig correspondence file',dest="target_ref_file")
optional.add_argument('-m', choices=['eval','bit','ident'], help='metric to use to select best matches',dest="metric", default="eval")
optional.add_argument('--rescale-metric', dest='metricR', action='store_true', help='divide metric value by length of hit region', default=False)

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
    interstich = vars(args)["interstich"]
    # allow_ovlp = vars(args)["allow_ovlp"]
    hit_ovlp = vars(args)["hit_ovlp"]
    ctg_ovlp = vars(args)["ctg_ovlp"]
    recip_ovlp = vars(args)["recip_ovlp"]
    flanks = vars(args)["flanks"]
    output_dir = vars(args)["output"]
    bt = vars(args)["bt"]
    ac = vars(args)["ac"]
    metric = vars(args)["metric"]
    metricR = vars(args)["metricR"]

    print blastfilearg, evalue, filefolder, extractiontype, contignum, local_rec, interstich, flanks, dry_run, noq

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

def target_processor(inpdict, local_rec, metric, metricR, hit_overlap, recip_overlap, cols, debugfile):
    outdict = {}
    for targetkey, targetval in inpdict[1].items():
        #run local actual reciprocal check (check each hit of target only matches one query)
        if local_rec == "actual":
            targetval = actual_reciprocator(targetval, metric, metricR)
        #run hit overlapper and sticher 
        #(cluster hits by direction and overlap on query target dif)
        #(hit stich by cluster)
        tgt_proc_out = hit_sticher(targetval, metric, metricR, hit_overlap, cols, debugfile) #stiched subcontigs per query
        # sys.exit()
        #run local range reciprocal check
        if local_rec == "range":
            tgt_proc_out = range_reciprocator(tgt_proc_out, metric, metricR, recip_overlap, cols, debugfile) #check that each subcontig matches only 1 Q
        outdict[targetkey] = tgt_proc_out
    return outdict

#only run for blastx, tblastn, tblastx; should not matter for blastn?
#input is dictionary {linecounter: [target_f, target_r, target_b, query_f, query_r, query_b, float(row[10]), float(row[11]),float(row[2])]}
#return list of dictionaries
def strand_selector(inpdict, metric, metricR):
    #split into things per direction
    clusterF = {}
    clusterR = {}
    for hitkey, hitval in inpdict.items():
        if hitval[2] == hitval[5]:
            clusterF[hitkey] = hitval
        else:
            clusterR[hitkey] = hitval
    return [clusterF, clusterR]
    
#work on later
def actual_reciprocator(inpdict, metric, metricR):
    for targetkey, targetval in inpdict[1].items():
        for querykey, queryval in targetval.items():
            inpdict[0][querykey]

def hit_sticher(inpdict, metric, metricR, hit_overlap, cols, debugfile):
    returnlist = {} #all queries for the target go here
    for querykey, queryval in inpdict.items():
        clusters = strand_selector(queryval, metric, metricR)
        directions = [True, False]
        stiched_subcontigs = [] #stiched subcontigs for a query go here
        for cluster_index in range(2):
            direct = directions[cluster_index]
            #add conditions for no hits and just 1 hit
            if len(clusters[cluster_index]) > 1:
                messagefunc("running hit overlapper...", cols, debugfile)
                #this is a hitlist {hit index: [hit val], 0: [2@@@22]}
                hitdict = clusters[cluster_index]
                messagefunc("best direction: "+str(direct), cols, debugfile)
            
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
                                                    if refhitkey not in bucketdict: #otherwise it must have already gotten the bucket number
                                                        bucketdict[refhitkey] = nbuck
                                                        nbuck += 1
                                                    if hitkey not in bucketdict: #same here
                                                        bucketdict[hitkey] = nbuck
                                                        nbuck += 1
                                                    #keep both, asign different buckets
                                                #this is this case for subcontig splitter
                                                elif qovlp > hit_overlap:
                                                    #keep both but assign same bucket
                                                    ovlp = False
                                                    if refhitkey in bucketdict: #if one present, the other must have not been checked before
                                                        bucketdict[hitkey] = bucketdict[refhitkey] #put in ref bucket
                                                    elif hitkey in bucketdict:
                                                        bucketdict[refhitkey] = bucketdict[hitkey]
                                                    #no need to up the current bucket since only old buckets were worked with
                                            else:
                                                #keep only one, remove the other by removing from hitdict and bucketdict
                                                if tovlp == qovlp:
                                                    #print "both same overlap (compare with expansion of q and target), overlap, set bool to true"
                                                    # print "check frames of two hits, if in frame, then keep, else remove the worst"
                                                    startshift = abs(min(hitval[0], hitval[1])-min(hitdict[refhitkey][0], hitdict[refhitkey][1]))-1
                                                    #double check, frame only works for some options!
                                                    if startshift > 0 and startshift%3 == 0:
                                                        # print "in frame, keep"
                                                        hitdict[refhitkey] = [min(hitdict[refhitkey][0], hitval[0],hitdict[refhitkey][1], hitval[1]), max(hitdict[refhitkey][0], hitval[0],hitdict[refhitkey][1], hitval[1]), direct, min(hitdict[refhitkey][3], hitval[3],hitdict[refhitkey][4], hitval[4]), max(hitdict[refhitkey][3], hitval[3],hitdict[refhitkey][4], hitval[4]),(not direct),min(hitdict[refhitkey][6], hitval[6]),max(hitdict[refhitkey][7], hitval[7]),max(hitdict[refhitkey][8], hitval[8])]
                                                        if refhitkey not in bucketdict: #otherwise it must have already gotten the bucket number
                                                            bucketdict[refhitkey] = nbuck
                                                            nbuck += 1
                                                        del hitdict[hitkey]
                                                        if hitkey in bucketdict:
                                                            del bucketdict[hitkey]
                                                    else:
                                                        # print "not in frame, pick best"
                                                        hitdict, bucketdict, nbuck = remove_worst(hitdict, bucketdict, nbuck, refhitkey, hitkey, debugfile)
                                                else:
                                                    # print "target ovlp is larger than query, remove the worst, set bool to true"
                                                    hitdict, bucketdict, nbuck = remove_worst(hitdict, bucketdict, nbuck, refhitkey, hitkey, debugfile)
                                                ovlp = True
                                                break
                        if not ovlp:
                            messagefunc("no more overlaps", cols, debugfile)
                #this part will split hits into different subcontig clusters and
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
                            #add distance assessment later too
                            # mindist = 10000000 #very lame
                            if metric == 'eval':
                                # bscore = 1000
                                # for htnm in range(len(tempbuck)):
                                for htnm in tempbuck.keys():
                                    if bhit == -1 or tempbuck[htnm][6] < bscore:
                                        bhit = htnm
                                        bscore = tempbuck[htnm][6]
                            elif metric == 'bit':
                                # bscore = 0
                                # for htnm in range(len(tempbuck)):
                                for htnm in tempbuck.keys():
                                    if bhit == -1 or tempbuck[htnm][7] > bscore:
                                        bhit = htnm
                                        bscore = tempbuck[htnm][7]
                            elif metric == 'ident':
                                # bscore = 0
                                # for htnm in range(len(tempbuck)):
                                for htnm in tempbuck.keys():
                                    if bhit == -1 or tempbuck[htnm][8] > bscore:
                                        bhit = htnm
                                        bscore = tempbuck[htnm][8]
                            sbctg.append(tempbuck[bhit]) #add best to subcontig
                            del bucketdict[bhit] #remove the best from further assessments
                            del hitdict[bhit] #same
                        elif len(tempbuck) == 1:
                            bhit = tempbuck.keys()[0]
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
                for sbctg in subcontigs:
                    #RUN STICHER and return margins and also sequence of regions and gaps
                    messagefunc("running hit sticher...", cols, debugfile)
                    #for a hit
                    median_coords = {}
                    start_coords = {}
                    end_coords = {}
                    start_target = {}
                    end_target = {}
                    eval_hit = {}
                    bit_hit = {}
                    ident_hit = {}
                    for chunk in range(len(sbctg)):
                        median_coords[chunk] = median([sbctg[chunk][3],sbctg[chunk][4]])
                        start_coords[chunk] = min(sbctg[chunk][3],sbctg[chunk][4])
                        end_coords[chunk] = max(sbctg[chunk][3],sbctg[chunk][4])
                        start_target[chunk] = min(sbctg[chunk][0],sbctg[chunk][1])
                        end_target[chunk] = max(sbctg[chunk][0],sbctg[chunk][1])
                        eval_hit[chunk] = sbctg[chunk][6]
                        if metricR:
                            bit_hit[chunk] = sbctg[chunk][7]/abs(sbctg[chunk][0]-sbctg[chunk][1])
                            ident_hit[chunk] = sbctg[chunk][8]/abs(sbctg[chunk][0]-sbctg[chunk][1])
                        else:
                            bit_hit[chunk] = sbctg[chunk][7]
                            ident_hit[chunk] = sbctg[chunk][8]
                        #populate subcontig eval = min eval
                        eval_sbctg = min(eval_hit.values())
                        #populate subcontig bit = max bit
                        bit_sbctg = max(bit_hit.values())
                        #populate subcontig ident = max ident
                        ident_sbctg = max(ident_hit.values())
                        #ranges
                        start_target_subctg = min(start_target.values())
                        end_target_subctg = max(end_target.values())
                        start_query_subctg = min(start_coords.values())
                        end_query_subctg = max(end_coords.values())
                    gapstart = 0
                    # stiched_subcontigs = [[direct],[scores],[ranges],[hits: [range, score], gap, [range, score], gap ...]]
                    stiched_subcontig = []
                    stiched_subcontig.append([direct]) #add direction
                    stiched_subcontig.append([eval_sbctg, bit_sbctg, ident_sbctg]) #add scores
                    stiched_subcontig.append([start_target_subctg, end_target_subctg, start_query_subctg, end_query_subctg]) #add ranges
                    #add per hit information
                    stiched_hits = []
                    for key in sorted(median_coords, key=lambda x: median_coords[x]):
                        if gapstart > 0:
                            stiched_hits.append(start_coords[key]-gapstart)
                        stiched_hits.append([start_target[key], end_target[key], start_coords[key], end_coords[key], eval_hit[key], bit_hit[key], ident_hit[key]])
                        gapstart = end_coords[key]
                    stiched_subcontig.append(stiched_hits)
                    stiched_subcontigs.append(stiched_subcontig)
            
            if len(clusters[cluster_index]) == 1:
                #for a hit
                stiched_subcontigs = []
                median_coords = {}
                start_coords = {}
                end_coords = {}
                start_target = {}
                end_target = {}
                eval_hit = {}
                bit_hit = {}
                ident_hit = {}
                sbctg = clusters[cluster_index].values()
                chunk = 0 #since only one hit
                start_coords[chunk] = min(sbctg[chunk][3],sbctg[chunk][4])
                end_coords[chunk] = max(sbctg[chunk][3],sbctg[chunk][4])
                start_target[chunk] = min(sbctg[chunk][0],sbctg[chunk][1])
                end_target[chunk] = max(sbctg[chunk][0],sbctg[chunk][1])
                eval_hit[chunk] = sbctg[chunk][6]
                if metricR:
                    bit_hit[chunk] = sbctg[chunk][7]/abs(sbctg[chunk][0]-sbctg[chunk][1])
                    ident_hit[chunk] = sbctg[chunk][8]/abs(sbctg[chunk][0]-sbctg[chunk][1])
                else:
                    bit_hit[chunk] = sbctg[chunk][7]
                    ident_hit[chunk] = sbctg[chunk][8]
                stiched_subcontig = []
                stiched_subcontig.append([direct]) #add direction
                stiched_subcontig.append([eval_hit[chunk], bit_hit[chunk], ident_hit[chunk]]) #add scores
                stiched_subcontig.append([start_target[chunk], end_target[chunk], start_coords[chunk], end_coords[chunk]]) #add ranges
                stiched_subcontig.append([[start_target[chunk], end_target[chunk], start_coords[chunk], end_coords[chunk], eval_hit[chunk], bit_hit[chunk], ident_hit[chunk]]]) #add the single hit
                stiched_subcontigs.append(stiched_subcontig)

        for subcont_index in range(len(stiched_subcontigs)):
            returnlist[querykey+"@$"+str(subcont_index)] = stiched_subcontigs[subcont_index]
    return returnlist


def range_reciprocator(inpdict, metric, metricR, recip_overlap, cols, debugfile):
    returnlist = {} #all queries for the target go here
    querylist = inpdict.keys() #list of all queries
    for querykey, queryval in inpdict.items(): #queries for a given target
        cond = True
        ref_query_range = inpdict[querykey][2]
        ref_query_scores = inpdict[querykey][1]
        for key in querylist: #running loop over other queries
            if key != querykey: #all other queries
                current_query_range = inpdict[key][2]
                current_query_scores = inpdict[key][1]
                if getOverlap([current_query_range[0],current_query_range[1]],[ref_query_range[0],ref_query_range[1]]) > recip_overlap:
                    if metric == "eval":
                        #first check eval
                        if current_query_scores[0] < ref_query_scores[0]:
                            cond = False
                            messagefunc(" has better eval hit to "+key, cols, debugfile)
                            break
                        #then bit
                        elif current_query_scores[0] == ref_query_scores[0]:
                            if current_query_scores[1] > ref_query_scores[1]:
                                cond = False
                                messagefunc(" has better ref_query_scores[1] hit to "+key, cols, debugfile)
                                break
                            #then keep but warn
                            elif current_query_scores[1] == ref_query_scores[1]:
                                wrn = "warning, target "+" at query "+querykey+" has equal hits to query "+key+", saved for both!"
                                warninglist.append(wrn)
                                messagefunc(wrn, cols, debugfile)
                                messagefunc("match to current query: ref_query_scores[0] "+str(ref_query_scores[0])+", bitmax "+str(ref_query_scores[1]), cols, debugfile)
                                messagefunc("match to "+key+": "+",".join(map(str, current_query_scores)), cols, debugfile)
                    elif metric == "bit":
                        #first check bit
                        if current_query_scores[1] > ref_query_scores[1]:
                            cond = False
                            messagefunc(" has better ref_query_scores[1] hit to "+key, cols, debugfile)
                            break
                        #then check eval
                        elif current_query_scores[1] == ref_query_scores[1]:
                            if current_query_scores[0] < ref_query_scores[0]:
                                cond = False
                                messagefunc(" has better eval hit to "+key, cols, debugfile)
                                break
                            #then keep but warn
                            elif current_query_scores[0] == ref_query_scores[0]:
                                wrn = "warning, target "+" at query "+querykey+" has equal hits to query "+key+", saved for both!"
                                warninglist.append(wrn)
                                messagefunc(wrn, cols, debugfile)
                                messagefunc("match to current query: ref_query_scores[0] "+str(ref_query_scores[0])+", bitmax "+str(ref_query_scores[1]), cols, debugfile)
                                messagefunc("match to "+key+": "+",".join(map(str, current_query_scores)), cols, debugfile)
                    elif metric == "ident":
                        #first check eval
                        if current_query_scores[2] > ref_query_scores[2]:
                            cond = False
                            messagefunc(" has better ident hit to "+key, cols, debugfile)
                            break
                        #then bit
                        elif current_query_scores[2] == ref_query_scores[2]:
                            if current_query_scores[1] > ref_query_scores[1]:
                                cond = False
                                messagefunc(" has better ref_query_scores[1] hit to "+key, cols, debugfile)
                                break
                            #then keep but warn
                            elif current_query_scores[1] == ref_query_scores[1]:
                                wrn = "warning, target "+" at query "+querykey+" has equal hits to query "+key+", saved for both!"
                                warninglist.append(wrn)
                                messagefunc(wrn, cols, debugfile)
                                messagefunc("match to current query: ident "+str(ref_query_scores[2])+", bitmax "+str(ref_query_scores[1]), cols, debugfile)
                                messagefunc("match to "+key+": "+",".join(map(str, current_query_scores)), cols, debugfile)    
        if cond: #only return good items
            returnlist[querykey] = queryval
    return returnlist
    
def query_processor(inpdict, rec_dict, target_ref, metric, metricR, contignum, contig_overlap, interstich, cols, debugfile):
    returndict = {}
    # 1 reformat target table as query table
    # 2 run reference reciprocation: for each query each target must match back to query transcript from reference
    query_dict = reformat_dict(inpdict, rec_dict, target_ref, metric, metricR, cols, debugfile) #do steps 1 and 2
    # 3 for each query rank targets by scores
    # 4 run sticher, i.e. fill in size of query with contigs in best order, set aside, continue
    for querykey, queryval in query_dict.items():
        if len(queryval) > 1:
            targetlist = rank_targets(queryval, metric, metricR) # step 3
            #need to do ranking again??
            stiched_targets = contig_sticher(targetlist, metric, metricR, contig_overlap, interstich, cols, debugfile) # step 4
        else:
            reformatted_queryval = [queryval.keys()[0], queryval.values()[0]]
            stiched_targets = [[reformatted_queryval[1][1],[reformatted_queryval]]]
        # 5 get top [contignum] contigs from step 4, output
        if len(stiched_targets) > contignum and contignum > 0:
            stiched_targets = stiched_targets[:contignum]
        returndict[querykey] = stiched_targets
    return returndict

# returnlist[querykey+"_"+str(subcont_index)] = stiched_subcontigs[subcont_index]
def reformat_dict(inpdict, rec_dict, target_ref, metric, metricR, cols, debugfile):
    returnlist = {}
    for targetkey, targetval in inpdict.items():
        # reciprocal search section
        for querykey, queryval in targetval.items():
            querykeynew, queryindex = querykey.split("@$")
            if rec_dict != None:
                if reference_reciprocator(querykey, rec_dict, target_ref, metric, cols, debugfile):
                    if querykeynew in returnlist:
                        returnlist[querykeynew][targetkey+"_"+queryindex] = queryval
                    else:
                        returnlist[querykeynew] = {targetkey+"_"+queryindex: queryval}
            else:
                if querykeynew in returnlist:
                    returnlist[querykeynew][targetkey+"_"+queryindex] = queryval
                else:
                    returnlist[querykeynew] = {targetkey+"_"+queryindex: queryval}
    return returnlist


def reference_reciprocator(query, rec_dict, target_ref, metric, metricR, cols, debugfile):
    returnlist = {}
    cond = True
    # check out refence matches:
    if query in target_ref[0]: #checking best contig for Q in reference
        messagefunc("reference check", cols, debugfile)
        best_ref_name = []
        best_ref_val = None
        for target_ref_name, target_ref_val in target_ref[0][query].items():
            # messagefunc(target_ref_name, cols, debugfile)
            ref_ranks_temp = compute_ranks(target_ref_val, metricR)
            # messagefunc(" ".join([str(x) for x in ref_ranks_temp]), cols, debugfile)
            if len(best_ref_name) == 0 and best_ref_val == None:
                best_ref_name = [target_ref_name]
                best_ref_val = ref_ranks_temp
            else:
                if metric == "eval":
                    if ref_ranks_temp[0] < best_ref_val[0]:
                        best_ref_name = [target_ref_name]
                        best_ref_val = ref_ranks_temp
                    elif ref_ranks_temp[0] == best_ref_val[0]:
                        best_ref_name.append(target_ref_name)
                        best_ref_val = ref_ranks_temp
                elif metric == "bit":
                    if ref_ranks_temp[1] > best_ref_val[1]:
                        best_ref_name = [target_ref_name]
                        best_ref_val = ref_ranks_temp
                    elif ref_ranks_temp[1] == best_ref_val[1]:
                        best_ref_name.append(target_ref_name)
                        best_ref_val = ref_ranks_temp
                elif metric == "ident":
                    if ref_ranks_temp[2] > best_ref_val[2]:
                        best_ref_name = [target_ref_name]
                        best_ref_val = ref_ranks_temp
                    elif ref_ranks_temp[2] == best_ref_val[2]:
                        best_ref_name.append(target_ref_name)
                        best_ref_val = ref_ranks_temp
        msg = "best ref: "+",".join(best_ref_name)+"; "+" ".join([str(x) for x in best_ref_val])
        messagefunc(msg, cols, debugfile)
        # check out sample to reference matches
        if target in rec_dict[0]: #checking for best contig for Q in current sample
            messagefunc("sample check", cols, debugfile)
            best_rec_name = None
            best_rec_val = None
            for rec_target, rec_hits in rec_dict[0][target].items():
                # messagefunc(rec_target, cols, debugfile)
                rec_ranks_temp = compute_ranks(rec_hits, metricR)
                if best_rec_name == None and best_rec_val == None:
                    best_rec_name = rec_target
                    best_rec_val = rec_ranks_temp
                else:
                    if metric == "eval":
                        if rec_ranks_temp[0] < best_rec_val[0]:
                            best_rec_name = rec_target
                            best_rec_val = rec_ranks_temp
                    elif metric == "bit":
                        if rec_ranks_temp[1] > best_rec_val[1]:
                            best_rec_name = rec_target
                            best_rec_val = rec_ranks_temp
                    elif metric == "ident":
                        if rec_ranks_temp[2] > best_rec_val[2]:
                            best_rec_name = rec_target
                            best_rec_val = rec_ranks_temp
            msg = "best target: "+best_rec_name+", "+" ".join([str(x) for x in best_rec_val])
            messagefunc(msg, cols, debugfile)
            if best_rec_name not in best_ref_name:
                cond = False
                messagefunc("reciprocator: target "+target+" removed from query "+query+": reciprocal condition violated", cols, debugfile)
        else:
            msg = "no matches to ref in target [strange, possibly queries for forward and reciprocal search differ]"
            messagefunc(msg, cols, debugfile)
            sys.exit()
    else:
        msg = "no matches to query in ref [should not happen, possibly queries for forward and reciprocal search differ]"
        messagefunc(msg, cols, debugfile)
        sys.exit()
    return cond

# {target:[[direction],[scores],[ranges],[[hit range and score], gap, etc]]}
# {target:[ 0[direction], 1[0eval, 1bit, 2ident], 2[ranges], 3[[hit range and score], gap, etc]]}
def rank_targets(inpdict, metric, metricR):
    returnlist = []
    if metric == "eval":
        metrickey = 0
    elif metric == "bit":
        metrickey = 1
    if metric == "ident":
        metrickey = 2
    #iterate through sorted target dict
    if metrickey == 0:
        for targetkey in sorted(inpdict, key=lambda x: inpdict[x][1][metrickey],reverse=True):
            returnlist.append([targetkey, inpdict[targetkey]])
    else:
        for targetkey in sorted(inpdict, key=lambda x: inpdict[x][1][metrickey], reverse=False):
            returnlist.append([targetkey, inpdict[targetkey]])
    return returnlist
    # [target, [ 0[direction], 1[0eval, 1bit, 2ident], 2[ranges], 3[[hit range and score], gap, etc]]]


def contig_sticher(inplist, metric, metricR, contig_overlap, interstich, cols, debugfile):
    messagefunc("running contig sticher...", cols, debugfile)
    if metric == "eval":
        metrickey = 0
    elif metric == "bit":
        metrickey = 1
    if metric == "ident":
        metrickey = 2
    ctg_counter = 0
    returndict = {}
    returnlist = []
    removed_contigs = []
    while len(removed_contigs) < len(inplist):
        current_list = []
        evals = []
        bits = []
        idents = []
        for contig1 in inplist:
            if contig1 not in removed_contigs:
                if interstich:
                    if len(current_list) == 0:
                        current_list.append(contig1)
                        removed_contigs.append(contig1) #record processed target
                        evals.append(contig1[1][1][0])
                        bits.append(contig1[1][1][1])
                        idents.append(contig1[1][1][2])
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
                            evals.append(contig1[1][1][0])
                            bits.append(contig1[1][1][1])
                            idents.append(contig1[1][1][2])
                else:
                    messagefunc("overlapping disabled", cols, debugfile)
                    current_list.append(contig1)
                    removed_contigs.append(contig1) #record processed target
                    evals.append(contig1[1][1][0])
                    bits.append(contig1[1][1][1])
                    idents.append(contig1[1][1][2])
                    break
            
        #these are merged contig score calculation: improve!
        scores = [min(evals), sum(bits), max(idents)]
        
        messagefunc("number of contigs: "+str(len(current_list)), cols, debugfile)
        if len(current_list) > 1:
            #repeat the function from before?
            median_coords = {}
            start_coords = {}
            end_coords = {}
            target_dict = {}
            for target in current_list:
                messagefunc(target[0], cols, debugfile)
                median_coords[target[0]] = median([min(target[1][2][2],target[1][2][3]),max(target[1][2][2],target[1][2][3])])
                start_coords[target[0]] = min(target[1][2][2],target[1][2][3])
                end_coords[target[0]] = max(target[1][2][2],target[1][2][3])
                target_dict[target[0]] = target[1]
            gapstart = 0
            output = []
            for key in sorted(median_coords, key=lambda x: median_coords[x]):
                if gapstart > 0:
                    output.append(start_coords[key]-gapstart)
                # else:
                #     output.append(0)
                output.append([key,target_dict[key]])
                gapstart = end_coords[key]
            returndict[ctg_counter] = [scores, output]
            ctg_counter += 1
        else:
            returndict[ctg_counter] = [scores, [current_list[0]]]
            ctg_counter += 1
    if metrickey == 0:
        for returnkey in sorted(returndict, key=lambda x: returndict[x][0][metrickey],reverse=True):
            returnlist.append(returndict[returnkey])
    else:
        for returnkey in sorted(returndict, key=lambda x: returndict[x][0][metrickey],reverse=False):
            returnlist.append(returndict[returnkey])
    return returnlist
    
def reformat_table(inpdict):
    returnlist = {}
    for querykey, queryval in inpdict.items():
        # total[ contig[ score[ ], contigs[ contig [ name, info[ ]]]]]
        #        N         N-0       N-1    N-1-M   N-1-M-0
        for elemN in queryval: #elemN = N
            for elemM in elemN[1]:
                if type(elemM) is not int:
                    targetname = "".join(elemM[0].split("_")[:-1])
                    if targetname in returnlist:
                        returnlist[targetname].add(querykey)
                    else:
                        returnlist[targetname] = set([querykey])
    return returnlist


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

def remove_worst(hitdict1, bucketdict1, nbuck1, refhitkey, hitkey, debugfile):
    if hitdict1[refhitkey][6] < hitdict1[hitkey][6] or hitdict1[refhitkey][7] > hitdict1[hitkey][7] or abs(hitdict1[refhitkey][0]-hitdict1[refhitkey][1]) >= abs(hitdict1[hitkey][0]-hitdict1[hitkey][1]):
        print >> debugfile, hitdict1[hitkey], "deleted"
        if refhitkey not in bucketdict1: #otherwise it must have already gotten the bucket number
            bucketdict1[refhitkey] = nbuck1
            nbuck1 += 1
        del hitdict1[hitkey]
        if hitkey in bucketdict1:
            del bucketdict1[hitkey]
    else:
        print >> debugfile, hitdict1[refhitkey], "deleted"
        if hitkey not in bucketdict1: #otherwise it must have already gotten the bucket number
            bucketdict1[hitkey] = nbuck1
            nbuck1 += 1
        del hitdict1[refhitkey]
        if refhitkey in bucketdict1:
            del bucketdict1[refhitkey]
    return hitdict1, bucketdict1, nbuck1

def median(lst): #taken from https://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0


def get_sequence(inplist, seq, extractiontype, fls, trans_out1):
    finalseq = Seq("")
    seqlen = len(seq.seq)
    direct = inplist[0][0]
    scores = inplist[1]
    ranges = inplist[2]
    hits = inplist[3]
    if extractiontype == "a":
        for i in hits:
            if trans_out1:
                #translation
                if type(i) is not int:
                    start = min(i[0],i[1])-1
                    end = max(i[0],i[1])
                    # messagefunc("len to translate: "+str(len(seq.seq[start:end])), debugfile)
                    if direct:
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
                    if direct:
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
        for i in hits:
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
        #THIS IS INCORRECT
        start = min(hits[0][0],hits[0][1])-1
        end = max(hits[0][0],hits[0][1])
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
            if direct:
                finalseq += seq.seq[start:end].translate()
            else:
                finalseq += seq.seq[start:end].reverse_complement().translate()
        else:
            if direct:
                finalseq += seq.seq[start:end]
            else:
                finalseq += seq.seq[start:end].reverse_complement()
    if not direct and extractiontype != "a" and extractiontype != "s":
        finalseq = finalseq.reverse_complement()
    return finalseq

def dumper(inplist, extractiontype):
    finalseq = Seq("")
    for i in inplist:
        if type(i) is not int:
            finalseq += i[1]
        else:
            if i > 0 and extractiontype != "n":
                finalseq += Seq("N"*i)
            else:
                finalseq += Seq("N")
    return finalseq
#---------------------------------------------------------------------#
############################ main script ##############################
debugfile_generic = open("absx.log", "w")

#### investigate if this works in batch!
#get terminal window size
if not os.popen('stty size', 'r').read():
    cols = 100
else:
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
    if rec_search != None:
        rec_list = glob.glob(rec_search+"/*.blast")
elif filefolder == "S":
    blastlist = [blastfilearg]
    if not dry_run:
        translist = [targetf]
    if rec_search != None:
        rec_list = [rec_search]

if dry_run:
    messagefunc("dry run, no target files", cols, debugfile_generic, False)
else:
    messagefunc("list of target fasta files detected (mask *.fasta):", cols, debugfile_generic, False)
    for l in translist:
        messagefunc(l, cols, debugfile_generic)

if rec_search != None:
    target_ref = readblastfilefunc(target_ref_file, cols, debugfile_generic)
else:
    target_ref = None

###perhaps some of these are no longer needed
#debug vars
number = 0
numberset = set()
totalloci = 0
#parsing blast files
messagefunc("parsing blast files...", cols, debugfile_generic, False)

b1 = 0
for b in blastlist:
    b1 += 1 #counter
    messagefunc("target "+str(b1)+" out of "+str(len(blastlist)), cols, debugfile_generic, False)
    #set up sample debug file
    debugfile = open(b.split("/")[-1]+"_absx.log", "w")
    messagefunc("target log started: "+b, cols, debugfile_generic, False)
    #read alignment table
    if bt == "blast":
        output = readblastfilefunc(b, cols, debugfile) #output 0 is query, 1 is target
    else:
        output = readhmmerfilefunc(b, cols, debugfile)
    #read reciprocal alignment table, writing only as a debug really, commented out for now
    if rec_search != None:
        if rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast" in rec_list:
            rec_out = readblastfilefunc(rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast", cols, debugfile)
            # with open(rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast_qtable.tab", "w") as rec_qout:
            #     bltableout(rec_out[0],rec_qout, "query")
        elif len(rec_list) == 1:
            rec_out = readblastfilefunc(rec_list[0], cols, debugfile)
            # with open(rec_list[0]+"_qtable.tab", "w") as rec_qout:
            #     bltableout(rec_out[0], rec_qout, "query")
        else:
            print "problem with reciprocal search file"
            print rec_search.rstrip("/")+"/"+b.split("/")[-1]+"_reciprocal.blast", rec_list
            sys.exit()
    else:
        rec_out = None
    final_table = {}
    final_target_table = {}
    #run target processor
    target_table = target_processor(output, local_rec, metric, metricR, hit_ovlp, recip_ovlp, cols, debugfile)
    
    #run query processor
    final_table = query_processor(target_table, rec_out, target_ref, metric, metricR, contignum, ctg_ovlp, interstich, cols, debugfile)

    final_target_table = reformat_table(final_table)

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
        messagefunc("scanning the target fasta file...", cols, debugfile_generic, False)
        messagefunc("--------scanning the target---------", cols, debugfile)
        
        #get the transcriptome filename, matching blast filename
        if filefolder == "M":
            seqname = b[:-6].split("/")[-1]
            target_db_match = [target_db_m1 for target_db_m1 in translist if seqname in target_db_m1]
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
            seqname = translist[0]
            inputf = SeqIO.parse(seqname, "fasta")
            target_db_name = seqname.split("/")[-1]
        
        c1 = len(final_target_table)
        messagefunc("searching for contigs in: "+target_db_name+", total number of contigs: "+str(c1), cols, debugfile, False)
        if noq:
            target_set = set()
        for seq in inputf: #going over seqs in target file
            if seq.id in final_target_table: #looking up same seq in target file
                for qname in final_target_table[seq.id]: #checking it's queries
                    for sprcontig in range(len(final_table[qname])):
                        for t in range(len(final_table[qname][sprcontig][1])): #looking for target in the query table
                            if type(final_table[qname][sprcontig][1][t]) is not int:
                                targetname = "".join(final_table[qname][sprcontig][1][t][0].split("_")[:-1])
                                if targetname == seq.id: #found target in the query table
                                    if len(final_table[qname][sprcontig][1]) == 1:
                                        #extraction
                                        messagefunc(str(c1)+" EXTRACTING: contig "+final_table[qname][sprcontig][1][t][0]+", query "+qname, cols, debugfile)
                                        s1 = get_sequence(final_table[qname][sprcontig][1][t][1], seq, extractiontype, flanks, trans_out)
                                        print >> debugfile, "- EXTRACTING: final seq", s1[:10], "ranges", final_table[qname][sprcontig][1][t][1]
                                        if noq:
                                            if seq.id not in target_set:
                                                seqwritefunc(s1, qname,target_db_name, seq.id, noq, output_dir)
                                                target_set.add(seq.id)
                                        else:
                                            seqwritefunc(s1, qname,target_db_name, seq.id, noq, output_dir)
                                    else:
                                        s1 = get_sequence(final_table[qname][sprcontig][1][t][1], seq, extractiontype, flanks, trans_out)
                                        final_table[qname][sprcontig][1][t][1] = s1
                                        messagefunc(str(c1)+" BUCKET: contig "+final_table[qname][sprcontig][1][t][0]+", query "+qname, cols, debugfile)
                                        dump_bucket = True
                                        for elem1 in final_table[qname][sprcontig][1]:
                                            if type(elem1) is not int:
                                                if type(elem1[1]) is list:
                                                    dump_bucket = False
                                                    break
                                        if dump_bucket:
                                            # sys.exit()
                                            s1 = dumper(final_table[qname][sprcontig][1], extractiontype)
                                            messagefunc(str(c1)+" EXTRACTING: bucket "+qname+" dumped", cols, debugfile)
                                            print >> debugfile, "- EXTRACTING: final seq", s1[:10]#, "ranges", final_table[qname][sprcontig][1]
                                            seqwritefunc(s1, qname,target_db_name, "Merged_"+qname, noq,output_dir)
                                            # cleanup
                                            del final_table[qname][sprcontig][1][t]
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