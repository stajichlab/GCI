#!/usr/bin/env python3

from Bio import Phylo
import sys, os, gzip, csv, io

import itertools

""" Compute distances between nodes. Should go in a library I think
"""
def pairwise_distances(tree):
    taxa = tree.get_terminals()
    pairwise = {}
    for left in taxa:
        pairwise[left] = {}
        for right in taxa:
            if( left == right ):
                pairwise[left][right] = 0
            else:
                pairwise[left][right] = tree.distance(left,right)
    print("done calculating pairwise")
    return pairwise

""" This function computes a weighted gene conservation index.
Should probably be moved to a library
"""

def weighted_GCI(querysp,gene,gene_hits,pair_distances,tree):
    gci = 0.0
    count = 0
    for htaxon in gene_hits[querysp][gene]:
        if querysp not in pair_distances:
            pair_distances[querysp] = {}
        if htaxon not in pair_distances[querysp]:
            pair_distances[querysp][htaxon] = tree.distance(querysp,htaxon)
                
        distance = pair_distances[querysp][htaxon]
        gci += gene_hits[querysp][gene][htaxon][1]        
        count += 1
 #       print("distance is %s for q(%s,%s) to h(%s,%s)" %
 #             (distance,querysp,gene,htaxon, gene_hits[querysp][gene][htaxon][0]))

#        print(query,htaxa,distance,gene_hits[htaxa])
#    print("gci %s for %s %s" % (gci, querysp,gene))

    if gci > 0:
        return [(gci / count),count]
    else:
        return [gci,count]
    
folder   = 'phmmer_out'
treefile = 'species.tre'
if len(sys.argv) > 1:
    folder = sys.argv[1]

if len(sys.argv) > 2:
    treefile = sys.argv[2]

species_tree = Phylo.read(treefile, 'newick')

if not species_tree:
    print("No Tree in file %s" % treefile)
    exit()

#tip_pairwise_distances = pairwise_distances(species_tree)
tip_pairwise_distances = {}
#print(tip_pairwise_distances)

gene_profiles = {}

for (dirpath, dirnames, filenames) in os.walk(folder):
    for filename in filenames:
        with gzip.open(os.path.join(folder,filename) ) as fh:
            domtblio = io.TextIOWrapper(fh, newline="")
            linect = 0
            for line in domtblio:
                if line.startswith("#"):
                    continue
                else:
                    row = line.split(None,23) # whitespace split max 23 columns
                    if(len(row) == 0 ):
                        print("skipping row with no data in file %s at line %d" %
                              (filename,linect))
                        continue
                    genename   = row[0]
                    genelen    = row[2]
                    
                    genenamelst = genename.split('|')
                
                    hitsp      = genenamelst[0]
                    hitid      = genenamelst[1]
                
                    queryname  = row[3]
                    queryvec   = queryname.split('|')
                    qsp        = queryvec[0]
                    qid        = queryvec[1]

                    querylen   = row[5]

                    fullevalue = float(row[6])
                    fullscore  = float(row[7])
                    fullbias   = float(row[8])

                    hit_ct     = row[9]
                    hit_tot    = row[10]

                    if qsp not in gene_profiles:
                        gene_profiles[qsp] = {}

                    if qid not in gene_profiles[qsp]:
                        gene_profiles[qsp][qid] = {}

                    if (hitsp not in gene_profiles[qsp][qid] or
                        gene_profiles[qsp][qid][hitsp][1] < fullscore):

                        gene_profiles[qsp][qid][hitsp] = [hitid,fullscore,fullevalue]

                linect += 1

for sp in gene_profiles:
    with open("%s.tsv" % (sp), 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter="\t",
                               quoting=csv.QUOTE_MINIMAL)
        
        for gene in sorted(gene_profiles[sp]):
            GCI = weighted_GCI(sp,gene,gene_profiles,
                               tip_pairwise_distances,species_tree)
            
            csvwriter.writerow([sp,gene,GCI[0],GCI[1]])

    
