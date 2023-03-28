import argparse
import os
import random
import math
import csv
import re
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from Bio import Phylo
from functools import reduce

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="file name of your tree", type=str, default='16S.nwk',metavar='16S.nwk')
parser.add_argument("-sp",
                    help="spacer file", type=str, default='am_BaVu.spacer.summary',metavar='.spacer.summary')

################################################## Definition ########################################################
args = parser.parse_args()
Pgain = 0.6 # prbability of gain, pe
Ploss = 0.4 # prbability of loss
Pcorrect = {'1':Pgain, # number of gains times chance of gains
            '0':Ploss} # number of loss times chance of loss
changegenomename = True
################################################### Function ########################################################
def changename(genomename):
    if len(genomename) > 8:
        newgenomename = 'S' + genomename[0:4] + '_' + genomename[-4:]
    return newgenomename

def read_table(filename):
    annoall = {}
    nodename = {}
    i = 0
    for lines in open(filename, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        if lines.startswith('spacer_seq'):
            for samplename in lines_set[1:]:
                if changegenomename:
                    samplename = changename(samplename)
                nodename.setdefault(samplename,i)
                i += 1
        else:
            spacer = lines_set[0]
            annoall.setdefault(spacer, [])
            for pre_ab in lines_set[1:]:
                annoall[spacer].append(pre_ab)
    return [annoall,nodename]

def assign_internal_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            old_name = clade.name
            clade.name = '%d_%s' % (idx, clade.name)
        else:
            old_name = None
            clade.name = str(idx)
        if old_name:
            names[old_name] = clade.name
    return names

def mismatch(c,chars):
    return len([x for x in chars if c not in x])

def allmatch(c,chars):
    return len([x for x in chars if c == x])

def pars(chars):
    # clade terminal
    if not chars:
        return set()
    # internal clades
    # chars is all scores of clades under the internal clades
    # all_chars is all unique stored anno/scores in chars
    all_chars = set(reduce( lambda a,b: a.union(b), chars ))
    # scores are [a score, number of different scores]
    scores = [(c,mismatch(c,chars)*Pcorrect[c]) for c in all_chars] # need swap Pcorrect # number of gains times chance of gains, number of loss times chance of loss
    # min_score is the minimum number of different scores
    min_score = min(scores,key=lambda x:x[1])[1] # find consensus among scores
    # choose the scores of highest agreement (minimum number of different scores)
    # we set 0.5 as the cutoff
    # return {0.0} for more 0.0 than 1.0
    # return {1.0} for more 1.0 than 0.0
    # return {0.0, 1.0} for half 0.0 and half 1.0 (final result is 0.5)
    return set([x[0] for x in scores if not x[1] > min_score])

def pars_prob(chars):
    # clade terminal
    if not chars:
        return set()
    # internal clades
    # chars is all scores of clades under the internal clades
    # all_chars is all unique stored anno/scores in chars
    all_chars = set(reduce(lambda a, b: a.union(b), chars))
    # scores are [a score, number of same scores]
    scores = [float(c) * allmatch(c, chars) * Pcorrect[c] for c in all_chars]
    return np.sum(scores)/len(all_chars)

def down_pass(clade,data,anno):
    # set the anno/scores from bottom to up
    for x in clade.clades:
        down_pass(x,data,anno)
    # set anno to reference terminal clades
    if clade.name in anno:
        data[clade.name] = set(anno[clade.name])
        return
    # ignore unknown terminal clades
    if clade.is_terminal():
        return
    # infer traits for internal clades that have reference terminal clades below
    # based on all scores of clades under the internal clades
    # pass the scores of highest agreement (minimum number of different scores)
    chars = pars_prob([data[x.name] for x in clade.clades if x.name in data])
    if chars:
        data[clade.name] = chars
    return

def up_pass(parent_chars,clade,data,anno):
    # set the anno/scores from up to bottom
    # reference terminal clades
    if clade.name in anno:
        data[clade.name] = set(anno[clade.name])
        if clade.is_terminal():
            return
    else:
        # correct scores for internal clades using parent scores
        # based on the scores of parent and itself
        if clade.name in data:
            data[clade.name] = pars_prob([data[clade.name],parent_chars])
        # infer traits for internal clades that have no reference terminal clade below
        # infer traits for unknown terminal clades
        # pass the scores of parent
        else:
            data[clade.name] = parent_chars
    for x in clade.clades:
        up_pass(data[clade.name],x,data,anno)
    return


################################################### Programme #######################################################
workingDir, filename = os.path.split(args.t)
treefile = args.t
spacerfile = args.sp

# read annotation files of traits in references (0.0 or 1.0)
# get numericIDs for each name
annoall,nodename = read_table(spacerfile)
tree = Phylo.read(treefile, 'newick')
# assign names to internal nodes
internal_names = assign_internal_names(tree)
# switch annotations to new names
nodename_list = [names for names in nodename]
new_names = dict([ (v,k) for (k,v) in internal_names.items() ]) #789_123->123 #123->abc.fasta
f1=open(spacerfile + '.temp.parsi.txt','w')
for spacer in annoall:
    anno = dict()
    i = 0
    for pre_ab in annoall[spacer]:
        anno.setdefault(nodename_list[i],pre_ab)
        i += 1
    anno = dict([(internal_names[k],v) for (k,v) in anno.items() ])
    # data stores the scores of internal and terminal clades
    data = {}
    # set the scores of internal and reference terminal clades
    down_pass(tree.clade,data,anno)# set the anno/scores from bottom to up
    up_pass(set(),tree.clade,data,anno)# set the anno/scores from up to bottom
    internals = [(x,data[x]) for x in data.keys()]
    # output the mean score of the scores with the highest agreement (0.0, 0.5 or 1.0)
    for k,v in internals:
        if k in new_names:#leaves
            if new_names[k] in nodename and k not in anno:
                f1.write('%s\t%s\t%s\n'%(spacer,nodename[new_names[k]],str(np.mean(list(map(int, list(v)))))))
            elif k not in anno:#reference
                f1.write('%s\t%s\t%s\n' % (spacer, k, str(np.mean(list(map(int, list(v)))))))
        elif k == '0':#root nodes
            f1.write('%s\t%s\t%s\n' % (spacer, 'root', str(np.mean(list(map(int, list(v)))))))
f1.close()

# add sum file to the original file
spacer_new = dict()
newsample = []
for lines in open(spacerfile + '.temp.parsi.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    spacer_new.setdefault(lines_set[0],[])
    if lines_set[1] not in newsample:
        newsample.append(lines_set[1])
    spacer_new[lines_set[0]].append(lines_set[2])
alloutput = []
for lines in open(spacerfile,'r'):
    lines = lines.split('\n')[0]
    lines_set = lines.split('\t')
    if lines.startswith('spacer_seq'):
        lines += '\t%s\n'%('\t'.join(newsample))
    else:
        lines += '\t%s\n'%('\t'.join(spacer_new.get(lines_set[0],[])))
    alloutput.append(lines)
f1=open(spacerfile + '.parsi.txt','w')
f1.write(''.join(alloutput))
f1.close()