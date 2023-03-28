import argparse
import os
from Bio import Phylo
from itertools import combinations
from datetime import datetime
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t",
                    help="file name of your tree", type=str, default='16S.nwk',metavar='16S.nwk')
parser.add_argument("-sp",
                    help="spacer file", type=str, default='am_BaVu.spacer.summary',metavar='.spacer.summary')

################################################## Definition ########################################################
args = parser.parse_args()
Pgain = 1 # penalty of gain
Ploss = 1 # penalty of loss
Pnochange = 0 # penalty of no change
changegenomename = False
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

def anno_ref(clade,data,anno):
    # set the anno/scores
    for x in clade.clades:
        anno_ref(x, data, anno)
    # set anno to reference terminal clades
    if clade.name in anno:
        data[clade.name] = set(anno[clade.name])
        return
    # ignore unknown terminal clades
    if clade.is_terminal():
        return
    return

def set_up_tree(clade, data, anno,multiscore_set):
    if clade.name in anno and clade.is_terminal():
        # reference terminal clades
        return
    if clade.name not in anno:
        # internal nodes n scores
        multiscore_set.append(clade.name)
    for x in clade.clades:
        set_up_tree(x, data, anno, multiscore_set)
    return

def exhaust_tree(clade, data, anno, temp_tree,multiscore_set,changeset):
    if clade.name in anno:
        # reference terminal clades, 1 score
        temp_tree[clade.name] = anno[clade.name]
        if clade.is_terminal():
            return
    else:
        # internal nodes n scores
        # pick up a score, change one per tree
        if clade.name not in multiscore_set:
            temp_tree[clade.name] = list(data[clade.name])[0]
        else:
            if clade.name in [multiscore_set[i] for i in changeset]:
                temp_tree[clade.name] = '1'
            else:
                temp_tree[clade.name] = '0'
    for x in clade.clades:
        exhaust_tree(x, data, anno, temp_tree,multiscore_set,changeset)
    return

def gain_loss(parent_status,child_status, score):
    if child_status == parent_status:
        score += Pnochange
    elif parent_status == '1' and child_status == '0':
        score += Ploss
    elif parent_status == '0' and child_status == '1':
        score += Pgain
    return score

def score_tree(parent_status,bestscore,temp_tree,clade,score = 0):
    # compute score for gain and loss
    score = gain_loss(parent_status, temp_tree[clade.name], score)
    if clade.is_terminal():
        return score
    # already a worse score
    if bestscore > -1 and score > bestscore:
        return score
    for x in clade.clades:
        score = score_tree(temp_tree[clade.name], bestscore, temp_tree, x, score)
    return score

def check_max_likelihood(bestscore,besttree,temp_tree,clade):
    score = score_tree(temp_tree[clade.name], bestscore, temp_tree, clade, 0)
    if bestscore == -1:
        bestscore = score
        besttree = temp_tree
    elif score < bestscore:
        bestscore = score
        besttree = temp_tree
    return [besttree,bestscore]

def exhaust_combination(multiscore_set):
    perm_all = []
    for i in range(1,len(multiscore_set)+1):
        perm_all += [i for i in combinations(range(0, len(multiscore_set)), i)]
    return perm_all

def max_likelihood(clade, data, anno):
    # set up all nodes with > 1 candidate scores and empty score set
    # from root to bottom
    multiscore_set = []
    set_up_tree(clade, data, anno,multiscore_set)
    #print('process %s combinations'%(multiscore_set))
    # store best tree
    besttree = {}
    bestscore = -1
    # exhaust tree with one change at a time
    perm_all = exhaust_combination(multiscore_set)
    perm_all.append('')
    for changeset in perm_all: # more combination
        temp_tree = {}
        exhaust_tree(clade, data, anno, temp_tree,multiscore_set,changeset)
        besttree, bestscore = check_max_likelihood(bestscore,besttree,temp_tree,clade)
    return [besttree, bestscore]


################################################### Programme #######################################################
workingDir, filename = os.path.split(args.t)
treefile = args.t
spacerfile = args.sp
print('start parsimony reconstruction',datetime.now())
# read annotation files of spacer presence/absence in references (0.0 or 1.0)
# get numericIDs for each name
annoall,nodename = read_table(spacerfile)
tree = Phylo.read(treefile, 'newick')
# assign names to internal nodes
internal_names = assign_internal_names(tree)
# switch annotations to new names
nodename_list = [names for names in nodename]
new_names = dict([ (v,k) for (k,v) in internal_names.items() ]) #789_123->123 #123->abc.fasta
f1=open(spacerfile + '.temp.parsi.txt','w')
spacer_num = 0
for spacer in annoall:
    anno = dict()
    i = 0
    for pre_ab in annoall[spacer]:
        anno.setdefault(nodename_list[i],pre_ab)
        i += 1
    anno = dict([(internal_names[k],v) for (k,v) in anno.items() ])
    # data stores the scores of terminal clades
    data = {}
    # set the scores of reference terminal clades
    anno_ref(tree.clade,data,anno)
    data,panelty_score = max_likelihood(tree.clade,data,anno)
    internals = [(x,data[x]) for x in data.keys()]
    # output the best tree
    for k,v in internals:
        if k in new_names:#leaves
            if new_names[k] in nodename and k not in anno:
                f1.write('%s\t%s\t%s\t%s\n'%(spacer,nodename[new_names[k]],v,panelty_score))
            elif k not in anno:#reference
                f1.write('%s\t%s\t%s\t%s\n' % (spacer, k, v,panelty_score))
        else:#internal nodes
            f1.write('%s\t%s\t%s\t%s\n' % (spacer, k, v,panelty_score))
    spacer_num += 1
    if spacer_num%10 == 0:
        print('processed %s spacers'%(spacer_num))
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
print('finished parsimony reconstruction',datetime.now())