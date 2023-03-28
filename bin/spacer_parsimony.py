import argparse
import os, glob
from Bio import Phylo
from datetime import datetime
import math
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-sp",
                    help="spacer folder", type=str,
                    default='/scratch/users/anniz44/genomes/crispr_MG/summary/',metavar='.')
parser.add_argument("-dmrca",
                    help="dmrca file for each lineage", type=str,
                    default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/alldmrca.txt',metavar='.')
parser.add_argument("-tree",
                    help="tree folder", type=str,
                    default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/tree/',metavar='.')
parser.add_argument("-treefasta",
                    help="tree fasta folder", type=str,
                    default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/final/',metavar='.')
parser.add_argument("-pgain",
                    help="penalty for gain", type=float, default=1.30102999566,
                    metavar='0 for gain penalty evaluation, or a positive number for spacer ancentral reconstruction (e.g. 1.79588)')
parser.add_argument("-pnochange",
                    help="penalty for no changes", type=float, default=0.04575749056,
                    metavar='0 for no change penalty evaluation, or a positive number for spacer ancentral reconstruction (e.g. 1.79588)')

################################################## Definition ########################################################
args = parser.parse_args()
if args.pgain == 0:
    probability_gain = [math.pow(10,x/10.0)/math.pow(10,4) for x in range(21,31)]# not including 1
else:
    probability_gain = [math.pow(10,-1*args.pgain)]
print(probability_gain)
changegenomename = True
empty_case = ('-', 'N')
min_no_strains_for_eachspacer = 2 # at least 2 strains to have the same spacer
################################################### new class #########################################################
__metaclass__ = type

class newnode:
    # create a class to store node
    def init(self,node):
        self.name = node
        self.parent = []
        self.children = []
        self.anno = ''
        self.score = {'1': 0,# gain,lose
                    '0': 0}
    def addparent(self,parentnode):
        self.parent.append(parentnode)
    def addchild(self,childnodes):
        self.children+=childnodes
    def addanno(self,anno):
        self.anno = anno
    def compute_allscore(self,treelist):
        if self.children != []:
            for child in self.children:
                childnode = treelist[child]
                if childnode.anno == '1':
                    self.score['0'] += childnode.score['1'] + Pgain
                    self.score['1'] += childnode.score['1'] + Pnochange
                elif childnode.anno == '0':
                    self.score['0'] += childnode.score['0'] + Pnochange
                    self.score['1'] += childnode.score['0'] + Ploss
                else:
                    # anno
                    self.score['0'] += min(childnode.score['1'] + Pgain, childnode.score['0'] + Pnochange)
                    self.score['1'] += min(childnode.score['1'] + Pnochange, childnode.score['0'] + Ploss)
    def compute_bestscore(self,treelist,besttree):
        besttree[self.name] = '1'
        if self.score['0'] < self.score['1']:
            besttree[self.name] = '0'
        if self.children != []:
            for child in self.children:
                childnode = treelist[child]
                if childnode.anno == '1':
                    self.score['0'] += childnode.score['1'] + Pgain
                    self.score['1'] += childnode.score['1'] + Pnochange
                elif childnode.anno == '0':
                    self.score['0'] += childnode.score['0'] + Pnochange
                    self.score['1'] += childnode.score['0'] + Ploss
                else:
                    # anno unknown
                    self.score['0'] += min(childnode.score['1'] + Pgain, childnode.score['0'] + Pnochange)
                    self.score['1'] += min(childnode.score['1'] + Pnochange, childnode.score['0'] + Ploss)

    def getscore(self,anno,temp_tree):
        total_score = 0
        total_nor_score = 0 # normalized agains Pgain + Ploss
        if self.children != []:
            for child in self.children:
                annochild = temp_tree[child]
                total_score += compute_score(anno, annochild)
                total_nor_score += compute_normalized_score(anno, annochild)
        return [total_score, total_nor_score]


################################################### Function ########################################################
def changename(genomename):
    if len(genomename) > 8:
        genomename = 'S' + genomename[0:4] + '_' + genomename[-4:]
    else:
        genomename = 'S' + genomename
    return genomename

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
            no_strains_with_thisspacer = 0
            spacer = lines_set[0]
            if 'NNN' not in spacer:
                # not ambiguous
                annoall.setdefault(spacer, [])
                for pre_ab in lines_set[1:]:
                    annoall[spacer].append(pre_ab)
                    if int(pre_ab) == 1:
                        no_strains_with_thisspacer += 1
                if no_strains_with_thisspacer < min_no_strains_for_eachspacer:
                    # at least X strains to have this spacer
                    annoall.pop(spacer)
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

def anno_tree(clade,anno,multiscore_set,newtree_all):
    # set the anno/scores
    newnode_temp = newtree_all.get(clade.name,newnode())
    if clade.name not in newtree_all:
        newnode_temp.init(clade.name)
    newnode_temp.addchild([x.name for x in clade.clades])
    if clade.name not in anno:
        # not annotated nodes
        multiscore_set.append(clade.name)
    newtree_all.setdefault(clade.name, newnode_temp)
    # set anno to reference terminal nodes
    if clade.name in anno:
        newnode_temp.addanno(anno[clade.name])
        return
    # return at terminal clades
    if not clade.is_terminal():
        for x in clade.clades:
            newnode_temp = newtree_all.get(x.name, newnode())
            if x.name not in newtree_all:
                newnode_temp.init(x.name)
            newnode_temp.addparent(clade.name)
            newtree_all.setdefault(x.name, newnode_temp)
            anno_tree(x, anno,multiscore_set,newtree_all)
    return

def tree_compute_allscore(newtree_all):
    # bottom to up
    allnode = list(newtree_all.keys())
    allnode.reverse()
    for node in allnode:
        newnode_temp = newtree_all.get(node, newnode())
        newnode_temp.compute_allscore(newtree_all)
        newtree_all[node] = newnode_temp
        ##print(node,newnode_temp.parent,newnode_temp.children,newnode_temp.score)
    return newtree_all

def compute_score(parentanno,childanno):
    if parentanno == '1' and childanno == '0':
        return Ploss
    elif parentanno == '0' and childanno == '1':
        return Pgain
    return Pnochange

def compute_normalized_score(parentanno,childanno):
    if parentanno == '1' and childanno == '0':
        return Ploss
    elif parentanno == '0' and childanno == '1':
        return Pgain
    return Pnochange

def compute_allchanges(besttree,anno,new_names):
    total_changeset = [0, 0, 0]  # nochange, gain, loss
    reference = [k for k in besttree if new_names.get(k,'None') == 'Srefer']
    if reference == []:
        reference = '0'
    else:
        reference = reference[0]
    root_anno = besttree[reference]
    for k in besttree:
        if k in anno:
            total_changeset[compute_gainloss(root_anno,besttree[k])] += 1
    return total_changeset

def compute_gainloss(parentanno,childanno):
    if parentanno == '1' and childanno == '0':
        return 2 #loss
    elif parentanno == '0' and childanno == '1':
        return 1 # gain
    return 0

def child_bestscore(besttree,node,newtree_all):
    newnode_temp = newtree_all.get(node, newnode())
    if newnode_temp.children == []:
        return
    for child in newnode_temp.children:
        newnode_temp_child = newtree_all.get(child, newnode())
        if newnode_temp_child.anno == '':
            child_0 = compute_score(besttree[node], '0') + newnode_temp_child.score['0']
            child_1 = compute_score(besttree[node], '1') + newnode_temp_child.score['1']
            besttree[child] = '1'
            if child_0 < child_1:
                besttree[child] = '0'
        else:
            besttree[child] = newnode_temp_child.anno
        child_bestscore(besttree, child, newtree_all)
    return

def score_tree(newtree_all,temp_tree,score = 0):
    normalized_score = 0
    for node in newtree_all:
        newnode_temp = newtree_all.get(node, newnode())
        anno_current = temp_tree[node]
        # compute score for gain and loss
        tempresult = newnode_temp.getscore(anno_current,temp_tree)
        score += tempresult[0]
        normalized_score += tempresult[1]
    return [score, normalized_score]

def tree_compute_bestscore(newtree_all):
    # up to bottom
    besttree = {}
    # start from the root
    node = '0'
    newnode_temp = newtree_all.get(node, newnode())
    bestscore = min(newnode_temp.score['0'], newnode_temp.score['1'])
    besttree[node] = '1'
    # find the best score
    if newnode_temp.score['0'] < newnode_temp.score['1']:
        besttree[node] = '0'
    # find best score for the children
    child_bestscore(besttree, node, newtree_all)
    # double check bestscore
    score, normalized_score = score_tree(newtree_all,besttree,0)
    if '%.3f'%(bestscore) != '%.3f'%(score):
        print('wrong best tree: best score %s tree score %s'%(bestscore,score))
    return [besttree,bestscore,normalized_score]

def max_likelihood(newtree_all, anno):
    # set up all nodes with > 1 candidate scores and empty score set
    # from root to bottom
    multiscore_set = []
    anno_tree(tree.clade, anno,multiscore_set,newtree_all)
    # tree compute all scores
    newtree_all = tree_compute_allscore(newtree_all)
    # tree find the best score
    besttree,bestscore,normalized_score = tree_compute_bestscore(newtree_all)
    return [besttree, bestscore,normalized_score]

def SNP_seq(seq1, seq2):
    SNP_total = 0
    total_length_all = min(len(seq1),len(seq2))
    for i in range(0, total_length_all):
        if seq1[i] not in empty_case and seq2[i] not in empty_case:
            if seq1[i] != seq2[i]:
                # a SNP excluding empty/missing allele
                SNP_total += 1
    return SNP_total

def load_fasta(fastafile,allnodename,tree):
    # tree distance to SNP ratio
    samplename1 = ''
    record_seq1 = ''
    SNP_tree_distance = 0
    tree_dis_min = 0
    for lines in open(fastafile,'r'):
        if not lines.startswith(' '):
            lines_set = lines.split('\n')[0].split('    ')
            samplename = lines_set[0]
            record_seq = lines_set[1]
            if SNP_tree_distance == 0 and samplename in allnodename:
                if samplename1 != '':
                    try:
                        tree_dis = tree.distance(samplename1, samplename)
                        if tree_dis > 0:
                            if tree_dis_min == 0:
                                tree_dis_min = tree_dis
                            tree_dis_min = min(tree_dis_min,tree_dis)
                            SNP_pair = SNP_seq(record_seq, record_seq1)
                            if SNP_pair > 0:
                                SNP_tree_distance = SNP_pair / tree_dis
                                return SNP_tree_distance
                    except ValueError:
                        print('no tree dis for %s %s' % (samplename1, samplename))
                samplename1 = samplename
                record_seq1 = record_seq
            elif samplename in allnodename:
                print('%s not in %s'%(samplename,allnodename))
    # if no SNPs, treat minimum tree dis as 1 SNP
    if SNP_tree_distance == 0 and tree_dis_min > 0:
        SNP_tree_distance = 1/tree_dis_min
    return SNP_tree_distance

def load_dmrca(dmrca_file):
    dmrca_all = dict()
    for lines in open(dmrca_file, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        if not lines.startswith('donor_species'):
            dmrca_all.setdefault(lines_set[0],float(lines_set[2]))
    return dmrca_all

################################################### Programme #######################################################
# load dmrca
dmrca_all = load_dmrca(args.dmrca)

# find all spacer summary
spacerdir = args.sp
treedir = args.tree
fastadir = args.treefasta
#allspacerfile = glob.glob('%s/*.spacers.summary'%(spacerdir))
allspacerfile = glob.glob('%s/H26_BL*.spacers.summary'%(spacerdir))+\
                glob.glob('%s/H26_BA*.spacers.summary'%(spacerdir))+ \
                glob.glob('%s/af_EsCo*.spacers.summary' % (spacerdir)) + \
                glob.glob('%s/D84_BL*.spacers.summary' % (spacerdir)) +\
                glob.glob('%s/am_TuSa*.spacers.summary'%(spacerdir))+\
                glob.glob('%s/cx_BL*.spacers.summary'%(spacerdir))+ \
                glob.glob('%s/aa_FaSp*.spacers.summary' % (spacerdir)) + \
                glob.glob('%s/cx_EsCo*.spacers.summary' % (spacerdir))+ \
                glob.glob('%s/D571_PB*.spacers.summary' % (spacerdir))+ \
                glob.glob('%s/am_BaVu*.spacers.summary' % (spacerdir))
f3 = open('%s/all.spacerchange.sum.txt' % (spacerdir), 'w')
f3.write('Species\tDonor\tCL\tSpacer\tPnochange\tPgain\tNormalized_score\tNo_scores\tNo.gain\tNo.loss\tNo.nochanges\tdMRCA\n')
# process each spacer summary
for spacerfile in allspacerfile:
    donor = os.path.split(spacerfile)[-1].split('_')[0]
    species = os.path.split(spacerfile)[-1].split('.')[0].split('_')[1].replace('PB','PaDi')
    alltreefile = glob.glob('%s/%s*%s*tree'%(treedir,species,donor))
    for treefile in alltreefile:
        try:
            # read tree
            tree = Phylo.read(treefile, 'newick')
            print('process %s'%(treefile))
            # get numericIDs for each name
            # read annotation files of spacer presence/absence in references (0.0 or 1.0)
            annoall,nodename = read_table(spacerfile)
            # load SNP to tree ratio
            SNP_tree_distance = load_fasta(os.path.join(fastadir, os.path.split(treefile)[-1].split('.out.tree')[0]),
                                           nodename, tree)
            print(SNP_tree_distance, 'No. snps')
            if SNP_tree_distance > 0:
                # with SNPs
                CL = treefile.split('_clustercluster')[1].split('.donor')[0]
                CLall = os.path.split(treefile)[-1].split('.all.parsi.fasta.out.tree')[0]
                # compute dmrca
                dmrca = dmrca_all.get(CLall, 0)
                print('start parsimony reconstruction for %s dmrca %s' % (CLall, dmrca), datetime.now())
                # assign names to internal nodes
                internal_names = assign_internal_names(tree)
                # switch annotations to new names
                nodename_list = [names for names in nodename]
                new_names = dict([(v, k) for (k, v) in internal_names.items()])  # 789_123->123 #123->abc.fasta
                f1 = open(spacerfile + '.CL%s.temp.parsi.txt' % (CL), 'w')
                f2 = open(spacerfile + '.CL%s.pgain.txt' % (CL), 'a')
                #f2.write('Spacer\tPnochange\tPgain\tNormalized_score\tNo_scores\n')
                spacer_num = 0
                for spacer in annoall:
                    anno = dict()
                    i = 0
                    for pre_ab in annoall[spacer]:
                        anno.setdefault(nodename_list[i], pre_ab)
                        i += 1
                    anno = dict([(internal_names[k], v) for (k, v) in anno.items() if k in internal_names])
                    # find best penalty
                    minchanges = -1
                    Pgainbest = 0
                    Pnochangebest = 0
                    pgain_output = []
                    for Pgain_pro in probability_gain:
                        Pgain = -math.log(Pgain_pro, 10)
                        if args.pnochange == 0:
                            probability_nochange = [1-math.pow(10, -x) for x in range(1,10) if 1-math.pow(10, -x) + Pgain_pro < 1]
                        else:
                            probability_nochange = [math.pow(10, -1 * args.pnochange)]
                        for Pnochange_pro in probability_nochange:
                            Pnochange = -math.log(Pnochange_pro, 10)
                            Ploss = -math.log((1 - Pnochange_pro - Pgain_pro), 10)
                            newtree_all = dict()
                            besttree, bestscore, normalized_score = max_likelihood(newtree_all, anno)
                            internals = [(x, besttree[x]) for x in besttree.keys()]
                            pgain_output.append('%s\t%s\t%s\t%s\t%.3f\n' % (spacer, Pnochange, Pgain, normalized_score, bestscore))
                            if minchanges == -1:
                                minchanges = normalized_score
                                Pgainbest = Pgain
                                Pnochangebest = Pnochange
                            elif normalized_score < minchanges:
                                minchanges = normalized_score
                                Pgainbest = Pgain
                                Pnochangebest = Pnochange
                    f2.write(''.join(pgain_output))
                    # set the scores of reference terminal clades
                    Pgain = Pgainbest
                    Pnochange = Pnochangebest
                    newtree_all = dict()
                    besttree, bestscore, normalized_score = max_likelihood(newtree_all, anno)
                    internals = [(x, besttree[x]) for x in besttree.keys()]
                    # compute total SNPs by tree distance
                    # SNP = int(SNP * SNP_tree_distance)
                    # compute gain or loss
                    total_changeset = compute_allchanges(besttree, anno, new_names)
                    f3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\n' % (species, donor, CL, spacer,Pnochange,
                                                                                 Pgain, normalized_score, bestscore,
                                                                                 total_changeset[1], total_changeset[2],
                                                                                 total_changeset[0],
                                                                                 dmrca))
                    # output the best tree: spacer, tree nodes, abs/pre, penalty
                    for k, v in internals:
                        if k in new_names:  # leaves
                            if k in anno:
                                f1.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (spacer, new_names[k], v, normalized_score, Pgain,Pnochange))
                            elif k in new_names:
                                f1.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (spacer, new_names[k], v, normalized_score, Pgain, Pnochange))
                            else:  # reference
                                f1.write('%s\tI%s\t%s\t%s\t%s\t%s\n' % (spacer, k, v, normalized_score, Pgain, Pnochange))
                        else:  # internal nodes
                            f1.write('%s\tI%s\t%s\t%s\t%s\t%s\n' % (spacer, k, v, normalized_score, Pgain, Pnochange))
                    spacer_num += 1
                    if spacer_num % 25 == 0:
                        print('processed %s spacers' % (spacer_num))
                f1.close()
                f2.close()
                # add sum file to the original file
                spacer_new = dict()
                newsample = []
                if len(probability_gain) == 1:
                    # found best score
                    for lines in open(spacerfile + '.CL%s.temp.parsi.txt' % (CL), 'r'):
                        lines_set = lines.split('\n')[0].split('\t')
                        spacer_new.setdefault(lines_set[0], [])
                        if lines_set[1] not in newsample:
                            newsample.append(lines_set[1])
                        spacer_new[lines_set[0]].append(lines_set[2])
                    alloutput = []
                    for lines in open(spacerfile, 'r'):
                        lines = lines.split('\n')[0]
                        lines_set = lines.split('\t')
                        if lines.startswith('spacer_seq'):
                            lines += '\t%s\n' % ('\t'.join(newsample))
                            alloutput.append(lines)
                        elif lines_set[0] in spacer_new:
                            lines += '\t%s\n' % ('\t'.join(spacer_new.get(lines_set[0], [])))
                            alloutput.append(lines)
                    f1 = open(spacerfile + '.CL%s.parsi.txt' % (CL), 'w')
                    f1.write(''.join(alloutput))
                    f1.close()
                os.system('rm %s' % (spacerfile + '.CL%s.temp.parsi.txt' % (CL)))
                print('finished parsimony reconstruction', datetime.now())
        except ValueError:
            pass # no tree no SNPs
f3.close()