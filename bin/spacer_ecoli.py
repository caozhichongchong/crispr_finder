import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json
input_folder = '/scratch/users/anniz44/genomes/crispr_MG/ecoli/CRISPRCasFinder/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/ecoli/summary/'
crispr_level_cutoff = 1  # all spacers

def spacer_check(results_crisprcasfinder):
    # only crisprcasfinder
    ecoli_name = os.path.split(results_crisprcasfinder)[-1].split('.fasta')[0]
    allfilename.append(ecoli_name)
    with open(results_crisprcasfinder) as f:
        data = json.load(f)
    for subseq in data['Sequences']:
        for crispr in subseq['Crisprs']:
            if crispr['Evidence_Level'] >= crispr_level_cutoff:
                for spacer in crispr['Regions']:
                    if spacer['Type'] == 'Spacer':
                        record_seq = spacer['Sequence']
                        # complementary reverse + sort
                        record_seq2 = str(Seq(record_seq).reverse_complement())
                        record_list = [record_seq, record_seq2]
                        record_list.sort()
                        record_seq = record_list[0]
                        # store seq
                        spacer_dict[record_seq] = spacer_dict.get(record_seq, set())
                        spacer_dict[record_seq].add(ecoli_name)

def merge_spacer():
    # merge spacer if one is in another
    allspaders = spacer_dict.keys()
    spacer_length = dict()
    for spacer in allspaders:
        lenspacer = len(spacer)
        spacer_length[lenspacer] = spacer_length.get(lenspacer, set())
        spacer_length[lenspacer].add(spacer)
    alllength = list(spacer_length.keys())
    alllength.sort()
    total_spacer_length = len(alllength)
    spacer_needmerge = []
    for i in range(0,total_spacer_length-1):
        for j in range(i+1,total_spacer_length):
            for spacer1 in spacer_length[alllength[i]]:
                spacer_needmerge += [[spacer1,spacer2] for spacer2 in spacer_length[alllength[j]] if spacer1 in spacer2]
    print(spacer_needmerge)
    remove_spacer = set()
    for spacer_pair in spacer_needmerge:
        spacer1, spacer2 = spacer_pair
        spacer_dict[spacer2].update(spacer_dict[spacer1])
        remove_spacer.add(spacer1)
    print(remove_spacer)
    for spacer in remove_spacer:
        spacer_dict.pop(spacer)

def output_spacer():
    allsummary = []
    allsummary.append('spacer\t%s\n'%('\t'.join(allfilename)))
    for spacer in spacer_dict:
        ecoli_with_spacer = spacer_dict[spacer]
        allsummary.append('%s\t%s\n'%(spacer,'\t'.join([str(i in ecoli_with_spacer) for i in allfilename])))
    f1 = open('%s/ecoli.spacer.sum.txt'%(output_folder),'w')
    f1.write(''.join(allsummary))
    f1.close()


spacer_dict = dict()
allfilename = []
allresult = glob.glob('%s/*.json'%(input_folder))
allresult.sort()
print(allresult)
for resultfile in allresult:
    spacer_check(resultfile)

merge_spacer()
output_spacer()