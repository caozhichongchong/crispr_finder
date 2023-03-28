import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json
from copy import deepcopy
input_folder = '/scratch/users/anniz44/genomes/crispr_MG/CRISPRCasFinder/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
filename_format = '_final'
crispr_level_cutoff = 2  # only level 2-4 + with cas
spacer_target = '/scratch/users/anniz44/genomes/crispr_MG/summary/BL.widespread.spacer.sharephage.short.txt'
species_target = 'BL'

def load_target():
    alldonors = set()
    allspacers = set()
    for lines in open(spacer_target,'r'):
        lines_set = lines.split('\t')
        alldonors.add('%s/%s_%s*'%(input_folder,lines_set[5],species_target))
        allspacers.add(lines_set[3])
    alldonorslist = []
    for donor in alldonors:
        alldonorslist += glob.glob(donor)
    return [alldonorslist,allspacers]

def parse_crispr(results_crisprcasfinder,allspacers,allstructure,donor_species,donor_speciesshort_structure):
    crispr_ID = 0
    # only crisprcasfinder
    with open(results_crisprcasfinder) as f:
        data = json.load(f)
    crispr_set = dict()
    for subseq in data['Sequences']:
        for crispr in subseq['Crisprs']:
            if crispr['Evidence_Level'] >= crispr_level_cutoff and subseq['Cas'] != []:
                crispr_ID += 1
                spacer_ID = 0
                for spacer in crispr['Regions']:
                    if spacer['Type'] == 'Spacer':
                        record_seq = spacer['Sequence']
                        # complementary reverse + sort
                        record_seq2 = str(Seq(record_seq).reverse_complement())
                        record_list = [record_seq, record_seq2]
                        record_list.sort()
                        record_seq = record_list[0]
                        spacer_ID += 1
                        if record_seq in allspacers:
                            # spacers on the same crispr
                            allstructure.append('%s\tCRISPR_%s\tSpacer_%s\t%s\n'%(donor_species,crispr_ID,spacer_ID,record_seq))
                            crispr_set.setdefault(crispr_ID,[])
                            crispr_set[crispr_ID].append(record_seq)
    for crispr_ID in crispr_set:
        spacer_order = crispr_set[crispr_ID]
        spacer_order2 = deepcopy(spacer_order)
        spacer_order2.reverse()
        if spacer_order2 in allspacer_order:
            spacer_order = spacer_order2
        if spacer_order not in donor_speciesshort_structure:
            donor_speciesshort_structure.append(spacer_order)
            allspacer_order.append(spacer_order)
    return [allstructure,donor_speciesshort_structure]

allspacer_order = []
alldonorslist,allspacers = load_target()
allstructure = []
allstructure.append('donor_species\tCRISPR_ID\tSpacer_ID\tSpacer\n')
alldonor_speciesshort_structure = dict()
for results_crisprcasfinder in alldonorslist:
    donor_species = '_'.join(os.path.split(results_crisprcasfinder)[-1].split('_')[:3])
    donor_speciesshort = '_'.join(os.path.split(results_crisprcasfinder)[-1].split('_')[:2])
    donor_speciesshort_structure = alldonor_speciesshort_structure.get(donor_speciesshort,[])
    allstructure,donor_speciesshort_structure = parse_crispr(results_crisprcasfinder, allspacers, allstructure, donor_species,donor_speciesshort_structure)
    alldonor_speciesshort_structure[donor_speciesshort] = donor_speciesshort_structure

allspacers=list(allspacers)
allspacers.sort()
alldonor_speciesshort_structuresum = []
alldonor_speciesshort_structuresum2 = ['donor_species\t%s\n'%('\t'.join([i[0:8] for i in allspacers]))]
for donor_speciesshort in alldonor_speciesshort_structure:
    donor_speciesshort_structure = alldonor_speciesshort_structure[donor_speciesshort]
    for spacer_order in donor_speciesshort_structure:
        alldonor_speciesshort_structuresum.append('%s\t%s\n'%(donor_speciesshort,'\t'.join(spacer_order)))
        spacer_order_temp = ['']*len(allspacers)
        for i in range(0,len(spacer_order)):
            spacer_order_temp[allspacers.index(spacer_order[i])]=str(i)
        alldonor_speciesshort_structuresum2.append('%s\t%s\n'%(donor_speciesshort,'\t'.join(spacer_order_temp)))

f1 = open('%s/all.BL.spacer.structure.txt' % (output_folder), 'w')
f1.write(''.join(allstructure))
f1.close()

f1 = open('%s/all.BL.spacer.structure.sum.txt' % (output_folder), 'w')
f1.write(''.join(alldonor_speciesshort_structuresum))
f1.close()

f1 = open('%s/all.BL.spacer.structure.summat.txt' % (output_folder), 'w')
f1.write(''.join(alldonor_speciesshort_structuresum2))
f1.close()