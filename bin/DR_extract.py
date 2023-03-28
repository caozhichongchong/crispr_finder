import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json
input_folder = '/scratch/users/anniz44/genomes/crispr_MG/CRISPRCasFinder/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/DR/'
target = ['D30_BL_19','H09_BL_17','H02_BL_23']
crispr_level_cutoff = 2  # only level 2-4 + with cas
filename_format = '_final'

def order_seq(seqset,filename):
    newoutput = []
    oldID = 0
    for seq in seqset:
        newoutput.append('>%s_%s\n%s\n'%(filename,oldID,seq))
        oldID += 1
    return newoutput

def add_name_crispr(filename,results_crisprcasfinder,allDR,allSP):
    with open(results_crisprcasfinder) as f:
        data = json.load(f)
    allSP_strain = set()
    allDR_strain = set()
    for subseq in data['Sequences']:
        for crispr in subseq['Crisprs']:
            if crispr['Evidence_Level'] >= crispr_level_cutoff and subseq['Cas'] != []:
                for spacer in crispr['Regions']:
                    if spacer['Type'] == 'Spacer':
                        record_seq = spacer['Sequence']
                        allSP_strain.add(record_seq)
                    elif spacer['Type'] == 'DR':
                        record_seq = spacer['Sequence']
                        allDR_strain.add(record_seq)
    allSP += order_seq(allSP_strain,filename)
    allDR += order_seq(allDR_strain, filename)

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

allDR = []
allSP = []
allspacerfiles=glob.glob('%s/*.json'%(input_folder))
for files in allspacerfiles:
    filefolder, filename = os.path.split(files)
    filename = filename.split(filename_format)[0]
    if filename in target:
        add_name_crispr(filename, files,allDR,allSP)

output_fasta(allDR,'%s/allDR.fasta'%(output_folder))
output_fasta(allSP,'%s/allSP.fasta'%(output_folder))