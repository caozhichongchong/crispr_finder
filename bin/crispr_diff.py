import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

crispr_sharing_cutoff = 1 # crispr shows in at least 1 genome
input_snp_fasta = glob.glob('/scratch/users/anniz44/genomes/crispr_MG/summary/*.spacers.summary')
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'

alloutput = '%s/all.crispr.pair.sum'%(output_folder)
f1 = open(alloutput,'w')
f1.write('genome1\tgenome2\tcrispr_diss\n')
f1.close()

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'a')
    f1.write(''.join(outputlist))
    f1.close()

def SNP_seq(seq1, seq2):
    SNP_total = 0
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
    return SNP_total

for fastafile in input_snp_fasta:
    donor_species = os.path.split(fastafile)[-1].split('_')[0]
    print('procession %s'%(donor_species))
    allseq = dict()
    allgenome = []
    for lines in open(fastafile):
        lines_set = lines.split('\n')[0].split('\t')
        if not lines.startswith('spacer_seq'):
            allgenomewithcrispr = sum([int(j) for j in lines_set[1:]])
            print(allgenomewithcrispr,lines_set[1:])
            if allgenomewithcrispr >= crispr_sharing_cutoff:
                for i in range(0,len(allgenome)):
                    allseq[allgenome[i]] += lines_set[i+1]
        else:
            allgenome = lines_set[1:]
            for genome in allgenome:
                allseq.setdefault(genome, '')
    print(allseq)
    allresults = []
    for i in range(0,len(allgenome)-1):
        genome1 = allgenome[i]
        seq1 = allseq[genome1]
        for j in range(i,len(allgenome)):
            genome2 = allgenome[j]
            seq2 = allseq[genome2]
            allresults.append('%s\t%s\t%s\n'%(genome1,genome2,SNP_seq(seq1, seq2)))
    output_fasta(allresults, alloutput)
