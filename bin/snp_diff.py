import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

input_snp_fasta = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/*.all.parsi.fasta')
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'

alloutput = '%s/all.SNP.pair.sum'%(output_folder)
f1 = open(alloutput,'w')
f1.write('cluster\tgenome1\tgenome2\tSNPs\n')
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
    species = os.path.split(fastafile)[-1].split('_')[0]
    donor = os.path.split(fastafile)[-1].split('.donor.')[1].split('.all.parsi.fasta')[0]
    CP = os.path.split(fastafile)[-1].split('.all.parsi.fasta')[0]
    print('procession %s %s'%(species,donor))
    allseq = dict()
    allgenome = []
    for lines in open(fastafile):
        if not lines.startswith(' ') and not lines.startswith('Srefer'):
            record_id = lines.split('    ')[0]
            record_seq = lines.split('    ')[1].split('\n')[0]
            genome = record_id.split('_')[-1]
            allseq.setdefault(genome,record_seq)
            allgenome.append(genome)
    allresults = []
    for i in range(0,len(allgenome)-1):
        genome1 = allgenome[i]
        seq1 = allseq[genome1]
        if donor.startswith('H') or donor.startswith('D') or donor.startswith('P'):
            newgenome1 = '%s_%s_%s' % (donor, species, genome1)
        else:
            newgenome1 = '%s_%s_g%s'%(donor,species,genome1)
        for j in range(i,len(allgenome)):
            genome2 = allgenome[j]
            seq2 = allseq[genome2]
            if donor.startswith('H') or donor.startswith('D') or donor.startswith('P'):
                newgenome2 = '%s_%s_%s' % (donor,species, genome2)
            else:
                newgenome2 = '%s_%s_g%s' % (donor, species, genome2)
            allresults.append('%s\t%s\t%s\t%s\n'%(CP,newgenome1,newgenome2,SNP_seq(seq1, seq2)))
    output_fasta(allresults, alloutput)
