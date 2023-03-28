import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/DR/'

def load_DR(DR,DR_set,Donorspecies_set,DR_Donorspeciesset):
    for record in SeqIO.parse(DR, 'fasta'):
        record_seq = str(record.seq)
        Donorspecies_set[donor_species].add(record_seq)
        DR_set.setdefault(record_seq,set())
        DR_set[record_seq].add(species)
        DR_Donorspeciesset.setdefault(record_seq, set())
        DR_Donorspeciesset[record_seq].add(donor_species)
    return [DR_set,Donorspecies_set,DR_Donorspeciesset]

allDR = glob.glob('%s/*.DR.fasta'%(output_folder))
print(len(allDR))
DR_set = dict()
DR_Donorspeciesset = dict()
Donorspecies_set = dict()
for DR in allDR:
    filename = os.path.split(DR)[-1]
    donor_species = filename.split('.')[0].replace('PaDi','PB')
    species = donor_species.split('_')[1]
    Donorspecies_set.setdefault(donor_species,set())
    DR_set, Donorspecies_set,DR_Donorspeciesset = load_DR(DR,DR_set,Donorspecies_set,DR_Donorspeciesset)

allsum = ['DR\tNum_species\tSpecies\n']
for record_seq in DR_set:
    specieset = DR_set[record_seq]
    allsum.append('%s\t%s\t%s\n'%(record_seq,len(specieset),';'.join(list(specieset))))

f1 = open('%s/allDRspecies.sum'%(output_folder),'w')
f1.write(''.join(allsum))
f1.close()

allsum = ['DR\tNum_species\tDonor_species\n']
for record_seq in DR_Donorspeciesset:
    specieset = DR_Donorspeciesset[record_seq]
    allsum.append('%s\t%s\t%s\n'%(record_seq,len(specieset),';'.join(list(specieset))))

f1 = open('%s/allDRdonorspecies.sum'%(output_folder),'w')
f1.write(''.join(allsum))
f1.close()

allsum = ['Donor_species\tNum_DR\tDR\n']
for donor_species in Donorspecies_set:
    record_seqset = Donorspecies_set[donor_species]
    allsum.append('%s\t%s\t%s\n'%(donor_species,len(record_seqset),';'.join(list(record_seqset))))

f1 = open('%s/alldonorspeciesDR.sum'%(output_folder),'w')
f1.write(''.join(allsum))
f1.close()
