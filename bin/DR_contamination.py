import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json

output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/DR/'
contamination_cutoff = 2 # consider contamination if DR shared by <= 2 genomes in a species
CRISPRCasFinder_folder = '/scratch/users/anniz44/genomes/crispr_MG/CRISPRCasFinder/'
SP_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
crispr_level_cutoff = 2  # only level 2-4 + with cas
filename_format = '_final'

def load_DR(DR,DR_species_count,species):
    for record in SeqIO.parse(DR, 'fasta'):
        record_seq = str(record.seq)
        record_id = str(record.id)
        DR_species_count.setdefault(record_seq,{})
        DR_species_count[record_seq].setdefault(species,[])
        DR_species_count[record_seq][species].append(record_id.split('__')[0])
    return DR_species_count

def sum_DR(DR_species_count):
    species_set_DR = {}
    for DRseq in DR_species_count:
        DRseq_species = DR_species_count[DRseq]
        if len(DRseq_species) > 1:
            print('DR %s shared across species'%(DRseq))
            # shared by at least 2 species
            DRseq_species_count = []
            DRseq_species_dict = {}
            for species in DRseq_species:
                genomeset = DRseq_species[species]
                DRseq_prevalence = len(genomeset)
                DRseq_species_count.append(DRseq_prevalence)
                DRseq_species_dict.setdefault(DRseq_prevalence,[])
                DRseq_species_dict[DRseq_prevalence].append(species)
            DRseq_species_count.sort()
            for DRseq_prevalence in DRseq_species_count:
                if DRseq_prevalence <= contamination_cutoff:
                    # excluding DR_seq in genomesets
                    for species in DRseq_species_dict[DRseq_prevalence]:
                        genomeset = DRseq_species[species]
                        for genome in genomeset:
                            lineage = '_'.join(genome.split('_')[:-1])
                            species_set_DR.setdefault(lineage,[])
                            species_set_DR[lineage].append(DRseq)
                print(DRseq_prevalence,DRseq_species_dict[DRseq_prevalence])
    return species_set_DR

def add_name_crispr(lineage,DRseq_exclude,filename,results_crisprcasfinder,newID, allseq_donorspecies,alloutput_donorspecies,alloutput_donorspeciesDR, seq_dict):
    oldID = 0
    # only crisprcasfinder
    with open(results_crisprcasfinder) as f:
        data = json.load(f)
    allDR = set()
    No_crispr = 0
    for subseq in data['Sequences']:
        for crispr in subseq['Crisprs']:
            if crispr['Evidence_Level'] >= crispr_level_cutoff and subseq['Cas'] != []:
                No_crispr += 1
                DRseq = crispr['DR_Consensus']
                if DRseq not in DRseq_exclude:
                    # not in contaminated DR set
                    if DRseq not in allDR:
                        allDR.add(DRseq)
                        alloutput_donorspeciesDR.append('>%s__%s\n%s\n' % (filename, newID, DRseq))
                    for spacer in crispr['Regions']:
                        if spacer['Type'] == 'Spacer':
                            record_seq = spacer['Sequence']
                            # complementary reverse + sort
                            record_seq2 = str(Seq(record_seq).reverse_complement())
                            record_list = [record_seq, record_seq2]
                            record_list.sort()
                            record_seq = record_list[0]
                            if record_seq not in allseq_donorspecies:
                                allseq_donorspecies.add(record_seq)
                                alloutput_donorspecies.append('>%s__%s\n%s\n' % (lineage, newID, record_seq))
                                newID += 1
                            seq_dict.setdefault(record_seq, set())
                            seq_dict[record_seq].add(filename)
                            oldID += 1
    return [newID, allseq_donorspecies,alloutput_donorspecies,alloutput_donorspeciesDR,No_crispr]

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

def extract_crispr(species_set_DR):
    for lineage in species_set_DR:
        DRseq = species_set_DR[lineage] # contaminated DR set
        print('re-extract SPs for %s excluding %s' % (lineage, DRseq))
        allspacerfiles = glob.glob('%s/%s*.json' % (CRISPRCasFinder_folder,lineage))
        newID = 0
        seq_dict = dict()
        allseq_donorspecies = set()
        alloutput_donorspecies = []
        alloutput_donorspeciesDR = []
        filename_loci = dict()
        allfilename = []
        i = 0
        for files in allspacerfiles:
            filefolder, filename = os.path.split(files)
            filename = filename.split(filename_format)[0]
            newID, allseq_donorspecies, alloutput_donorspecies, alloutput_donorspeciesDR, No_crispr = add_name_crispr(lineage,DRseq,
                filename, files, newID, allseq_donorspecies, alloutput_donorspecies, alloutput_donorspeciesDR,
                 seq_dict)
            filename_loci.setdefault(filename, i)
            allfilename.append(filename)
            i += 1
        alloutput_donor = []
        alloutput_donor.append('spacer_seq\t' + '\t'.join(allfilename) + '\n')
        for record_seq in seq_dict:
            result_temp = ['0'] * len(allfilename)
            for filename in seq_dict[record_seq]:
                result_temp[filename_loci[filename]] = '1'
            alloutput_donor.append('%s\t%s\n' % (record_seq, '\t'.join(result_temp)))
        output_fasta(alloutput_donor, '%s/../%s.spacers.summary'%(output_folder,lineage))
        output_fasta(alloutput_donorspecies, '%s/../%s.spacers.fasta' % (output_folder, lineage))
        output_fasta(alloutput_donorspeciesDR, '%s/../%s.DR.fasta' % (output_folder, lineage))

# load all DR
allDR = glob.glob('%s/*.DR.fasta'%(output_folder))
print(len(allDR))
DR_species_count = dict()
for DR in allDR:
    filename = os.path.split(DR)[-1]
    donor_species = filename.split('.')[0].replace('PaDi','PB')
    species = donor_species.split('_')[1]
    # find same DR across species
    DR_species_count = load_DR(DR,DR_species_count,species)

# find DR from contamination
species_set_DR = sum_DR(DR_species_count)

# re-extract those species
extract_crispr(species_set_DR)