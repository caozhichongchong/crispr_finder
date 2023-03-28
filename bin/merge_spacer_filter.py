import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json
input_folder = '/scratch/users/anniz44/genomes/crispr_MG/CRISPRCasFinder/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
#genome_name_meta = '/scratch/users/anniz44/genomes/donor_species/selected_species/file.change.name.new.txt.addtime.addclonal.txt'
genome_name_meta_folder = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/*.genome.cluster.txt')
filename_format = '_final'
crispr_level_cutoff = 2  # only level 2-4 + with cas
extractfasta = True
try:
    os.mkdir(output_folder)
except IOError:
    pass

try:
    os.mkdir(output_folder + '/DR')
except IOError:
    pass

def add_name_crispr(filename,results_crisprcasfinder,newID, allseq_donorspecies,alloutput_donorspecies,alloutput_donorspeciesDR,outputlist, seq_dict):
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
                if extractfasta:
                    DRseq = crispr['DR_Consensus']
                    if DRseq not in allDR:
                        allDR.add(DRseq)
                        alloutput_donorspeciesDR.append('>%s__%s\n%s\n'%(filename,newID,DRseq))
                    for spacer in crispr['Regions']:
                        if spacer['Type'] == 'Spacer':
                            record_seq = spacer['Sequence']
                            # complementary reverse + sort
                            record_seq2 = str(Seq(record_seq).reverse_complement())
                            record_list = [record_seq, record_seq2]
                            record_list.sort()
                            record_seq = record_list[0]
                            outputlist.append('>%s__%s\n%s\n' % (filename, oldID, record_seq))
                            if record_seq not in allseq_donorspecies:
                                allseq_donorspecies.add(record_seq)
                                alloutput_donorspecies.append('>%s__%s\n%s\n' % (this_donor_species, newID, record_seq))
                                newID += 1
                            seq_dict.setdefault(record_seq, set())
                            seq_dict[record_seq].add(filename)
                            oldID += 1
    return [newID, allseq_donorspecies,alloutput_donorspecies,alloutput_donorspeciesDR,outputlist,No_crispr]

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

def load_clonal(genome_name_meta):
    clonal = dict()
    clonal_new = []
    for lines in open(genome_name_meta,'r'):
        if not lines.startswith('oldname'):
            lines_set = lines.split('\n')[0].split('\t')
            newname = lines_set[4].split('_')
            newname2 = '%s_%s_%s_%s'%(newname[0], newname[1],'_'.join(newname[3:]),newname[2])
            CL = lines_set[3].replace('_cluster',' ')
            clonal.setdefault(lines_set[1],CL)
            clonal_new.append('%s\t%s'%(newname2,lines))
    output_fasta(clonal_new, genome_name_meta + '.newname.txt')
    return clonal

def load_clonal_clonal_population(genome_name_meta,clonal):
    for lines in open(genome_name_meta,'r'):
        if not lines.startswith('species'):
            lines_set = lines.split('\n')[0].split('\t')
            species, genomename, cluster, nothing = lines_set[:4]
            CL = '%s %s'%(species,cluster.replace('cluster',''))
            clonal.setdefault(genomename, CL)
    return clonal

# all spacer file
allspacerfiles=glob.glob('%s/*.json'%(input_folder))
# add clonal name
clonal = dict()
for genome_name_meta in genome_name_meta_folder:
    clonal = load_clonal_clonal_population(genome_name_meta,clonal)
# spacer merge
alloutput = []
alloutputcrispr = ['donor_species\tgenome\tCL\tNo.crispr\n']
outputfilename = '%s/all.spacers.fasta'%(output_folder)
donor_species_files = dict()
for files in allspacerfiles:
    filefolder, filename = os.path.split(files)
    filename = filename.split(filename_format)[0]
    this_donor_species = filename.split('_')
    this_donor_species = '%s_%s' % (this_donor_species[0], this_donor_species[1])
    donor_species_files.setdefault(this_donor_species,set())
    donor_species_files[this_donor_species].add(files)

for this_donor_species in donor_species_files:
    filename_loci = dict()
    i = 0
    allfilename = []
    seq_dict = dict()
    alloutput_donorspecies = []
    alloutput_donorspeciesDR = []
    allseq_donorspecies = set()
    newID = 0
    for files in donor_species_files[this_donor_species]:
        filefolder, filename = os.path.split(files)
        filename = filename.split(filename_format)[0]
        newID, allseq_donorspecies,alloutput_donorspecies,alloutput_donorspeciesDR,alloutput,No_crispr = add_name_crispr(
            filename,files, newID,allseq_donorspecies,alloutput_donorspecies,alloutput_donorspeciesDR,alloutput,seq_dict)
        filename_loci.setdefault(filename,i)
        allfilename.append(filename)
        alloutputcrispr.append('%s\t%s\t%s\t%s\n' % (this_donor_species, filename, clonal.get(filename,'None'),No_crispr))
        i += 1
    alloutput_donor = []
    alloutput_donor.append('spacer_seq\t'+'\t'.join(allfilename)+'\n')
    for record_seq in seq_dict:
        result_temp = ['0']*len(allfilename)
        for filename in seq_dict[record_seq]:
            result_temp[filename_loci[filename]]='1'
        alloutput_donor.append('%s\t%s\n'%(record_seq,'\t'.join(result_temp)))
    output_fasta(alloutput_donor, '%s/%s.spacers.summary'%(output_folder,this_donor_species))
    output_fasta(alloutput_donorspecies, '%s/%s.spacers.fasta' % (output_folder, this_donor_species))
    output_fasta(alloutput_donorspeciesDR, '%s/DR/%s.DR.fasta' % (output_folder, this_donor_species))

#output_fasta(alloutput,outputfilename)
f1 = open('%s/allcrisprsize.txt'%(output_folder), 'w')
f1.write(''.join(alloutputcrispr))
f1.close()
