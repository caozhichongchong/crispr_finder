import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json
input_folder = '/scratch/users/anniz44/genomes/crispr_MG/crass_details/'
input_folder2 = '/scratch/users/anniz44/genomes/crispr_MG/crisprcasfinder_output/Result/crisprcasfinder_output/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
genome_name_meta = '/scratch/users/anniz44/genomes/donor_species/selected_species/file.change.name.new.txt.addtime.addclonal.txt'
filename_format = '_1.fastq'
filename_format2 = '_final'
#seq_limit1 = 27
#seq_limit2 = 43
crispr_level_cutoff = 2  # only level 2-4 + with cas

try:
    os.mkdir(output_folder)
except IOError:
    pass

def qualify_len(seq):
    if len(seq)>=seq_limit1 and len(seq)<=seq_limit2:
        return True

def changespacername(newname,oldgroup,oldID):
    group = int(newname.split('G')[1].split('SP')[0])
    if oldgroup!= 0 and oldID != 0:
        if group == oldgroup:
            # same group
            newID = oldID + 1
        else:
            # diff group
            newID = 1
    else:
        # first spacer
        newID = 1
    newname = 'G%sSP%s__%s' % (group, newID, newname)
    return [newname,group,newID]

def add_name_both(filename,fastafile,outputlist,allbacoutput, seq_dict, Spacers = False):
    oldgroup, oldID = [0,0]
    # both crassfinder and crisprcasfinder
    results_crisprcasfinder = glob.glob('%s/%s%s/result.json'%(input_folder2,filename,filename_format2))
    results_crisprcasfinderset = set()
    if results_crisprcasfinder != []:
        for lines in open(results_crisprcasfinder[0],'r'):
            if lines.startswith('\"Sequence\"'):
                    results_crisprcasfinderset.add(str(lines.split(': ')[1].split('\"')[1]))
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        if Spacers:
            record_id, oldgroup, oldID = changespacername(record_id,oldgroup, oldID)
        record_seq = str(record.seq)
        if Spacers:
            if record_seq not in results_crisprcasfinderset:
                allbacoutput.append('%s\t%s\n'%(filename,record_seq))
            else:
                outputlist.append('>%s__%s\n%s\n' % (filename, record_id, record_seq))
                seq_dict.setdefault(record_seq, set())
                seq_dict[record_seq].add(filename)
    return [outputlist,allbacoutput]

def add_name_crass(filename,fastafile,outputlist,allbacoutput, seq_dict, Spacers = False): #not used
    oldgroup, oldID = [0,0]
    # both crassfinder and crisprcasfinder
    results_crisprcasfinder = glob.glob('%s/%s%s/result.json'%(input_folder2,filename,filename_format2))
    results_crisprcasfinderset = set()
    if results_crisprcasfinder != []:
        with open(results_crisprcasfinder[0]) as f:
            data = json.load(f)
        for subseq in data['Sequences']:
            for crispr in subseq['Crisprs']:
                if crispr['Evidence_Level'] >= crispr_level_cutoff and subseq['Cas'] != []:
                    for spacer in crispr['Regions']:
                        if spacer['Type'] == 'Spacer':
                            record_seq = spacer['Sequence']
                            results_crisprcasfinderset.add(record_seq)
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        if Spacers:
            record_id, oldgroup, oldID = changespacername(record_id,oldgroup, oldID)
        record_seq = str(record.seq)
        if Spacers:
            #if qualify_len(record_seq):
            outputlist.append('>%s__%s\n%s\n' % (filename, record_id, record_seq))
            seq_dict.setdefault(record_seq, set())
            seq_dict[record_seq].add(filename)
    return [outputlist,allbacoutput]

def add_name_crispr(filename,fastafile,newID, allseq_donorspecies,alloutput_donorspecies,outputlist,allbacoutput, seq_dict, Spacers = False):
    oldID = 0
    # only crisprcasfinder
    results_crisprcasfinder = glob.glob('%s/%s%s/result.json'%(input_folder2,filename,filename_format2))
    if results_crisprcasfinder != []:
        with open(results_crisprcasfinder[0]) as f:
            data = json.load(f)
        for subseq in data['Sequences']:
            for crispr in subseq['Crisprs']:
                if crispr['Evidence_Level'] >= crispr_level_cutoff and subseq['Cas'] != []:
                    for spacer in crispr['Regions']:
                        if spacer['Type'] == 'Spacer':
                            record_seq = spacer['Sequence']
                            # complementary reverse + sort
                            record_seq2= str(Seq(record_seq).reverse_complement())
                            record_list = [record_seq,record_seq2]
                            record_list.sort()
                            record_seq = record_list[0]
                            #if qualify_len(record_seq):
                            outputlist.append('>%s__%s\n%s\n' % (filename, oldID, record_seq))
                            if record_seq not in allseq_donorspecies:
                                allseq_donorspecies.add(record_seq)
                                alloutput_donorspecies.append('>%s__%s\n%s\n' % (this_donor_species, newID, record_seq))
                                newID += 1
                            seq_dict.setdefault(record_seq, set())
                            seq_dict[record_seq].add(filename)
                            oldID += 1
    return [newID, allseq_donorspecies,alloutput_donorspecies,outputlist,allbacoutput]

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
            clonal.setdefault(lines_set[1],newname2)
            clonal_new.append('%s\t%s'%(newname2,lines))
    output_fasta(clonal_new, genome_name_meta + '.newname.txt')
    return clonal

allspacerfiles = glob.glob('%s/am_*/spacers.fasta'%(input_folder))
# add clonal name
clonal = load_clonal(genome_name_meta)
# spacer merge
alloutput = []
allbacoutput = []
outputfilename = '%s/all.spacers.fasta'%(output_folder)
donor_species_files = dict()
for files in allspacerfiles:
    filefolder, filename = os.path.split(files)
    filefolder, filename = os.path.split(filefolder)
    filename = filename.replace(filename_format, '')
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
    allseq_donorspecies = set()
    newID = 0
    for files in donor_species_files[this_donor_species]:
        filefolder, filename = os.path.split(files)
        filefolder, filename = os.path.split(filefolder)
        filename = filename.replace(filename_format, '')
        newID, allseq_donorspecies,alloutput_donorspecies,alloutput, allbacoutput = add_name_crispr(
            filename,files, newID,allseq_donorspecies,alloutput_donorspecies,alloutput,allbacoutput,seq_dict, True)
        #alloutput,allbacoutput = add_name_crass(filename,files, alloutput,allbacoutput,seq_dict, True)
        filename_loci.setdefault(filename,i)
        allfilename.append(filename)
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

output_fasta(alloutput,outputfilename)
