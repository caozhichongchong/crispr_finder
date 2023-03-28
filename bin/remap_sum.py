# sum remapping -> remap_sum.py
# mismatch + complement
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
input_folder = '/scratch/users/anniz44/genomes/crispr_MG/remapping/'
depth_cutoff = 1
def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

for mappingfiles in glob.glob('%s/*'%(input_folder)):
    this_donor_species = os.path.split(mappingfiles)[-1].split('.grep.txt')[0]
    seq_dict = dict()
    allfilename = []
    filename_loci = dict()
    i = 0
    for lines in open(mappingfiles,'r'):
        lines_set = lines.split('\n')[0].split(' ')
        record_seq, filename, depth = lines_set[0:3]
        depth = int(depth)
        filename = filename.split('_1.fastq')[0].split('_2.fastq')[0]
        if depth >= depth_cutoff:
            # pick up representativeness
            record_seq2 = str(Seq(record_seq).reverse_complement())
            record_list = [record_seq, record_seq2]
            record_list.sort()
            record_seq = record_list[0]
            seq_dict.setdefault(record_seq, set())
            seq_dict[record_seq].add(filename)
        elif depth > 0:
            print(lines)
        if filename not in allfilename:
            filename_loci.setdefault(filename, i)
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
    print(this_donor_species)
