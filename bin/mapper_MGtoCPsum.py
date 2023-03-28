import glob
import os
from Bio import SeqIO

outputfolder = '/scratch/users/anniz44/Metagenomes/BN10_MG/crispr/mapperresult/'
ref_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly'

def vcf_to_count(mappervcf):
    countmap = dict()
    for lines in open(mappervcf, 'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            depth = int(float(lines_set[4]))
            countmap.setdefault(depth,0)
            countmap[depth] += 1
    return countmap

def countmap_out(countmap, outputfile):
    allout = ['Depth\tNo.positions\n']
    total_position = 0
    for depth in countmap:
        num_position = countmap[depth]
        allout.append('%s\t%s\n'%(depth,num_position))
        total_position += num_position
    f1 = open(outputfile, 'w')
    f1.write(''.join(allout))
    f1.close()
    return total_position

def load_ref(reffasta):
    total_length = 0
    for record in SeqIO.parse(reffasta, 'fasta'):
        total_length += len(str(record.seq))
    return total_length


allmapperfile = glob.glob('%s/*.vcf'%(outputfolder))
allsum = ['sample\tspecies\tref_cluster\ttotal_mappingpositions\ttotal_refpositions\n']
for mapperfile in allmapperfile:
    countmap = vcf_to_count(mapperfile)
    samplefile = os.path.split(mapperfile)[-1].split('.vcf')[0]
    samplename, ref =samplefile.split('_1.fasta_')
    reference_genomefull = '%s/%s/%s.all.spades2.fasta'%(ref_folder, ref, ref)
    total_position = countmap_out(countmap, os.path.join(outputfolder,samplefile + '_sum.txt'))
    try:
        allsum.append('%s\t%s\t%s\t%s\t%s\n'%(samplename,ref.split('_')[0],ref,total_position,load_ref(reference_genomefull)))
    except FileNotFoundError:
        pass

f1 = open('%s/allmappingsum.txt'%(outputfolder), 'w')
f1.write(''.join(allsum))
f1.close()
