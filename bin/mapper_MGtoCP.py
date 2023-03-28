import glob
import os

input_script = '/scratch/users/anniz44/scripts/1MG/crispr/mapper/'
latest_mapper = glob.glob('/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/mapper-1*.jar')[0]
ref_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly'
metadata = '/scratch/users/anniz44/genomes/crispr_MG/summary/all.spacerchange.sum.all.txt'
metagenome_folder = '/scratch/users/anniz44/Metagenomes/BN10_MG/crispr/'
outputfolder = '/scratch/users/anniz44/Metagenomes/BN10_MG/crispr/mapperresult/'

try:
    os.mkdir(input_script)
except IOError:
    pass

try:
    os.mkdir(outputfolder)
except IOError:
    pass

donor_set = ['aa', 'af', 'am', 'bq', 'cx']
max_penalty = 0.1
def runmapper(database, files, files2, tempbamoutput):
    cmds = '/usr/bin/time -v java -Xms10g -Xmx10g -jar %s --max-penalty %s --num-threads 40 --reference %s --queries %s  --queries %s --out-vcf %s.vcf\n' % (
        latest_mapper, max_penalty, database, files, files2, tempbamoutput)
    return cmds

def readmetadata(metadata):
    species_select = dict()
    for lines in open(metadata,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        donor = lines_set[0]
        if donor in donor_set:
            lineage = lines_set[2]
            species_select.setdefault(donor,[])
            species_select[donor].append(lineage.replace(' ','_clustercluster'))
    print(species_select)
    for donor in species_select:
        reference_genomeset = species_select[donor]
        for reference_genome in reference_genomeset:
            reference_genomefull = '%s/%s/%s.all.spades2.fasta'%(ref_folder, reference_genome, reference_genome)
            metagenome_1set = glob.glob('%s/%s*_1.fasta'%(metagenome_folder, donor))
            for metagenome_1 in metagenome_1set:
                samplename = os.path.split(metagenome_1)[-1]
                metagenome_2 = '%s/%s' % (metagenome_folder, samplename.replace('_1.fasta','_2.fasta'))
                outputvcf = '%s/%s_%s'%(outputfolder,samplename,reference_genome)
                print(donor,reference_genomefull, metagenome_1, metagenome_2, outputvcf)
                cmds = runmapper(reference_genomefull, metagenome_1, metagenome_2, outputvcf)
                f1 = open(os.path.join(input_script, '%s_%s.sh' % (samplename,reference_genome)), 'w')
                f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                f1.close()
    f1 = open(os.path.join(input_script, '../allmapper.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script, '*.sh')):
        f1.write('jobmit %s %s small\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()

readmetadata(metadata)