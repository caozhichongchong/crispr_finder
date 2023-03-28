import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

input_folder = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/*')
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
output_script = '/scratch/users/anniz44/scripts/1MG/crispr/remap'

os.system('rm -r %s'%(output_script))
try:
    os.mkdir(output_script)
except IOError:
    pass

for donor_species_folder in input_folder:
    donor_species = os.path.split(donor_species_folder)[-1]
    database = '%s/%s.spacers.fasta' % (output_folder,donor_species)
    outputfile = '%s/../remapping/%s.grep.txt' % (output_folder, donor_species)
    allseq = set()
    i = 0
    cmds = '#!/bin/bash\n'
    try:
        f1 = open(database,'r')
        f1.close()
        for record in SeqIO.parse(database, 'fasta'):
            allseq.add(str(record.seq))
        for fq1 in glob.glob('%s/fastq/*_1.fastq'%(donor_species_folder)):
            fq2 = fq1.replace('_1.fastq','_2.fastq')
            for seq in allseq:
                cmds += '/scratch/users/anniz44/scripts/1MG/crispr/count.sh %s %s >> %s\n' % (seq, fq1, outputfile)
                cmds += '/scratch/users/anniz44/scripts/1MG/crispr/count.sh %s %s >> %s\n' % (seq, fq2, outputfile)
                seq_rc = str(Seq(seq).reverse_complement())
                cmds += '/scratch/users/anniz44/scripts/1MG/crispr/count.sh %s %s >> %s\n' % (seq_rc, fq1, outputfile)
                cmds += '/scratch/users/anniz44/scripts/1MG/crispr/count.sh %s %s >> %s\n' % (seq_rc, fq2, outputfile)
                i += 1
                if (i%1000) == 0:
                    f1 = open(os.path.join(output_script, '%s.%s.remap.sh'%(donor_species,int(i/1000))), 'w')
                    f1.write(cmds)
                    f1.close()
                    cmds = '#!/bin/bash\n'
        f1 = open(os.path.join(output_script, '%s.%s.remap.sh' % (donor_species, int(i / 1000))), 'a')
        f1.write(cmds)
        f1.close()
    except FileNotFoundError:
        print(donor_species_folder)

f1 = open(os.path.join(output_script, '../allremap.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(output_script, '*.sh')):
    f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (output_script, '../allremap.sh'))
