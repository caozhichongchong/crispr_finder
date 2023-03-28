import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

#IMG_VR annotate online too https://img.jgi.doe.gov/cgi-bin/vr/main.cgi
plasmid_db = '/scratch/users/anniz44/genomes/NCBI_Genome/newplasmid/*.fna'
#MG_phage = '/scratch/users/anniz44/genomes/crispr_MG/MG_phage/*/final.contigs.fa'
phage_dbmap = '/scratch/users/anniz44/scripts/database/virus/phage.ENA.details.txt'
phage_db = '/scratch/users/anniz44/scripts/database/virus/phage.ENA.fasta'
IMG_db = '/scratch/users/anniz44/scripts/database/virus/IMG_VR_2020-10-12_5.1/IMGVR_all_nucleotides.subset.fna'
IMG_dbmap = '/scratch/users/anniz44/scripts/database/virus/IMG_VR_2020-10-12_5.1/IMGVR_all_Host_information.short.tsv'
spacer_fasta = '/scratch/users/anniz44/genomes/crispr_MG/summary/all.spacers.fasta'
selected_species = ['BA','BL','PaDi']
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/annotation/'
input_script = '/scratch/users/anniz44/scripts/1MG/crispr/annotation/'


def usearch(databaselist, fasta):
    output_name = '%s/spacer_select' % (output_folder)
    for database in databaselist:
        database_name = os.path.split(database)[-1]
        cmds = ''
        try:
            f1 = open('%s.nhr' % (database), 'r')
        except IOError:
            cmds += ('makeblastdb -in %s -dbtype nucl\n' % (database))
        cmds += 'blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 10 -num_threads 40 -qcov_hsp_perc 90 -perc_identity 90\n' % (
            database, fasta, '%s.%s.txt' % (output_name, database_name))
        f1 = open(os.path.join(input_script, '%s.sh' % (database_name)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
        f1.close()

# annotate
if False:
    selected_sp = set()
    for record in SeqIO.parse(spacer_fasta, 'fasta'):
        record_id = str(record.id)
        if any(species in record_id for species in selected_species):
            seq = str(record.seq)
            selected_sp.add('>%s\n%s\n'%(seq,seq))

    f1 = open(spacer_fasta +'.selected','w')
    f1.write(''.join(list(selected_sp)))
    f1.close()
    spacer_fasta = spacer_fasta +'.selected'
    try:
        os.mkdir(input_script)
    except IOError:
        pass
    try:
        os.mkdir(output_folder)
    except IOError:
        pass
    usearch(glob.glob(plasmid_db) + glob.glob(phage_db) + glob.glob(IMG_db),spacer_fasta)
    # all scripts
    f1 = open(os.path.join(input_script, '../allannotation.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script, '*.sh')):
        try:
            f2 = open('%s.err'%(sub_scripts),'r')
        except IOError:
            f1.write('jobmit %s\n' % (sub_scripts))

    f1.close()
    print('please run %s/%s' % (input_script, '../allannotation.sh'))

# sum annotation
def load_blast(blastresult,tag):
    for lines in open(blastresult,'r'):
        spacer = lines.split('\t')[0]
        ref = lines.split('\t')[1]
        allspacer.setdefault(spacer,set())
        allspacer[spacer].add('%s\t%s'%(tag,ref))


allplasmidoutput = glob.glob('%s/*plasmid*.txt'%(output_folder))
allphageoutput = glob.glob('%s/*ENA*.txt'%(output_folder))[0]
allIMGoutput = glob.glob('%s/*IMGVR*.txt'%(output_folder))[0]
allspacer = dict()
for plasmidoutput in allplasmidoutput:
    load_blast(plasmidoutput, 'plasmid')
load_blast(allphageoutput, 'phage')
load_blast(allIMGoutput, 'phage')

alloutput = []
alloutput.append('spacer\tannotation\n')
for spacer in allspacer:
    annoresult = list(allspacer[spacer])
    annoresult.sort()
    alloutput.append('%s\t%s\n'%(spacer,';'.join(annoresult)))

f1 = open(spacer_fasta +'.selected.annosum.details.txt','w')
f1.write(''.join(alloutput))
f1.close()


