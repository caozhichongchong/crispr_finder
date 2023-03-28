import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
# make sure you have model folder in input_script
input_fastq = '/scratch/users/anniz44/Metagenomes/BN10_MG/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/MG_phage/'
input_script = '/scratch/users/anniz44/scripts/1MG/crispr/MG_phage/'
IMG_db = '/scratch/users/anniz44/scripts/database/virus/IMG_VR_2020-10-12_5.1/IMGVR_all_nucleotides.subset.fna'

#fq_format = '_1.fq'
#fq_format2 = '_2.fq'
fq_format = '_1.fasta'
fq_format2 = '_2.fasta'

try:
    os.mkdir(input_script)
except IOError:
    pass

try:
    os.mkdir(output_folder)
except IOError:
    pass

allfastq = glob.glob('%s/*%s'%(input_fastq,fq_format))
for fastq_1 in allfastq:
    fastq_2 = fastq_1.replace(fq_format,fq_format2)
    filename = os.path.split(fastq_1)[-1].split(fq_format)[0]
    output_file = '%s/%s/'%(output_folder,filename)
    cmd = ''
    try:
        os.mkdir(output_file)
    except IOError:
        pass
    #cmd = '#export HDF5_USE_FILE_LOCKING=\'FALSE\'\n'
    if False:
        try:
            f2  = open('%s/%s.fasta'%(output_file,filename),'r')
        except IOError:
            # megahit
            #cmd += 'megahit -1 %s -2 %s -o %s --min-contig-len 5000 -t 40 -m 0.99\n'%(fastq_1,fastq_2,output_file)
            cmd += 'megahit -r %s %s -o %s --out-prefix %s --min-contig-len 5000 -t 40 -m 0.99\n' % (fastq_1,fastq_2, output_file,filename)
            cmd += 'mv %s/megahit_out/final.contigs.fa %s/%s.fasta\n'%(input_script,output_file,filename)
            cmd += 'rm -r %s/megahit_out/\n' %(input_script)
    # phage seq from MG
    #cmd += 'cd %s/\n'%(input_script)
    #cmd += 'python3 /scratch/users/anniz44/bin/pro/MARVEL/marvel_bins.py -i %s -t 40\n'%(output_file)
    fasta = '%s/%s.fasta'%(output_folder,filename)
    cmd += 'source activate virsorter\n'
    cmd += '/scratch/users/anniz44/bin/pro/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl --diamond -f %s --db 2 --wdir %s/%s_phage/ --ncpu 40 --data-dir /scratch/users/anniz44/bin/pro/VirSorter/virsorter-data/\n'%(fasta,output_folder,filename)
    cmd += 'blastn -db %s -query %s -out %s.IMG.txt -outfmt 6 -max_target_seqs 10 -num_threads 40 -qcov_hsp_perc 90 -perc_identity 90\n' % (
        IMG_db,fasta , fasta)
    #cmd += 'predict-metagenome %s/final.contigs.fa\n'%(output_file)
    f1 = open(os.path.join(input_script, '%s.sh'%(filename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
    f1.write('export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n')
    f1.write(cmd)
    f1.close()

# all scripts
if True:
    f1 = open(os.path.join(input_script, '../allMG_phage.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script, '*.sh')):
        try:
            f2 = open('%s.err'%(sub_scripts),'r')
        except IOError:
            f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

    f1.close()
    print('please run %s/%s' % (input_script, '../allMG_phage.sh'))
