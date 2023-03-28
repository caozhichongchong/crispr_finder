import glob,os

allfasta = glob.glob('/scratch/users/anniz44/Metagenomes/BN10_MG/assembly/bins/*.fa')
outputdir = '/scratch/users/anniz44/genomes/crispr_MG/CRISPRCasFinderMeta/'
scriptdir = '/scratch/users/anniz44/scripts/1MG/crispr/'
scriptdir_crispr = '%s/CRISPRCasFinderscriptMeta/'%(scriptdir)

try:
    os.mkdir(outputdir)
except IOError:
    pass

os.system('rm -r %s'%(scriptdir_crispr))

try:
    os.mkdir(scriptdir_crispr)
except IOError:
    pass

try:
    os.mkdir(scriptdir_crispr + '/secondshell/')
except IOError:
    pass

try:
    os.mkdir(scriptdir_crispr + '/singularity/')
except IOError:
    pass

try:
    os.mkdir(scriptdir_crispr + '/assembly/')
except IOError:
    pass

try:
    os.mkdir(scriptdir_crispr + '/binning/')
except IOError:
    pass

def runspades(file1,file2,output_name):
    temp_output = output_name.split('.fasta')[0]
    cmds = '%s -1 %s -2 %s -o %s --threads 40 --memory 220 --only-assembler\n' % \
            ('spades.py',file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -r %s\n' % (temp_output)
    return cmds

def runbinning(assembly):
    cmds = 'source activate concoct_env\n'

def runderep(donor):
    cmds = 'py37\nmodule add c3ddb/glibc/2.14\ndRep dereplicate output_directory -g path/to/genomes/%s*.fasta'%(donor)

def runCRISPRCasFinder(fasta, i):
    try:
        f1 = open(fasta.replace('.tar.gz',''),'r')
    except IOError:
        os.system('#tar -xvf %s'%(fasta))
    fasta = fasta.replace('.tar.gz','')
    fastaname = os.path.split(fasta)[-1]
    tempoutput = '%s/%s' % (outputdir, fastaname)
    try:
        f1 = open('%s.result.json'%(tempoutput),'r')
    except IOError:

        cmds = '#!/bin/bash\nsource ~/.bashrc\n'
        cmds += 'export PATH=/scratch/users/anniz44/bin/pro/vmatch-2.3.1-Linux_x86_64-64bit:$PATH\n/usr/bin/perl /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CRISPRCasFinder.pl -so /scratch/users/anniz44/bin/pro/CRISPRCasFinder/sel392v2.so -cf /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CasFinder-2.0.3 -drpt /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G -out %s_folder -in %s\n' % (
            tempoutput, fasta
        )
        cmds += 'mv %s_folder/result.json %s.result.json\n' % (tempoutput, tempoutput)
        cmds += 'rm -r %s_folder/\n' % (tempoutput)
        secondshell = '%s/%s.sh' % (scriptdir_crispr + '/secondshell/', fastaname)
        f1 = open(secondshell, 'w')
        f1.write(cmds)
        f1.close()
        os.system('chmod u+x %s' % (secondshell))
        singularityshell = '%s/%s.sh' % (scriptdir_crispr + '/singularity/', i%(250))
        cmds = 'mkdir %s_temp\ncd %s_temp\n' % (tempoutput,tempoutput)
        cmds += 'singularity exec -B /scratch/users/anniz44/ -H /home/anniz44/ /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CrisprCasFinder.simg %s\n' % (
            secondshell)
        cmds += 'rm -r %s_temp/\n' % (tempoutput)
        f1 = open(singularityshell, 'a')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\nmodule add c3ddb/singularity/3.5.2\n')
        f1.write(cmds)
        f1.close()
        i += 1
    return i

# sum all sigularity codes in c7
fscript = open(os.path.join(scriptdir, 'allruncrispr.sh'), 'w')
fscript.write('#!/bin/bash\nsource ~/.bashrc\n')
# generate code for each fasta
singularityshellall = set()
i = 0
for fasta in allfasta:
    i = runCRISPRCasFinder(fasta, i)

allsingularityshell = glob.glob('%s/*.sh'%(scriptdir_crispr + '/singularity/'))
for singularityshell in allsingularityshell:
    fscript.write('jobmit %s %s c7\n' % (singularityshell, os.path.split(singularityshell)[-1]))
fscript.close()
print('please run: sh allruncrispr.sh')

