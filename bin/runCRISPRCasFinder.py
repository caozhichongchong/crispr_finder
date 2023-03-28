import glob,os

#allfasta = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/allfasta/*/fasta/*.fasta')
allfasta = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/allfasta/*/*.fasta')
outputdir = '/scratch/users/anniz44/genomes/crispr_MG/CRISPRCasFinder/'
scriptdir = '/scratch/users/anniz44/scripts/1MG/crispr/'
scriptdir_crispr = '%s/CRISPRCasFinderscript/'%(scriptdir)

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


def runCRISPRCasFinder(fasta, i):
    fastaname = os.path.split(fasta)[-1]
    tempoutput = '%s/%s' % (outputdir, fastaname)
    try:
        f1 = open('%s.result.json'%(tempoutput),'r')
    except IOError:
        cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\nexport PATH=/scratch/users/anniz44/bin/pro/vmatch-2.3.1-Linux_x86_64-64bit:$PATH'
        cmds += '\n/usr/bin/perl /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CRISPRCasFinder.pl -so /scratch/users/anniz44/bin/pro/CRISPRCasFinder/sel392v2.so -cf /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CasFinder-2.0.3 -drpt /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G -out %s_folder -in %s\n' % (
            tempoutput, fasta
        )
        cmds += 'mv %s_folder/result.json %s.result.json\n' % (tempoutput, tempoutput)
        cmds += 'rm -r %s_folder/\n' % (tempoutput)
        secondshell = '%s/%s.sh' % (scriptdir_crispr + '/secondshell/', fastaname)
        f1 = open(secondshell, 'w')
        f1.write(cmds)
        f1.close()
        os.system('chmod u+x %s' % (secondshell))
        singularityshell = '%s/%s.sh' % (scriptdir_crispr + '/singularity/', i%25)
        if singularityshell not in singularityshellall:
            f1 = open(singularityshell, 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\nmodule add c3ddb/singularity/3.5.2\n')
            f1.close()
            fscript.write('jobmit %s %s c7\n' % (singularityshell, os.path.split(singularityshell)[-1]))
            singularityshellall.add(singularityshell)
        cmds = 'mkdir %s_temp\ncd %s_temp\n' % (tempoutput,tempoutput)
        cmds += 'singularity exec -B /scratch/users/anniz44/ -H /home/anniz44/ /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CrisprCasFinder.simg %s\n' % (
            secondshell)
        cmds += 'rm -r %s_temp/\n' % (tempoutput)
        f1 = open(singularityshell, 'a')
        f1.write(cmds)
        f1.close()
        i += 1
    return i


# sum all codes in c7
fscript = open(os.path.join(scriptdir, 'allruncrispr.sh'), 'w')
fscript.write('#!/bin/bash\nsource ~/.bashrc\n')
# generate code for each fasta
singularityshellall = set()
i = 0
for fasta in allfasta:
    i = runCRISPRCasFinder(fasta, i)

fscript.close()
print('please run: sh allruncrispr.sh')
