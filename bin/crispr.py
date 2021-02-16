################################################### SET PATH ########################################################
# collect all WGS
import os,glob
# move formatted fastq
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species'

for folders in glob.glob('%s/round2/*'%(WGS_dir)):
    folders_short = folders.split('_cluster')[0].replace('round2','round1')
    os.system('mv %s/fastq/* %s/fastq/'%(folders,folders_short))

for folders in glob.glob('%s/round3/*'%(WGS_dir)):
    folders_short = folders.split('_cluster')[0].replace('round3','round1')
    os.system('mv %s/fastq/* %s/fastq/'%(folders,folders_short))

# move unformatted fastq
WGS_dir2 = '/scratch/users/anniz44/genomes/GHM/newgenomes/BN10_WG'
for files in glob.glob('%s/*/*.concat.trim.fastq'%(WGS_dir2)):
    folders_short,fastq_name = os.path.split(files)
    donor = fastq_name.split('_')[0]
    folders_short = os.path.split(folders_short)[-1]
    genome_folder_name_set = folders_short.split('_')
    newgenome_folder_name = '%s_%s%s%s' % (donor, genome_folder_name_set[0][0:2],
                                        genome_folder_name_set[1][0:1].upper(), genome_folder_name_set[1][1:2])
    newgenome_folder_name = newgenome_folder_name.replace('BiAd', 'BA').replace('BiLo', 'BL')
    newgenome_folder = '%s/round1/%s'%(WGS_dir,newgenome_folder_name)
    try:
        os.mkdir(newgenome_folder)
    except IOError:
        pass
    os.system('mv %s %s/'%(files,newgenome_folder))

# changename
changename_file = '%s/file.change.name.txt'%(WGS_dir)
changename = dict()
for lines in open(changename_file,'r'):
    lines_set = lines.split('\t')
    genome_folder_name_set = lines_set[1].split('_')
    changename.setdefault(lines_set[0],'%s_%s%s%s_g%.4d' % (genome_folder_name_set[0], genome_folder_name_set[1][0:2],
                                        genome_folder_name_set[2][0:1].upper(), genome_folder_name_set[2][1:2],
                          int(genome_folder_name_set[-1].split('g')[1])))

for folders in glob.glob('%s/round1/*'%(WGS_dir)):
    allfastq = glob.glob('%s/*R1.concat.trim.fastq'%(folders))
    if allfastq!=[]:
        allin = glob.glob('%s/fastq/*_1.fastq'%(folders))
        i = 0
        if allin != []:
            for allinfastq in allin:
                i = max(i,int(os.path.split(allinfastq)[-1].split('_g')[-1].split('_1.fastq')[0]))
        i += 1
        folders_short = os.path.split(folders)[-1]
        try:
            os.mkdir(folders + '/fastq/')
        except IOError:
            pass
        for fastq in allfastq:
            oldname = os.path.split(fastq)[-1].split('R1.concat.trim.fastq')[0]
            newname = '%s_g%.4d'%(folders_short,i)
            changename.setdefault(oldname,newname)
            os.system('mv %s %s/fastq/%s_1.fastq'%(fastq,folders,newname))
            os.system('mv %s %s/fastq/%s_2.fastq' % (fastq.replace('R1.concat.trim.fastq',
                                                                   'R2.concat.trim.fastq'), folders, newname))
            i += 1


# output changename
changename_output = set()
for oldname in changename:
    newname = changename[oldname].replace('BiAd', 'BA').replace('BiLo', 'BL')
    changename[oldname] = newname
    changename_output.add('%s\t%s\t\n'%(oldname,newname))

f1 = open('%s/file.change.name.new.txt'%(WGS_dir),'w')
f1.write(''.join(list(changename_output)))
f1.close()

# changename fasta
genome_dir1 = ('/scratch/users/anniz44/genomes/GHM/newgenomes/BN10_Genome_Library_07052018/')
allfasta = glob.glob('%s/*_final.scaffolds.fasta'%(genome_dir1))
allshortname = ['_prokka.faa','_prokka.faa.add','_prokka.gff']
for fasta in allfasta:
    folder, oldname = os.path.split(fasta)
    oldname = oldname.split('_final.scaffolds.fasta')[0]
    newname = changename.get(oldname,'None')
    if newname != 'None':
        newfolder = '_'.join(newname.split('_')[0:2])
        newfolder = '%s/round1/%s/fasta/'%(WGS_dir,newfolder)
        try:
            os.mkdir(newfolder)
        except IOError:
            pass
        os.system('mv %s %s/%s_final.scaffolds.fasta'%(fasta,newfolder,
                                                   newname))
        for shortname in allshortname:
            os.system('mv %s %s/%s%s' % (fasta.replace('_final.scaffolds.fasta',
                                                                       shortname), newfolder,
                                                         newname,shortname))

# check fastq
output_script = '/scratch/users/anniz44/scripts/1MG/virus/runspades/'
try:
    os.mkdir(output_script)
except IOError:
    pass

def runspades(file1,file2,temp_output,output_name):
    cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
    cmds += 'spades.py --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            (file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -r %s\n' % (temp_output)
    cmds += 'prokka --kingdom Bacteria --outdir %s --locustag Bacter %s\n' % \
                 (temp_output, output_name)
    cmds += 'mv %s/*.gff %s\n'%(temp_output,output_name.replace('_final.scaffolds.fasta','_prokka.gff'))
    cmds += 'mv %s/*.faa %s\n' % (temp_output, output_name.replace('_final.scaffolds.fasta', '_prokka.faa'))
    cmds += 'rm -r %s\n' % (temp_output)
    return cmds

for folders in glob.glob('%s/round1/*'%(WGS_dir)):
    allfastq = glob.glob('%s/fastq/*_1.fastq'%(folders))
    for fastq in allfastq:
        fastqname = os.path.split(fastq)[-1].split('_1.fastq')[0]
        try:
            os.mkdir('%s/fasta/'%(folders))
        except IOError:
            pass
        fasta = glob.glob('%s/fasta/%s_final.scaffolds.fasta'%(folders,fastqname))
        fasta2 = glob.glob('%s/fasta/%s_prokka.gff' % (folders, fastqname))
        if fasta == [] or fasta2 == []:
            print('missing fasta %s in %s/fasta/%s*'%(fastqname,folders,fastqname))
            cmds = runspades(fastq, fastq.replace('_1.fastq','_2.fastq'), '%s/fasta/prokka_%s'%(folders,fastqname),
                             '%s/fasta/%s_final.scaffolds.fasta'%(folders,fastqname))
            f1 = open('%s/%s.sh'%(output_script,fastqname),'w')
            f1.write(cmds)
            f1.close()

f1 = open(os.path.join(output_script, '../allspades.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(output_script, '*.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (output_script, '../allspades.sh'))

################################################### END ########################################################
################################################### SET PATH ########################################################
# run crass on WGS
# rerun crass
import os,glob
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details/'
output_script = '/scratch/users/anniz44/scripts/1MG/virus/crass_MG/'
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'

try:
    os.mkdir(output_script)
except IOError:
    pass

try:
    os.mkdir(output_folder + '/../')
except IOError:
    pass

try:
    os.mkdir(output_folder)
except IOError:
    pass

allfastq = glob.glob('%s/am_*/fastq/*_1.fastq'%(WGS_dir))
for fastq in allfastq:
    filedir, filename = os.path.split(fastq)
    fastq2 = fastq.replace('_1.fastq','_2.fastq')
    output_file = os.path.join(output_folder,filename)
    f1 = open('%s/%s' % (output_script, filename) + '.sh', 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    f1.write('/scratch/users/anniz44/bin/pro/crass-0.3.12/src/crass/crass -o ' + output_file +
             ' -l 4 ' + fastq + ' ' + fastq2 + '\n')
    f1.write(
        '/scratch/users/anniz44/bin/pro/crisprtools/src/crisprtools extract -s  %s/crass.crispr  > %s/spacers.fasta\n' % (
            output_file,
            output_file))
    f1.write(
        '/scratch/users/anniz44/bin/pro/crisprtools/src/crisprtools extract -d  %s/crass.crispr  > %s/DR.fasta\n' % (
            output_file,
            output_file))
    f1.close()

f1 = open(os.path.join(output_script, '../allcrass.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(output_script, '*.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (output_script, '../allcrass.sh'))


################################################### END ########################################################
################################################### SET PATH ########################################################
# run crisprfinder on assembly
import os,glob
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details/'
output_script = '/scratch/users/anniz44/scripts/1MG/virus/crass_assembly/'
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
merge_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/fastawithcrispr/'

try:
    os.mkdir(merge_dir)
except IOError:
    pass

allfastqwithcrispr = glob.glob('%s/*_1.fastq'%(output_folder))
for fastqfile in allfastqwithcrispr:
    fasta = glob.glob('%s/am_*/fasta/%s' % (WGS_dir,
                                            os.path.split(fastqfile)[-1].replace('_1.fastq','_final.scaffolds.fasta')))
    os.system('cp %s %s/'%(fasta[0],merge_dir))

try:
    os.mkdir(output_script)
except IOError:
    pass

try:
    os.mkdir(output_folder + '/../')
except IOError:
    pass

try:
    os.mkdir(output_folder)
except IOError:
    pass


allfastqwithcrispr = glob.glob('%s/*_1.fastq'%(output_folder))
for fastqfile in allfastqwithcrispr:
    fasta = glob.glob('%s/am_*/fasta/%s' % (WGS_dir,
                                            os.path.split(fastqfile)[-1].replace('_1.fastq','_final.scaffolds.fasta')))[0]
    filename = os.path.split(fasta)[-1]
    output_file = os.path.join(output_folder,filename)
    f1 = open('%s/%s' % (output_script, filename) + '.sh', 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nmodule add c3ddb/singularity/3.5.2\n')
    f1.write('singularity exec -B /scratch/users/anniz44/bin/pro/CRISPRCasFinder/ /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CrisprCasFinder.simg perl /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CRISPRCasFinder.pl -so /scratch/users/anniz44/bin/pro/CRISPRCasFinder/sel392v2.so '+
             '-cf /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CasFinder-2.0.3 -drpt /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G '+
             '-out %s -in %s'%(output_file,fasta))
    f1.close()

f1 = open(os.path.join(output_script, '../allcrisprfinder.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(output_script, '*.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (output_script, '../allcrisprfinder.sh'))

################################################### END ########################################################
################################################### SET PATH ########################################################
# remove empty folder
os.system('#sh DR.sh')
os.system('#rsync -a *.fasta details/')

import os,glob
input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details'
try:
    os.mkdir('%s/nocrispr2'%(input_folder))
except IOError:
    pass

for files in glob.glob('%s/*'%(input_folder)):
    if 'nocrispr2' not in files:
        try:
            f1 = open('%s/crass.crispr'%(files),'r')
            if int(os.path.getsize('%s/spacers.fasta'%(files)))==0:
                os.system('rm -r %s/nocrispr2/%s/*' % (input_folder, os.path.split(files)[-1]))
                os.system('mv %s %s/nocrispr2/'%(files,input_folder))
                print(files)
        except IOError:
            os.system('rm -r %s/nocrispr2/%s/*' % (input_folder,os.path.split(files)[-1]))
            os.system('mv %s %s/nocrispr2/'%(files,input_folder))
            print(files)


################################################### END ########################################################
################################################### SET PATH ########################################################
# merge crispr spacers
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details/'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/'
filename_format = '_1.fastq'
seq_limit1 = 27
seq_limit2 = 43

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

def add_name(fastafile,outputlist,Spacers = False):
    filefolder, filename = os.path.split(fastafile)
    filefolder, filename = os.path.split(filefolder)
    filename = filename.replace(filename_format,'')
    oldgroup, oldID = [0,0]
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        if Spacers:
            record_id, oldgroup, oldID = changespacername(record_id,oldgroup, oldID)
        record_seq = str(record.seq)
        if not Spacers or qualify_len(record_seq):
            outputlist.append('>%s__%s\n%s\n'%(filename,record_id,record_seq))
    return outputlist

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

allspacerfiles = glob.glob('%s/*/spacers.fasta'%(input_folder))
allDRfiles = glob.glob('%s/*/DR.fasta'%(input_folder))
# spacer merge
alloutput = []
outputfilename = '%s/all.spacers.fasta'%(output_folder)
for files in allspacerfiles:
    alloutput = add_name(files, alloutput, True)

output_fasta(alloutput,outputfilename)
# cluster spacer seq
for cutoff in [0.95,0.9,0.85,0.8]:
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.%s.fasta -uc %s.cluster.%s.uc -threads %s\n'
                       % ('usearch', outputfilename, cutoff, outputfilename,cutoff,
                          outputfilename,cutoff, 40))
    os.system(cmd_cluster)

# DR merge
alloutput = []
outputfilename = '%s/all.DR.fasta'%(output_folder)
for files in allDRfiles:
    alloutput = add_name(files, alloutput, False)

output_fasta(alloutput,outputfilename)

################################################### END ########################################################
################################################### SET PATH ########################################################
# map spacers to phage genomes
import os,glob

input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/'
database = '/scratch/users/anniz44/scripts/database/virus/phage.ENA.fasta'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/anno_virusMG/'

try:
    os.mkdir(input_script)
except IOError:
    pass

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = ''
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds = args.bw + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, files, files2, args.sam, min(40, args.t),
            tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam, min(40, args.t),
            tempbamoutput)
        cmds += 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
            tempbamoutput, tempbamoutput)
        cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            args.bcf, min(40, args.t), database,
            tempbamoutput, args.bcf, min(40, args.t), vcfoutput)
        cmds += '#%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
            args.bcf, vcfoutput, vcfoutput)
    return cmds

def usearch(database,fasta):
    # phage db
    cmds = ('blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 100 -num_threads 40 -qcov_hsp_perc 90\n'%(database,fasta,fasta+ '.phage.txt'))
    # phage db
    cmds += ('blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 100 -num_threads 40 -perc_identity 90\n' % (
        database, fasta, fasta + '.phage.all.txt'))
    # nt db
    cmds += ('#blastn -db nt -remote -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 5 -qcov_hsp_perc 90 -perc_identity 90 -word_size 20\n' % (
    fasta, fasta + '.phage.nt.txt'))
    return cmds

def usearch_self(fasta):
    try:
        f1 = open('%s.nhr'%(fasta),'r')
    except IOError:
        os.system('makeblastdb -in %s -dbtype nucl'%(fasta))
    cmds = ('blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 1000 -num_threads 40 -qcov_hsp_perc 90 -perc_identity 90 \n'%(fasta,fasta,fasta+ '.selfblast.txt'))
    return cmds

allspacerfiles = glob.glob('%s/*.fasta.cluster.*.fasta'%(input_folder))
for spacerfile in allspacerfiles:
    cmds = usearch(database, spacerfile)
    f1 = open(os.path.join(input_script, '%s.sh'%(os.path.split(spacerfile)[-1])), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()

allspacerfiles = glob.glob('%s/all.spacers.fasta'%(input_folder))
for spacerfile in allspacerfiles:
    cmds = usearch_self(spacerfile)
    f1 = open(os.path.join(input_script, '%s.sh'%(os.path.split(spacerfile)[-1])), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# map spacers to bacterial genomes -> after runspades
import os,glob

input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details/'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/anno_virusam/'
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'

try:
    os.mkdir(input_script)
except IOError:
    pass

def usearch(database,fasta):
    # bacterial db
    try:
        f1 = open('%s.nhr'%(database),'r')
    except IOError:
        os.system('makeblastdb -in %s -dbtype nucl'%(database))
    try:
        f1 = open(fasta+ '.bac.txt','r')
        cmds = ''
    except IOError:
        cmds = ('blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 10 -num_threads 40 -qcov_hsp_perc 90 -perc_identity 90\n'%(database,fasta,fasta+ '.bac.txt'))
    return cmds

allspacerfiles = glob.glob('%s/*/spacers.fasta'%(input_folder))
for spacerfile in allspacerfiles:
    filefolder, filename = os.path.split(spacerfile)
    filefolder, filename = os.path.split(filefolder)
    genomefile = glob.glob('%s/%s_%s/fasta/%s'%(WGS_dir,filename.split('_')[0],filename.split('_')[1],
                                                filename.replace('_1.fastq',
                                                                 '_final.scaffolds.fasta')))
    if genomefile!=[]:
        cmds = usearch(genomefile[0], spacerfile)
        f1 = open(os.path.join(input_script, '%s.sh'%(os.path.split(filename)[-1])), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
        f1.close()
    else:
        print('missing file ','%s/%s_%s/fasta/%s'%(WGS_dir,filename.split('_')[0],filename.split('_')[1],
                                                filename.replace('_1.fastq',
                                                                 '_final.scaffolds.fasta')))

# all scripts
f1 = open(os.path.join(input_script, '../allanno_virusam.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script, '*_1.fastq.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (input_script, '../allanno_virusam.sh'))

################################################### END ########################################################
################################################### SET PATH ########################################################
# summarize map spacers to bacterial genomes
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details//'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
filename_format = '_1.fastq'
distance_spacer = 5000
seq_limit1 = 27
seq_limit2 = 43

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

def add_name(inputfile,outputlist):
    fastafile = inputfile.split('.bac.txt')[0]
    filefolder, filename = os.path.split(inputfile)
    filefolder, filename = os.path.split(filefolder)
    filename = filename.replace(filename_format,'')
    oldgroup, oldID = [0,0]
    newname_map = dict()
    crispr_contig = dict()
    # new name
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        newrecord_id, oldgroup, oldID = changespacername(record_id,oldgroup, oldID)
        record_seq = str(record.seq)
        if qualify_len(record_seq):
            newname_map.setdefault(record_id,'%s__%s'%(filename,newrecord_id))
    # remove crispr parts on contig
    for lines in open(inputfile, 'r'):
        lines_set = lines.split('\t')
        contig = lines_set[1]
        locus = int(lines_set[8])
        crispr_contig.setdefault(contig,[])
        crispr_contig[contig].append(locus)
    for contig in crispr_contig:
        alllocus = crispr_contig[contig]
        alllocus.sort()
        alllocus_crispr = set()
        for i in range(1,len(alllocus)):
            if abs(alllocus[i]-alllocus[i-1])<=distance_spacer:
                alllocus_crispr.add(alllocus[i])
                alllocus_crispr.add(alllocus[i-1])
        if len(alllocus_crispr)>1:
            crispr_contig[contig] = [min(alllocus_crispr),max(alllocus_crispr)]
        else:
            crispr_contig[contig] = [0,0]
    # filter true self-targeting
    for lines in open(inputfile,'r'):
        record_id = lines.split('\t')[0]
        newrecord_id = newname_map.get(record_id)
        if newrecord_id != None:
            contig = lines_set[1]
            locus = int(lines_set[8])
            if contig not in crispr_contig or locus > crispr_contig[contig][1] or locus < crispr_contig[contig][0]:
                outputlist.append('%s\t%s'%(newrecord_id,lines))
    return outputlist

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

allspacerfiles = glob.glob('%s/*/spacers.fasta.bac.txt'%(input_folder))
alloutput = []
outputfilename = '%s/all.spacers.fasta.bac.txt'%(output_folder)
for files in allspacerfiles:
    alloutput = add_name(files, alloutput)

output_fasta(alloutput,outputfilename)

################################################### END ########################################################
################################################### SET PATH ########################################################
# annotate map bac seq
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details//'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'

filename_format = '_1.fastq'
outputfilename = '%s/all.spacers.fasta.bac.txt'%(output_folder)
allbac = dict()
def mapping_process(allmapping):
    Contig_all = dict()
    for mapping in allmapping:
        newname, oldname, contig, identity, hitlength = mapping[0:5]
        contig = contig.replace('|','_')
        hit_left, hit_right = mapping[9:11]
        mid_hit = (float(hit_left) + float(hit_right))/2
        e_value = mapping[11]
        bit_score = mapping[12].split('\n')[0]
        Contig_all.setdefault(contig,[[],[]])
        Contig_all[contig][0].append(mid_hit)
        newline = '%s\t%s\t%s\t%s\t%s\t%s\t%s'%(newname,contig,mid_hit,bit_score,identity, hitlength,e_value)
        Contig_all[contig][1].append(newline)
    return Contig_all

# read all bac mapping
for lines in open(outputfilename,'r'):
    lines_set = lines.split('\t')
    genomes = lines_set[0].split('__')[0]
    allbac.setdefault(genomes,[])
    allbac[genomes].append(lines_set)

# annotate bac mapping
Output = []
for genomes in allbac:
    allmapping = allbac[genomes]
    gff_file = glob.glob('%s/%s_%s/fasta/%s_prokka.gff' % (WGS_dir, genomes.split('_')[0], genomes.split('_')[1],
                                                genomes))
    gff_file = gff_file[0]
    Contig_all = mapping_process(allmapping)
    for lines in open(gff_file,'r'):
        if not lines.startswith('#') and not lines.startswith('>'):
            try:
                lines_set = lines.split('\t')
                contig = lines_set[0]
                hit_left, hit_right = lines_set[3:5]
                hit_left = float(hit_left)
                hit_right = float(hit_right)
                if contig in Contig_all:
                    allmid = Contig_all[contig][0]
                    for i in range(0,len(allmid)):
                        if (hit_left < allmid[i] and  allmid[i] < hit_right) or (hit_left > allmid[i] and  allmid[i] > hit_right):
                            print(contig, hit_left, hit_right,allmid[i])
                            annotation = lines_set[8]
                            try:
                                name_gene = annotation.split('Name=')[1].split(';')[0]
                            except IndexError:
                                name_gene = ''
                            try:
                                gene_name = annotation.split('gene=')[1].split(';')[0]
                            except IndexError:
                                gene_name = ''
                            try:
                                product = annotation.split('product=')[1].split('\n')[0]
                            except IndexError:
                                product = ''
                            try:
                                note = annotation.split('note=')[1].split(';')[0]
                            except IndexError:
                                note = ''
                            Contig_all[contig][1][i] += '\t%s\t%s\t%s\t%s'%(name_gene,gene_name,product,note)
                            print(Contig_all[contig][1][i] )
            except ValueError:
                pass
    for contig in Contig_all:
        for anno in Contig_all[contig][1]:
            Output.append(anno)

f1 = open(outputfilename + '.anno.txt','w')
f1.write('newname\tcontig\tmid_hit\tbit_score\tidentity\thitlength\te_value\tname_gene\tgene_name\tproduct\tnote\n')
f1.write('\n'.join(Output))
f1.close()
# excluding annotation as CRISPR
os.system('grep -v \"CRISPR\" %s.anno.txt > %s.anno.nocrispr.txt'%(outputfilename,outputfilename))
################################################### END ########################################################
################################################### SET PATH ########################################################
# phage seq from MG
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

input_fastq = '/scratch/users/anniz44/Metagenomes/BN10_MG/fastq/'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/MG_phage/'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/MG_phage/'

fq_format = '_1.fq'
fq_format2 = '_2.fq'

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
    filename = os.path.split(fastq_1)[-1]
    output_file = '%s/%s/'%(output_folder,filename)
    cmd = 'export HDF5_USE_FILE_LOCKING=\'FALSE\'\n'
    # megahit
    cmd += 'megahit -1 %s -2 %s -o %s --min-contig-len 1000 -t 40 -m 0.99\n'%(fastq_1,fastq_2,output_file)
    # phage seq from MG
    cmd += 'predict-metagenome %s/final.contigs.fa\n'%(output_file)
    f1 = open(os.path.join(input_script, '%s.sh'%(filename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s'%(cmd))
    f1.close()

# all scripts
f1 = open(os.path.join(input_script, '../allMG_phage.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script, '*.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (input_script, '../allMG_phage.sh'))
################################################### END ########################################################
################################################### SET PATH ########################################################
# phage map spacers to phage in MG
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

MG_phage = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/MG_phage/*/final.contigs.fa'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
spacersfile = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/all.spacers.fasta.cluster.0.95.fasta'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/map_MGphage/'
phage_score = 0.5
try:
    os.mkdir(input_script)
except IOError:
    pass

def extract_phage(phage_file):
    qualify_phage = dict()
    for lines in open(phage_file+'.phage.txt','r'):
        if not lines.startswith('name'):
            lines_set = lines.split('\t')
            score = float(lines_set[-1].split('\n')[0])
            if score > phage_score:
                qualify_phage.setdefault(lines_set[0],score)
    qualify_phage_fa = set()
    for record in SeqIO.parse(phage_file, 'fasta'):
        record_id = str(record.id)
        if record_id in qualify_phage:
            qualify_phage_fa.add('>%s__%s\n%s\n'%(record_id,
                                                 qualify_phage[record_id], str(record.seq)))
    f1 = open(phage_file+'.phage.fa','w')
    f1.write(''.join(list(qualify_phage_fa)))
    f1.close()

def usearch(database,fasta):
    database_name = os.path.split(os.path.split(database)[0])[-1]
    try:
        f1 = open('%s.nhr'%(database),'r')
    except IOError:
        os.system('makeblastdb -in %s -dbtype nucl'%(database))
    cmds = 'blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 100 -num_threads 40 -qcov_hsp_perc 90\n'%(
        database,fasta,fasta+ '.MG.%s.phage.txt'%(database_name))
    cmds += 'blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 100 -num_threads 40 -perc_identity 90\n' % (
        database, fasta, fasta + '.all.MG.%s.phage.txt' %(database_name))
    return cmds

allMG_phage = glob.glob(MG_phage)
for phage_file in allMG_phage:
    extract_phage(phage_file)
    cmds = usearch(phage_file + '.phage.fa', spacersfile)
    f1 = open(os.path.join(input_script, '%s.sh'%(os.path.split(os.path.split(phage_file)[0])[-1])), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# phage map spacers to plasmid
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

plasmid_db = '/scratch/users/anniz44/genomes/NCBI_Genome/newplasmid/*.fna'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/details/'
spacersfile = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/anno_details/all.spacers.fasta.cluster.0.95.fasta'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/map_plasmid/'

try:
    os.mkdir(input_script)
except IOError:
    pass

def usearch(database,fasta):
    database_name = os.path.split(database)[-1]
    try:
        f1 = open('%s.nhr'%(database),'r')
    except IOError:
        os.system('makeblastdb -in %s -dbtype nucl'%(database))
    cmds = 'blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 100 -num_threads 40 -qcov_hsp_perc 90\n'%(
        database,fasta,fasta+ '.plasmid.%s.txt'%(database_name))
    cmds += 'blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 100 -num_threads 40 -perc_identity 90\n' % (
        database, fasta, fasta + '.all.plasmid.%s.txt' %(database_name))
    return cmds

allplasmid = glob.glob(plasmid_db)
for plasmid_file in allplasmid:
    cmds = usearch(plasmid_file, spacersfile)
    f1 = open(os.path.join(input_script, '%s.sh'%(os.path.split(plasmid_file)[-1])), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# prophage
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/prophage/'
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
filename_format = '_final.scaffolds.fasta'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/prophage/'

os.system('rm -r %s'%(input_script))
os.mkdir(input_script)
try:
    os.mkdir(output_folder)
except IOError:
    pass

def find_prophage(fasta,outputdir):
    cmds = 'phigaro -f %s -o %s -t 40 --save-fasta -d --not-open -c /scratch/users/anniz44/bin/pro/polly/.phigaro/config.yml -e tsv gff\n'%(fasta,outputdir)
    return cmds

allgenome = glob.glob('%s/*/fasta/*%s'%(WGS_dir,filename_format))
i = 0
for fasta in allgenome:
    filename = os.path.split(fasta)[-1]
    output_file = '%s/%s/'%(output_folder,filename)
    # megahit
    cmd = find_prophage(fasta,output_file)
    f1 = open(os.path.join(input_script, '%s.sh'%(int(i / 10))), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s'%(cmd))
    f1.close()
    i += 1

# all scripts
f1 = open(os.path.join(input_script, '../allprophage.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script, '*.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (input_script, '../allprophage.sh'))
################################################### END ########################################################
################################################### SET PATH ########################################################
# check prophage results and re-run spades + prophage
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/prophage/'
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
filename_format = '_final.scaffolds.fasta'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/prophage_rerun/'

try:
    os.mkdir(input_script)
except IOError:
    pass

def find_prophage(fasta,outputdir):
    cmds = 'phigaro -f %s -o %s -t 40 --save-fasta -d --not-open -c /scratch/users/anniz44/bin/pro/polly/.phigaro/config.yml -e tsv gff\n'%(fasta,outputdir)
    return cmds

def runspades(file1,file2,temp_output,output_name):
    cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
    cmds += 'spades.py --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            (file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -r %s\n' % (temp_output)
    cmds += 'prokka --kingdom Bacteria --outdir %s --locustag Bacter %s\n' % \
                 (temp_output, output_name)
    cmds += 'mv %s/*.gff %s\n'%(temp_output,output_name.replace('_final.scaffolds.fasta','_prokka.gff'))
    cmds += 'mv %s/*.faa %s\n' % (temp_output, output_name.replace('_final.scaffolds.fasta', '_prokka.faa'))
    cmds += 'rm -r %s\n' % (temp_output)
    return cmds

# check failed jobs
need_rerun = set()
for folders in glob.glob('%s/*%s'%(output_folder,filename_format)):
    output_file = glob.glob('%s/*.gff3'%(folders))
    if output_file == []:
        # failed job
        need_rerun.add(os.path.split(folders)[-1].split(filename_format)[0])

print(need_rerun)

for genome in need_rerun:
    fastq = glob.glob('%s/*/fastq/%s_1.fastq'%(WGS_dir,genome))[0]
    fastqname = os.path.split(fastq)[-1].split('_1.fastq')[0]
    folders = os.path.split(fastq)[0]
    folders = os.path.split(folders)[0]
    fasta = '%s/fasta/%s%s' % (folders, fastqname,filename_format)
    cmds = runspades(fastq, fastq.replace('_1.fastq', '_2.fastq'), '%s/fasta/prokka_%s' % (folders, fastqname),
                     fasta)
    output_file = '%s/%s/' % (output_folder, fastqname + filename_format)
    cmds += find_prophage(fasta, output_file)
    f1 = open('%s/%s.sh' % (input_script, fastqname), 'w')
    f1.write(cmds)
    f1.close()

# all scripts
f1 = open(os.path.join(input_script, '../allprophage_rerun.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script, '*.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (input_script, '../allprophage_rerun.sh'))

################################################### END ########################################################
################################################### SET PATH ########################################################
# map crispr to prophage
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/prophage/'
WGS_dir ='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
filename_format = '_final.scaffolds.fasta'
input_script = '/scratch/users/anniz44/scripts/1MG/virus/map_prophage/'
input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details//'

try:
    os.mkdir(input_script)
except IOError:
    pass

def usearch(database,fasta):
    # bacterial db
    try:
        f1 = open(fasta+ '.prophage.all.txt','r')
        cmds = ''
    except IOError:
        try:
            f1 = open('%s.nhr' % (database), 'r')
            cmds = ''
        except IOError:
            cmds = ('makeblastdb -in %s -dbtype nucl\n' % (database))
        cmds += ('#blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 10 -num_threads 40 -qcov_hsp_perc 90 -perc_identity 90\n'%(database,fasta,fasta+ '.prophage.txt'))
        cmds += (
                    'blastn -db %s -task blastn-short -query %s -out %s -outfmt 6 -max_target_seqs 10 -num_threads 40 -qcov_hsp_perc 90\n' % (
            database, fasta, fasta + '.prophage.all.txt'))
    return cmds

allspacerfiles = glob.glob('%s/*/spacers.fasta'%(input_folder))
for spacerfile in allspacerfiles:
    filefolder, filename = os.path.split(spacerfile)
    filefolder, filename = os.path.split(filefolder)
    filename = filename.replace('_1.fastq','_final.scaffolds.fasta')
    prophage_file = glob.glob('%s/%s/*.fasta'%(output_folder,filename))
    if prophage_file!=[]:
        cmds = usearch(prophage_file[0], spacerfile)
        f1 = open(os.path.join(input_script, '%s.sh'%(filename)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
        f1.close()
    else:
        print('missing file ','%s/%s/*.fasta'%(output_folder,filename))

# all scripts
f1 = open(os.path.join(input_script, '../allmapprophage.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script, '*.sh')):
    try:
        f2 = open('%s.err'%(sub_scripts),'r')
    except IOError:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s' % (input_script, '../allmapprophage.sh'))

################################################### END ########################################################
################################################### SET PATH ########################################################
# summarize map spacers to prophages
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details//'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
filename_format = '_1.fastq'
distance_spacer = 5000
seq_limit1 = 27
seq_limit2 = 43
kmer_cutoff = 16

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

def add_name(inputfile,outputlist,kmer = False):
    fastafile = inputfile.split('.prophage.')[0]
    filefolder, filename = os.path.split(inputfile)
    filefolder, filename = os.path.split(filefolder)
    filename = filename.replace(filename_format,'')
    oldgroup, oldID = [0,0]
    newname_map = dict()
    # new name
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        newrecord_id, oldgroup, oldID = changespacername(record_id,oldgroup, oldID)
        record_seq = str(record.seq)
        if qualify_len(record_seq):
            newname_map.setdefault(record_id,'%s__%s'%(filename,newrecord_id))
    # filter true self-targeting
    for lines in open(inputfile,'r'):
        record_id = lines.split('\t')[0]
        newrecord_id = newname_map.get(record_id)
        if newrecord_id != None and (not kmer or int(lines.split('\t')[3]) >= kmer_cutoff):
            outputlist.append('%s\t%s'%(newrecord_id,lines))
    return outputlist

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

allspacerfiles = glob.glob('%s/*/spacers.fasta.prophage.txt'%(input_folder))
alloutput = []
outputfilename = '%s/all.spacers.fasta.prophage.txt'%(output_folder)
for files in allspacerfiles:
    alloutput = add_name(files, alloutput)

output_fasta(alloutput,outputfilename)

allspacerfiles = glob.glob('%s/*/spacers.fasta.prophage.all.txt'%(input_folder))
alloutput = []
outputfilename = '%s/all.spacers.fasta.prophage.all.txt'%(output_folder)
for files in allspacerfiles:
    alloutput = add_name(files, alloutput, True)

output_fasta(alloutput,outputfilename)

################################################### END ########################################################
################################################### SET PATH ########################################################
# annotate phage and plasmid
import os,glob
import copy

input_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/details/'
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/details/'
database_map = '/scratch/users/anniz44/scripts/database/virus/phage.ENA.details.txt'
id_cutoff = 90
kmer_cutoff = 16
max_mismatch = 1

# load database
virus_anno = dict()
for lines in open(database_map,'r'):
    lines_set = lines.split('\t')
    virus_anno.setdefault(lines_set[0].split('.')[0],lines_set[-1].split('\n')[0])

# load all mapping results
allphage_cov90output = glob.glob('%s/*.fasta.phage.txt'%(output_folder))
allphage_id90output = glob.glob('%s/*.fasta.phage.all.txt'%(output_folder))
allMGphage_cov90output = glob.glob('%s/*.fasta.MG.*.phage.txt'%(output_folder))
allMGphage_id90output = glob.glob('%s/*.fasta.all.MG.*.phage.txt'%(output_folder))
allplasmid_cov90output = glob.glob('%s/*.fasta.plasmid.*.txt'%(output_folder))
allplasmid_id90output = glob.glob('%s/*.fasta.all.plasmid.*.txt'%(output_folder))

def filter_id(input_file,anno = True):
    qualified_lines = []
    for lines in open(input_file,'r'):
        lines_set = lines.split('\t')
        if float(lines_set[2])>=id_cutoff and int(lines_set[4]) <=max_mismatch:
            if anno:
                qualified_lines.append(anno_virus(lines_set))
            else:
                qualified_lines.append(lines)
    f1 = open(input_file + '.anno.txt', 'w')
    f1.write(''.join(qualified_lines))
    f1.close()

def filter_kmer(input_file,anno = True):
    qualified_lines = []
    for lines in open(input_file, 'r'):
        lines_set = lines.split('\t')
        if int(lines_set[3]) >= kmer_cutoff and int(lines_set[4]) <=max_mismatch:
            if anno:
                qualified_lines.append(anno_virus(lines_set))
            else:
                qualified_lines.append(lines)
    f1 = open(input_file + '.anno.txt', 'w')
    f1.write(''.join(qualified_lines))
    f1.close()

def anno_virus(lines_set):
    return '%s\t%s\t%s'%(lines_set[0].split('__')[0],
                         virus_anno[lines_set[1].split('.')[0]],'\t'.join(lines_set))

# annotate virus
for files in allphage_cov90output:
    filter_id(files)

for files in allphage_id90output:
    filter_kmer(files)

for files in allMGphage_cov90output:
    filter_id(files,False)

for files in allMGphage_id90output:
    filter_kmer(files,False)

for files in allplasmid_cov90output:
    filter_id(files,False)

for files in allplasmid_id90output:
    filter_kmer(files,False)

# sort and filter output by bitscore
def fiter_score(filename, col_spacer, col_anno, col_score, col_accession, tag):
    score_line = dict()
    for lines in open(filename,'r'):
        if not lines.startswith('Query'):
            lines_set = lines.split('\t')
            score = float(lines_set[col_score].split('\n')[0])
            anno = lines_set[col_anno].split('\n')[0].replace(':',' ').replace(';',' ')
            if col_accession!=False:
                accession = lines_set[col_accession].split('\n')[0].replace(':',' ').replace(';',' ')
                anno = '%s__%s'%(accession,anno)
            spacer = lines_set[col_spacer].split('\n')[0]
            if spacer not in score_line:
                score_line.setdefault(spacer,'%s:%s'%(score,anno))
            else:
                scoresetold = score_line[spacer]
                print(scoresetold)
                scoreold, annooldset = scoresetold.split(':')
                scoreold = float(scoreold)
                if score > scoreold:
                    # better hit
                    score_line[spacer] = '%s:%s'%(score,anno)
                elif score == scoreold:
                    if all(anno.split('phage')[0] != annoold.split('phage')[0] for annoold in annooldset.split(';')):
                        # new anno
                        score_line[spacer]=('%s:%s;%s' % (score, annooldset, anno))
                else:
                    # worse hit
                    pass
    newoutput = set()
    for spacer in score_line:
        scoresetold = score_line[spacer]
        scoreold, annooldset = scoresetold.split(':')
        newoutput.add('%s\t%s\t%s\t%s\t\n'%(spacer,scoreold,annooldset,tag))
    f1 = open(filename + '.filter.txt', 'w')
    f1.write(''.join(list(newoutput)))
    f1.close()

allphage_cov90output = glob.glob('%s/*.fasta.phage.txt.anno.txt'%(output_folder))
allphage_id90output = glob.glob('%s/*.fasta.phage.all.txt.anno.txt'%(output_folder))
allMGphage_cov90output = glob.glob('%s/*.fasta.MG.*.phage.txt.anno.txt'%(output_folder))
allMGphage_id90output = glob.glob('%s/*.fasta.all.MG.*.phage.txt.anno.txt'%(output_folder))
allIMG_blast = glob.glob('%s/all.spacers.fasta.cluster.*.IMG.virusdb.blast.txt'%(output_folder))
allplasmid_cov90output = glob.glob('%s/*.fasta.plasmid.*.txt.anno.txt'%(output_folder))
allplasmid_id90output = glob.glob('%s/*.fasta.all.plasmid.*.txt.anno.txt'%(output_folder))

for files in allphage_cov90output:
    fiter_score(files,2,1,-1,3,'phagedb_id90_hit90_mis1')

for files in allphage_id90output:
    fiter_score(files,2,1,-1,3,'phagedb_id90_kmer16_mis1')

for files in allMGphage_cov90output:
    files_name = os.path.split(files)[-1].split('_1.fq')[0].split('.MG.')[1]
    fiter_score(files,0,1,-1,False,'MGphage_id90_hit90_mis1_%s'%(files_name))

for files in allMGphage_id90output:
    files_name = os.path.split(files)[-1].split('_1.fq')[0].split('.MG.')[1]
    fiter_score(files,0,1,-1,False,'MGphage_id90_kmer16_mis1_%s'%(files_name))

for files in allIMG_blast:
    fiter_score(files,0,7,14,2,'IMGphage_id90_kmer16_mis1')

for files in allplasmid_cov90output:
    fiter_score(files,0,1,-1,False,'plasmid_id90_hit90_mis1')

for files in allplasmid_id90output:
    fiter_score(files,0,1,-1,False,'plasmid_id90_kmer16_mis1')

################################################### END ########################################################
################################################### SET PATH ########################################################
# same spacers analysis (seq similarity)
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
id_cutoff = 90
hit_length = 90
kmer_cutoff = 16
cluster_cutff = 0.9

def spacer_donor(spacer):
    genome = spacer.split('__')[0]
    donor = genome.split('_')[0]
    species = genome.split('_')[1]
    return [donor,species,genome]

def filter_id_kmer(input_file):
    # same species or not
    qualified_lines_samegenome = []
    qualified_lines_samespecies = []
    qualified_lines_diffspecies = []
    allspacer_set = set()
    for lines in open(input_file,'r'):
        lines_set = lines.split('\t')
        spacer1 = lines_set[0]
        spacer2 = lines_set[1]
        if spacer1 != spacer2:
            spacer_set = [spacer1, spacer2]
            spacer_set.sort()
            spacer_set = ':'.join(spacer_set)
            if spacer_set not in allspacer_set:
                allspacer_set.add(spacer_set)
                if float(lines_set[2])>=id_cutoff and int(lines_set[3]) >= kmer_cutoff:
                    donor1, species1, genome1 = spacer_donor(spacer1)
                    donor2, species2, genome2 = spacer_donor(spacer2)
                    if genome1 == genome2:
                        qualified_lines_samegenome.append(lines)
                    elif species1 == species2:
                        qualified_lines_samespecies.append(lines)
                    else:
                        qualified_lines_diffspecies.append(lines)
    f1 = open(input_file + '.filtered.samegenome.txt', 'w')
    f1.write(''.join(qualified_lines_samegenome))
    f1.close()
    f1 = open(input_file + '.filtered.diffgenome.samespecies.txt', 'w')
    f1.write(''.join(qualified_lines_samespecies))
    f1.close()
    f1 = open(input_file + '.filtered.diffgenome.diffspecies.txt', 'w')
    f1.write(''.join(qualified_lines_diffspecies))
    f1.close()

def extract_seq(fasta,blast_result):
    select_seq = set()
    output_seq = set()
    for lines in open(blast_result,'r'):
        lines_set = lines.split('\t')
        spacer1 = lines_set[0]
        spacer2 = lines_set[1]
        select_seq.add(spacer1)
        select_seq.add(spacer2)
    for record in SeqIO.parse(fasta, 'fasta'):
        record_id = str(record.id)
        if record_id in select_seq:
            output_seq.add('>%s\n%s\n'%(record_id,str(record.seq)))
    outputfilename = blast_result + '.fasta'
    f1 = open(outputfilename, 'w')
    f1.write(''.join(list(output_seq)))
    f1.close()
    # cluster spacer seq
    cmd_cluster = (
            '%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.%s.fasta -uc %s.cluster.%s.uc -threads %s\n'
            % ('usearch', outputfilename, cluster_cutff, outputfilename, cluster_cutff,
               outputfilename, cluster_cutff, 40))
    os.system(cmd_cluster)


# filter blast result into same genome, same species, diff species
self_blast = '%s/all.spacers.fasta.selfblast.txt'%(output_folder)
filter_id_kmer(self_blast)

# extract seq, cluster
self_blast_filter = glob.glob('%s/all.spacers.fasta.selfblast.txt.*'%(output_folder))
input_fasta = '%s/all.spacers.fasta'%(output_folder)
for blast_result in self_blast_filter:
    extract_seq(input_fasta, blast_result)

################################################### END ########################################################
################################################### SET PATH ########################################################
# summarize annotation phage, bac, prophage, plasmid
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
# map annotate results
def loaduc(uc_file):
    ucseq = dict()
    ucanno = dict()
    for lines in open(uc_file,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        if lines_set[9] != '*':
            coreid = lines_set[9]
            ucseq.setdefault(lines_set[8], coreid)
        else:
            coreid = lines_set[8]
            ucseq.setdefault(coreid, coreid)
        ucanno.setdefault(coreid, [])
    return [ucseq,ucanno]

def annotate_file(anno_file,ucanno):
    for lines in open(anno_file, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        coreid = lines_set[0]
        if coreid in ucanno:
            ucanno[coreid].append('\t'.join(lines_set[1:]))
        else:
            print('missing %s'%(coreid))
    return ucanno

def annotate_fasta(fasta,ucanno,ucseq,annotated):
    ucanno_output = set()
    for record in SeqIO.parse(fasta, 'fasta'):
        record_id = str(record.id)
        if record_id in ucanno:
            coreid = record_id
        elif record_id in ucseq:
            coreid = ucseq[record_id]
        else:
            coreid = ''
        if coreid in ucanno:
            if ucanno[coreid]!= []:
                for anno in ucanno[coreid]:
                    ucanno_output.add('%s\t%s\n' % (record_id,
                                                  anno))
                annotated.add(record_id)
    f1 = open('%s.allanno.txt' % (fasta), 'a')
    f1.write(''.join(list(ucanno_output)))
    f1.close()

def non_annotate(fasta,annotated):
    ucanno_output = set()
    for record in SeqIO.parse(fasta, 'fasta'):
        record_id = str(record.id)
        if record_id not in annotated:
            ucanno_output.add('%s\tNone\n' % (record_id))
    f1 = open('%s.allanno.txt' % (fasta), 'a')
    f1.write(''.join(list(ucanno_output)))
    f1.close()

# summarize annotation results for 0.95 cluster
all_load_anno = dict()
allanno = glob.glob('%s/anno_details/all.spacers.fasta.cluster.0.95.fasta*filter.txt'%(output_folder))
uc_file = '%s/anno_details/all.spacers.fasta.cluster.0.95.uc'%(output_folder)
for anno_file in allanno:
    ucseq, ucanno = loaduc(uc_file)
    ucanno = annotate_file(anno_file, ucanno)
    all_load_anno.setdefault(anno_file,[])
    all_load_anno[anno_file] = [ucseq,ucanno]

# summarize annotation results for 0.8 cluster (IMG only)
allanno = glob.glob('%s/anno_details/all.spacers.fasta.cluster.0.8.IMG.virusdb.blast.txt.filter.txt'%(output_folder))
uc_file = '%s/anno_details/all.spacers.fasta.cluster.0.8.uc'%(output_folder)
for anno_file in allanno:
    ucseq2, ucanno2 = loaduc(uc_file)
    ucanno2 = annotate_file(anno_file, ucanno2)
    all_load_anno.setdefault(anno_file,[])
    all_load_anno[anno_file] = [ucseq2,ucanno2]

# add bac anno
bac_anno = dict()
annolist = [7,8,9,10]
anno_file = '%s/anno_details/all.spacers.fasta.bac.txt.anno.nocrispr.txt'%(output_folder)
for lines in open(anno_file,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    spacer = lines_set[0]
    anno = ''
    if len(lines_set)>7:
        for i in annolist:
            if lines_set[i]!='':
                anno += '%s;'%(lines_set[i])
    if anno!='':
        anno = [lines_set[3],  anno, 'bac_id90_hit90_mis1','']
        bac_anno.setdefault(spacer,['\t'.join(anno)+'\t'])

all_load_anno.setdefault(anno_file, [])
all_load_anno[anno_file] = [dict(),bac_anno]

# add prophage all anno
bac_anno = dict()
anno_file = '%s/anno_details/all.spacers.fasta.prophage.all.txt'%(output_folder)
for lines in open(anno_file,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    spacer = lines_set[0]
    anno = [lines_set[-1], lines_set[2], 'prophage_id90_kmer16_mis1', '']
    bac_anno.setdefault(spacer, ['\t'.join(anno) + '\t'])

all_load_anno.setdefault(anno_file, [])
all_load_anno[anno_file] = [dict(),bac_anno]

# add prophage anno
bac_anno = dict()
anno_file = '%s/anno_details/all.spacers.fasta.prophage.txt'%(output_folder)
for lines in open(anno_file,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    spacer = lines_set[0]
    anno = [lines_set[-1], lines_set[2], 'prophage_id90_hit90_mis1', '']
    bac_anno.setdefault(spacer, ['\t'.join(anno) + '\t'])

all_load_anno.setdefault(anno_file, [])
all_load_anno[anno_file] = [dict(),bac_anno]

# annotate special spacers and all spacers
allfasta = glob.glob('%s/all.spacers.fasta.selfblast.txt.*.fasta'%(output_folder)) + \
           glob.glob('%s/all.spacers.fasta'%(output_folder))
os.system('rm %s*.fasta.allanno.txt'%(output_folder))
for fasta in allfasta:
    annotated = set()
    for anno_file in all_load_anno:
        ucseq, ucanno = all_load_anno[anno_file]
        annotate_fasta(fasta, ucanno, ucseq, annotated)
        print(len(annotated),anno_file)
    non_annotate(fasta,annotated)

# sort and filter output by bitscore
def fiter_score(filename, col_spacer, col_score):
    score_line = dict()
    newoutput = set()
    emptyout = set()
    for lines in open(filename,'r'):
        if 'None' not in lines:
            lines_set = lines.split('\n')[0].split('\t')
            score = float(lines_set[col_score])
            spacer = lines_set[col_spacer]
            if spacer not in score_line:
                score_line.setdefault(spacer,'%s::%s'%(score,lines))
            else:
                scoresetold = score_line[spacer]
                print(scoresetold)
                scoreold, lines_old = scoresetold.split('::')
                scoreold = float(scoreold)
                if score > scoreold:
                    # better hit
                    score_line[spacer] = '%s::%s'%(score,lines)
                elif score == scoreold:
                    score_line[spacer]=('%s::%s;;%s' % (score, lines_old, lines))
                else:
                    # worse hit
                    pass
        else:
            emptyout.add(lines)
    for spacer in score_line:
        scoresetold = score_line[spacer]
        scoreold, lines_all = scoresetold.split('::')
        for lines in lines_all.split(';;'):
            if '\n' not in lines:
                newoutput.add(lines + '\n')
            else:
                newoutput.add(lines)
    f1 = open(filename + '.filter.txt', 'w')
    f1.write(''.join(list(newoutput)))
    f1.write(''.join(list(emptyout)))
    f1.close()

allfasta = glob.glob('%s/all.spacers.fasta.selfblast.txt.*.fasta.allanno.txt'%(output_folder)) + \
           glob.glob('%s/all.spacers.fasta.allanno.txt'%(output_folder))
for fasta in allfasta:
    fiter_score(fasta, 0,1)
    #check output
    allset = set()
    for lines in open(fasta,'r'):
        allset.add(lines.split('\t')[0])
    for lines in open(fasta + '.filter.txt','r'):
        if lines.split('\t')[0] not in allset:
            print('missing in .filter.txt', lines)

################################################### END ########################################################
################################################### SET PATH ########################################################
# spacers in a genome targeting the same phage (ID)
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
id_cutoff = 90
kmer_cutoff = 16
cluster_cutff = 0.9

def spacer_donor(spacer):
    genome = spacer.split('__')[0]
    donor = genome.split('_')[0]
    species = genome.split('_')[1]
    return [donor,species,genome]

def split_donor(spacer_list,target,output_list):
    allgenome = set()
    allspecies = set()
    for spacer in spacer_list:
        donor, species, genome = spacer_donor(spacer)
        allgenome.add(genome)
        allspecies.add(species)
    if len(allspecies) > 1:
        tag = 'diffgenome.diffspecies'
    elif len(allgenome) > 1:
        tag = 'diffgenome.samespecies'
    else:
        tag = 'samegenome'
    for spacer in spacer_list:
        donor, species, genome = spacer_donor(spacer)
        output_list.add('%s\t%s\t%s\t%s\t%s\t\n'%(target,tag,species, genome,spacer))
    return output_list

def find_same_donor(anno_file):
    same_donor = dict()
    output_list = set()
    for lines in open(anno_file,'r'):
        if 'None' not in lines:
            lines_set = lines.split('\t')
            spacer,evalue,target_set,dataset,nonuse = lines_set[0:5]
            target_set = target_set.split(';')
            for target in target_set:
                target = '%s\t%s'%(target,dataset)
                same_donor.setdefault(target,set())
                same_donor[target].add(spacer)
    for target in same_donor:
        spacer_list = same_donor[target]
        if len(spacer_list) > 1:
            output_list = split_donor(spacer_list,target,output_list)
    f1 = open('%s.target.txt'%(anno_file),'w')
    f1.write('target_phage\tdataset\ttag\tspecies\tgenome\tspacer\t\n')
    f1.write(''.join(list(output_list)))
    f1.close()

# spacer targeting the same phage by anno
allanno_file = '%s/all.spacers.fasta.allanno.txt.filter.txt'%(output_folder)
find_same_donor(allanno_file)

################################################### END ########################################################
################################################### SET PATH ########################################################
# add time tag and clonal pop tag on genomes
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary/allspacers_withtimetag/'
time_meta = '/scratch/users/anniz44/scripts/1MG/virus/table_BN10_WGS_metadata_07052018.txt'
genome_name_meta = '/scratch/users/anniz44/genomes/donor_species/selected_species/file.change.name.new.txt'
clonal_pop_files = '/scratch/users/anniz44/genomes/pan-genome/roary/clonal_population/*.genome.cluster.txt'

# load time
time_tag = dict()
for lines in open(time_meta,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    oldname = lines_set[0]
    time = lines_set[4]
    time_tag.setdefault(oldname,time)

#load clonal
clonal_pop = dict()
for files in glob.glob(clonal_pop_files):
    for lines in open(files,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        species,newname,cluster = lines_set[0:3]
        clonal_pop.setdefault(newname,'%s_%s'%(species,cluster))

# load all genomes
newgenome_info = dict()
newgenome_infoout = set()
for lines in open(genome_name_meta,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    oldname = lines_set[0]
    newname = lines_set[1]
    newtime = time_tag.get(oldname,'None')
    clonal_info = clonal_pop.get(newname,'None')
    longname = '%s__Time%s__Cluster%s'%(newname,newtime,clonal_info)
    newgenome_infoout.add('%s\t%s\t%s\t%s\t%s\t\n'%(oldname,newname,newtime,clonal_info,longname))
    newgenome_info.setdefault(newname,longname)
    if newtime == 'None':
        print('missing time tag %s %s'%(newname,oldname))

f1 = open(genome_name_meta + '.addtime.addclonal.txt','w')
f1.write('oldname\tnewname\ttime_tag\tclonal_pop\tlong_newname\t\n')
f1.write(''.join(list(newgenome_infoout)))
f1.close()

# add tags to all results
for fastafile in glob.glob(output_folder + '/../notimetag/*.fasta'):
    fastafilename = os.path.split(fastafile)[-1]
    newoutput = set()
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        genomename = record_id.split('__')[0]
        if genomename in newgenome_info:
            record_id = record_id.replace(genomename + '__',newgenome_info[genomename] + '__')
        else:
            print('missing info for %s'%(genomename))
        newoutput.add('>%s\n%s\n'%(record_id, str(record.seq)))
    f1 = open(os.path.join(output_folder,fastafilename),'w')
    f1.write(''.join(list(newoutput)))
    f1.close()

for testfile in glob.glob(output_folder + '/../notimetag/*.txt'):
    testfilename = os.path.split(testfile)[-1]
    newoutput = set()
    if 'allanno' in testfilename:
        colnum = -1
        for lines in open(testfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            if colnum == -1:
                for i in range(0,len(lines_set)):
                    if lines_set[i].split('__')[0] in newgenome_info:
                        colnum = i
                        newoutput.add('%s\t%s' % (newgenome_info[lines_set[i].split('__')[0]],lines))
                        break
                if colnum == -1:
                    newoutput.add('newgenome\t%s'%(lines))
                    print(lines)
            else:
                newoutput.add('%s\t%s' % (newgenome_info[lines_set[colnum].split('__')[0]], lines))
    else:
        for lines in open(testfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            newoutput.add('%s\t%s\t%s' % (newgenome_info[lines_set[0].split('__')[0]], newgenome_info[lines_set[1].split('__')[0]],lines))
    f1 = open(os.path.join(output_folder,testfilename),'w')
    f1.write(''.join(list(newoutput)))
    f1.close()

id_cutoff = 90
kmer_cutoff = 16

for testfile in glob.glob(output_folder + '/../details/all.spacers.fasta.selfblast.txt'):
    testfilename = os.path.split(testfile)[-1]
    newoutput = set()
    allspacer_set = set()
    for lines in open(testfile,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        spacer1 = lines_set[0]
        spacer2 = lines_set[1]
        if spacer1 != spacer2:
            spacer_set = [spacer1, spacer2]
            spacer_set.sort()
            spacer_set = ':'.join(spacer_set)
            if spacer_set not in allspacer_set:
                allspacer_set.add(spacer_set)
                if float(lines_set[2]) >= id_cutoff and int(lines_set[3]) >= kmer_cutoff:
                    newoutput.add('%s\t%s\t%s' % (newgenome_info[lines_set[0].split('__')[0]], newgenome_info[lines_set[1].split('__')[0]],lines))
    f1 = open(os.path.join(output_folder,testfilename),'w')
    f1.write(''.join(list(newoutput)))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# merge same spacers of strains of the same time point and/or same clonal pop
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
genome_name_meta = '/scratch/users/anniz44/genomes/donor_species/selected_species/file.change.name.new.txt.addtime.addclonal.txt'
# find all unique spacers
samespacer = dict()
for record in SeqIO.parse(os.path.join(output_folder,'allspacers_withtimetag/all.spacers.fasta'), 'fasta'):
    record_id = str(record.id)
    record_seq = str(record.seq)
    samespacer.setdefault(record_seq,set())
    samespacer[record_seq].add(record_id)
# merge and rename spacers of same clonal pop (CP)
spacer_rename = []
newspacer_count = dict()
selected_spacer = dict()
selected_spacer_fasta_rep = []
Seq_num = 0
for record_seq in samespacer:
    allgenome = samespacer[record_seq]
    Seq_num += 1
    allCP = dict()
    for spacer in allgenome:
        CP = '_'.join(spacer.split('_')[0:2]) + '__' + '__'.join(spacer.split('__G')[0].split('__')[1:])
        genomename = spacer.split('__')[0]
        if CP in allCP:
            spacer_newname = allCP[CP]
        else:
            newspacer_count.setdefault(CP, 0)
            newspacer_count[CP] += 1
            spacer_newname = 'SPG%s:%s__SP%s' % (Seq_num, CP, newspacer_count[CP])
            allCP.setdefault(CP, spacer_newname)
        spacer_rename.append('SPG%s\t%s\t%s\t%s\t%s\t%s\t\n' % (Seq_num,spacer_newname, spacer, genomename, CP, record_seq))
        selected_spacer.setdefault(spacer, ['SPG%s'%(Seq_num),spacer_newname])
    selected_spacer_fasta_rep.append('>SPG%s\n%s\n' % (Seq_num, record_seq))

f1 = open(os.path.join(output_folder, 'all.spacers.newname.txt'), 'w')
f1.write('spacergroup\tnewspacer\tspacer\tgenome\tCP\tseq\t\n')
f1.write(''.join(spacer_rename))
f1.close()

f1 = open(os.path.join(output_folder, 'all.spacers.fasta'), 'w')
f1.write(''.join(selected_spacer_fasta_rep))
f1.close()
# annotate each SPG (spacer group)
anno_out = set()
annotated = set()
for lines in open(os.path.join(output_folder,'allspacers_withtimetag/all.spacers.fasta.allanno.txt.filter.txt'),'r'):
    lines_set = lines.split('\n')[0].split('\t')
    if 'None' not in lines_set:
        oldname2, oldname1, bitscore, accession, dataset = lines_set[0:5]
        oldname = oldname2 + '__' + '__'.join(oldname1.split('__')[1:])
        SPG, spacer_newname = selected_spacer[oldname]
        anno_out.add('%s\t%s\t%s\t%s\t\n'%(SPG,bitscore,accession,dataset))
        annotated.add(SPG)

for i in range(1,Seq_num+1):
    SPG = 'SPG%s'%(i)
    if SPG not in annotated:
        anno_out.add('%s\tNone\t\t\t\n'%(SPG))

anno_out = list(anno_out)
anno_out.sort()
f1 = open(os.path.join(output_folder, 'all.spacers.fasta.allanno.txt.filter.txt'), 'w')
f1.write(''.join(anno_out))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# spacers mapping results
'/scratch/users/anniz44/scripts/1MG/virus/spacers_mapping.py'
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
genome_name_meta = '/scratch/users/anniz44/genomes/donor_species/selected_species/file.change.name.new.txt.addtime.addclonal.txt'

os.system('rm %s/*target.txt %s/*selfblast*'%(output_folder,output_folder))
# read metadata
selected_spacer = dict()
for lines in open('%s/all.spacers.newname.txt'%(output_folder)):
    lines_set = lines.split('\n')[0].split('\t')
    SPG, spacer_newname,oldname = lines_set[0:3]
    selected_spacer.setdefault(oldname,[SPG, spacer_newname])

# spacers mapping results
def anno_spacer(newname,oldname):
    newname = newname.split(':')[1]
    donor_species = newname.split('__')[0]
    CP = newname.split('__SP')[0]
    cluster = newname.split('__SP')[0].split('__Cluster')[1]
    genome = oldname.split('__Time')[0]
    return [donor_species,CP,cluster,genome]

def anno_spacer_pair(oldname1,oldname2):
    comset = []
    SPG1, spacer_newname1 = selected_spacer[oldname1]
    SPG2, spacer_newname2 = selected_spacer[oldname2]
    donor_species1, CP1, cluster1, genome1 = anno_spacer(spacer_newname1, oldname1)
    donor_species2, CP2, cluster2, genome2 = anno_spacer(spacer_newname2, oldname2)
    if donor_species1 == donor_species2:
        comset.append('samespecies')
        if cluster1 == cluster2:
            comset.append('samecluster')
            if CP1 == CP2:
                comset.append('sameCP')
                if genome1 == genome2:
                    comset.append('samegenome')
                else:
                    comset.append('diffgenome')
            else:
                comset.append('diffCP')
        else:
            comset.append('diffcluster')
    else:
        comset.append('diffspecies')
    comset_result = '.'.join(comset)
    return [comset_result,'\t'.join([spacer_newname1,spacer_newname2,genome1,genome2])]

def output_anno(allanno_out):
    for comset_result in allanno_out:
        anno_out = allanno_out[comset_result]
        anno_out = list(anno_out)
        anno_out.sort()
        f1 = open(os.path.join(output_folder, 'all.spacers.fasta.selfblast.txt.%s.txt' % (comset_result)), 'a')
        f1.write(''.join(anno_out))
        f1.close()

allanno_out = dict()
i = 0
for lines in open(os.path.join(output_folder,'allspacers_withtimetag/all.spacers.fasta.selfblast.txt'),'r'):
    lines_set = lines.split('\n')[0].split('\t')
    genome1, genome2, spacer1, spacer2 = lines_set[0:4]
    oldname1 = genome1 + '__' + '__'.join(spacer1.split('__')[1:])
    oldname2 = genome2 + '__' + '__'.join(spacer2.split('__')[1:])
    comset_result , anno_details = anno_spacer_pair(oldname1, oldname2)
    allanno_out.setdefault(comset_result,set())
    allanno_out[comset_result].add('%s\t%s\t%s\t%s\t\n' % (anno_details,
                                         lines_set[4],lines_set[5],lines_set[6]))
    i += 1
    if (i%1000) == 0:
        output_anno(allanno_out)
        allanno_out = dict()
        print('output %s lines'%(i))

# spacers same target
def anno_spacer(newname,oldname):
    newname = newname.split(':')[1]
    donor_species = newname.split('__')[0]
    CP = newname.split('__SP')[0]
    cluster = newname.split('__SP')[0].split('__Cluster')[1]
    genome = oldname.split('__Time')[0]
    return [donor_species,CP,cluster,genome]

def anno_spacer_set(allspacers,dataset,target,allanno_out):
    anno_set = dict()
    comset = []
    donor_speciesset = set()
    CPset = set()
    clusterset = set()
    genomeset = set()
    for oldname1 in allspacers:
        SPG1, spacer_newname1 = selected_spacer[oldname1]
        donor_species1, CP1, cluster1, genome1 = anno_spacer(spacer_newname1, oldname1)
        donor_speciesset.add(donor_species1)
        CPset.add(CP1)
        clusterset.add(cluster1)
        genomeset.add(genome1)
        anno_set.setdefault(oldname1,[spacer_newname1,genome1])
    if len(donor_speciesset)  == 1 :
        comset.append('samespecies')
        if len(clusterset)  == 1:
            comset.append('samecluster')
            if len(CPset)  == 1:
                comset.append('sameCP')
                if len(genomeset)  == 1:
                    comset.append('samegenome')
                else:
                    comset.append('diffgenome')
            else:
                comset.append('diffCP')
        else:
            comset.append('diffcluster')
    else:
        comset.append('diffspecies')
    comset_result = '.'.join(comset)
    allanno_out.setdefault(comset_result,set())
    for oldname1 in anno_set:
        allanno_out[comset_result].add('%s\t%s\t%s\t\n'%(target,'\t'.join(anno_set[oldname1]),
                                                                               dataset))
    return allanno_out

def output_anno(allanno_out):
    for comset_result in allanno_out:
        anno_out = allanno_out[comset_result]
        anno_out = list(anno_out)
        anno_out.sort()
        f1 = open(os.path.join(output_folder, 'all.spacers.fasta.allanno.txt.filter.txt.%s.target.txt' % (comset_result)), 'a')
        f1.write(''.join(anno_out))
        f1.close()

alltarget = dict()
for lines in open(os.path.join(output_folder, 'allspacers_withtimetag/all.spacers.fasta.allanno.txt.filter.txt'), 'r'):
    if 'None' not in lines:
        lines_set = lines.split('\n')[0].split('\t')
        genome1, spacer1 = lines_set[0:2]
        oldname1 = genome1 + '__' + '__'.join(spacer1.split('__')[1:])
        alltargetset = lines_set[3].split(';')
        dataset = lines_set[4]
        for target in alltargetset:
            if target!= '':
                alltarget.setdefault(target,[set(),dataset])
                alltarget[target][0].add(oldname1)

allanno_out = dict()
i = 0
for target in alltarget:
    allspacers,dataset = alltarget[target]
    if len(allspacers) > 1:
        allanno_out = anno_spacer_set(allspacers,dataset,target,allanno_out)
    i += 1
    if (i%1000) == 0:
        output_anno(allanno_out)
        allanno_out = dict()
        print('output %s lines'%(i))

# remove duplicated strains in the same CP
os.system('mv %s/*sameCP.diffgenome* %s/details'%(output_folder,output_folder))

################################################### END ########################################################
################################################### SET PATH ########################################################
# spacers dynamics in each species
'/scratch/users/anniz44/scripts/1MG/virus/summarize_spacer.py'
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

Tag = {'diffspecies':1,
       'samespecies.diffcluster':2,
       'samespecies.samecluster.diffCP':3,
       'samespecies.samecluster.sameCP.samegenome':4}
output_folder = '/scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/summary//'
# summarize all unique spacers in species
def anno_newspacer(newspacer):
    newspacer =newspacer.split(':')[1]
    donor_species = newspacer.split('__')[0]
    CP = newspacer.split('__SP')[0]
    cluster = CP.split('__Cluster')[1]
    return [donor_species,CP,cluster]

allspacer = dict()
allCP = dict()
CP_cluster = dict()
for lines in open('%s/all.spacers.newname.txt'%(output_folder),'r'):
    if not lines.startswith('spacergroup'):
        lines_set = lines.split('\n')[0].split('\t')
        newspacer = lines_set[1]
        donor_species, CP, cluster = anno_newspacer(newspacer)
        allspacer.setdefault(donor_species,set())
        allspacer[donor_species].add(CP)
        allCP.setdefault(CP,[set(),set(),set(),set(),set()])
        allCP[CP][0].add(newspacer)
        CP_cluster.setdefault(CP,cluster)

# summarize all spacers in CP
def anno_newspacer(newspacer):
    newspacer =newspacer.split(':')[1]
    donor_species = newspacer.split('__')[0]
    CP = newspacer.split('__SP')[0]
    cluster = CP.split('__Cluster')[1]
    return [donor_species,CP,cluster]
alloutput = set()
for files in glob.glob('%s/all.spacers.fasta.selfblast.txt.*'%(output_folder)):
    tag = os.path.split(files)[-1].split('all.spacers.fasta.selfblast.txt.')[1].split('.txt')[0]
    allCPpair = dict()
    for lines in open(files,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        newspacer1 = lines_set[0]
        newspacer2 = lines_set[1]
        donor_species1, CP1, cluster1 = anno_newspacer(newspacer1)
        donor_species2, CP2, cluster2 = anno_newspacer(newspacer2)
        CPpair1 = '%s\t%s'%(CP1,CP2)
        CPpair2 = '%s\t%s'%(CP2,CP1)
        allCPpair.setdefault(CPpair1,set())
        allCPpair[CPpair1].add(newspacer1)
        allCP[CP1][Tag[tag]].add(newspacer1)
        allCP[CP2][Tag[tag]].add(newspacer2)
        if CPpair1!=CPpair2:
            allCPpair.setdefault(CPpair2, set())
            allCPpair[CPpair2].add(newspacer2)
    for CPpair in allCPpair:
        alloutput.add('%s\t%s\t%s\t\n'%(tag,CPpair,len(allCPpair[CPpair])))

f1 = open(os.path.join(output_folder, 'all.CP.selfblast.sum'), 'w')
f1.write('Tag\tCP1\tCP2\tNo.spacers\t\n')
f1.write(''.join(list(alloutput)))
f1.close()


alloutput = []
for donor_species in allspacer:
    for CP in allspacer[donor_species]:
        alloutput.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,CP_cluster[CP],CP,
                                               len(allCP[CP][0]),len(allCP[CP][1]),len(allCP[CP][2]),len(allCP[CP][3]),len(allCP[CP][4]),
                                                                   len([i for i in allCP[CP][0]
                                                                        if i not in allCP[CP][1] and i not in allCP[CP][2]
                                                                        and i not in allCP[CP][3] and i not in allCP[CP][4]
                                                                        ])
                                               ))

f1 = open(os.path.join(output_folder, 'all.donor.species.sum'), 'w')
f1.write('donor_species\tcluster\tCP\tNo.allspacers\tNo.diffspecies\tNo.samespecies.diffcluster\tNo.samespecies.samecluster.diffCP\tNo.samespecies.samecluster.sameCP.samegenome\tNo.uniquespacers\t\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# notused
'wget --post-file="/scratch/users/anniz44/genomes/GHM/newgenomes/GMC_Genome_Library_01312020/0262AD_0917_052_A3_final.scaffolds.fasta" "http://phaster.ca/phaster_api?contigs=1" -O /scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/prophage/0262AD_0917_052_A3_final.scaffolds.fasta.prophage.contig1.txt'
'wget --post-file="/scratch/users/anniz44/genomes/GHM/newgenomes/GMC_Genome_Library_01312020/0262AD_0917_052_A3_final.scaffolds.fasta" "http://phaster.ca/phaster_api" -O /scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/prophage/0262AD_0917_052_A3_final.scaffolds.fasta.prophage.txt'
'wget "http://phaster.ca/phaster_api?acc=ZZ_2ab1b1f8b4" -O /scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/prophage/0262AD_0917_052_A3_final.scaffolds.fasta.prophage.contig1.txt'
'wget "http://phaster.ca/phaster_api?acc=ZZ_d72aace991" -O /scratch/users/anniz44/genomes/GHM/newgenomes/crispr_MG/prophage/0262AD_0917_052_A3_final.scaffolds.fasta.prophage.txt'
# blastx on MAC
'top | grep 91731'



