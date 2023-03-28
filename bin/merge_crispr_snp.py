import os,glob

output_folder_SNP = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
output_folder_crispr = glob.glob('/scratch/users/anniz44/genomes/crispr_MG/summary/')
SNPdifffile = '%s/all.SNP.pair.sum'%(output_folder_SNP)

def output_fasta(outputlist, outputfilename):
    f1 = open(outputfilename, 'a')
    f1.write(''.join(outputlist))
    f1.close()

for output_folder in output_folder_crispr:
    crisprdifffile = '%s/all.crispr.pair.sum'%(output_folder)
    alloutput = '%s/all.SNP.crispr.pair.sum'%(output_folder)
    f1 = open(alloutput,'w')
    f1.write('donor_species\tgenome1\tgenome2\tcrispr\tSNPs\tCP\n')
    f1.close()
    crispr = dict()
    for lines in open(crisprdifffile,'r'):
        genome1,genome2,crisprdiff = lines.split('\n')[0].split('\t')[:3]
        # order genomes
        genomelist = [genome1,genome2]
        genomelist.sort()
        crispr.setdefault('\t'.join(genomelist),[crisprdiff,'None','None'])
    for lines in open(SNPdifffile,'r'):
        CP,genome1,genome2,SNPdiff = lines.split('\n')[0].split('\t')[:4]
        # order genomes
        genomelist = [genome1,genome2]
        genomelist.sort()
        try:
            crispr['\t'.join(genomelist)][1] = SNPdiff
            crispr['\t'.join(genomelist)][2] = CP
        except KeyError:
            pass
    crispr_output = []
    for genomelist in crispr:
        crisprdiff,SNPdiff,CP=crispr[genomelist]
        try:
            crispr_output.append('%s\t%s\t%s\t%s\t%s\n'%(genomelist.split('_')[1],genomelist,crisprdiff,SNPdiff,CP))
        except IndexError:
            print(genomelist)
    output_fasta(crispr_output, alloutput)
