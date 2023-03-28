import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json
# run spacer_BLHGT_extract.py and mafft_BLHGT.sh
input_file = '/scratch/users/anniz44/genomes/crispr_MG/summary/HGT/all.BL.widespreadspacer.neighbour.fasta.align'
input_file2 = '/scratch/users/anniz44/genomes/crispr_MG/summary/HGT/all.BL.widespreadspacer.neighbour.fasta'

def load_originaltarget():
    HGTseqall = dict()
    alllen = set()
    for record in SeqIO.parse(input_file2, 'fasta'):
        lenseq = len(str(record.seq))
        HGTseqall.setdefault(lenseq,[])
        HGTseqall[lenseq].append(str(record.id))
        alllen.add(lenseq)
    alllen = list(alllen)
    alllen.sort(reverse=True)
    for i in range(0,5):
        print('max seq %s %s'%(alllen[i],HGTseqall[alllen[i]]))

def load_target():
    print('start loading multi seq align %s'%(input_file))
    HGTseqall = []
    for record in SeqIO.parse(input_file, 'fasta'):
        HGTseqall.append(str(record.seq))
    return HGTseqall

def computeSNP(allseq,i):
    allalleles_count = {'a':0,'t':0,'g':0,'c':0}
    allalleles = [x[i] for x in allseq]
    allcount = []
    for allele in allalleles_count:
        allelecount = allalleles.count(allele)
        allalleles_count[allele] = allelecount
        allcount.append(allelecount)
    allelecountN = allalleles.count('n')
    allcount.sort(reverse=True)
    totalcount = sum(allcount) + allelecountN
    if totalcount == 0:
        print(allcount,allelecountN,allalleles)
    # return how manys SNPs differs from consensus (remove '-')
    return '%s\t%s\t%s\n'%((allcount[0] + allelecountN)/totalcount,totalcount,len([x for x in allcount if x!=0])-1)

def compare_across_HGT(HGTseqall):
    total_length = len(HGTseqall[0])
    print('start processing multi seq align for %s loci'%(total_length))
    alloutput = ['locus\tMAF\tcoverage\tSNPs\n']
    for i in range(0,total_length):
        alloutput.append('%s\t%s'%(i,computeSNP(HGTseqall,i)))
        if i%10000 == 0:
            print('finished %s loci'%(i))
    f1 = open(input_file + '.SNPcout.txt','w')
    f1.write(''.join(alloutput))
    f1.close()

load_originaltarget()
HGTseqall = load_target()
compare_across_HGT(HGTseqall)
