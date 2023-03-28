import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

refdir = '/scratch/users/anniz44/Metagenomes/BN10_MG/assembly/derep/'
genomedir = '/scratch/users/anniz44/Metagenomes/BN10_MG/assembly/bins/'
scriptdir = '/scratch/users/anniz44/scripts/1MG/crispr/fastANI/'
outdirANI = '/scratch/users/anniz44/Metagenomes/BN10_MG/assembly/bins_derep_map/'

def donor_genome(genomelist):
    donor_dict = dict()
    for genome in genomelist:
        donor = os.path.basename(genome)[:2]
        donor_dict.setdefault(donor,[])
        donor_dict[donor].append(genome)
    return donor_dict

def outputlist(genomelist,outname):
    f1 = open(outname,'w')
    f1.write('\n'.join(genomelist))
    f1.close()

def outputscript(reflistout,genomelistout,scriptout,aniout):
    f1 = open(scriptout, 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
    f1.write('fastANI --ql %s --rl %s -o %s'%(genomelistout,reflistout,aniout))
    f1.close()

allgenome = glob.glob('%s/*.fa'%(genomedir))
allref = glob.glob('%s/*.fa'%(refdir))

allgenomedonor = donor_genome(allgenome)
allrefdonor = donor_genome(allref)

for donor in allgenomedonor:
    genomelist = allgenomedonor[donor]
    reflist = allrefdonor.get(donor,[])
    if len(reflist) > 0:
        # output genome list
        reflistout = '%s/%s.ref'%(scriptdir,donor)
        genomelistout = '%s/%s.genome'%(scriptdir,donor)
        outputlist(genomelist, genomelistout)
        outputlist(reflist, reflistout)
        # scripts to run fastANI
        scriptout = '%s/%s.sh'%(scriptdir,donor)
        aniout  = '%s/%s.out'%(outdirANI,donor)
        outputscript(reflistout, genomelistout, scriptout, aniout)
