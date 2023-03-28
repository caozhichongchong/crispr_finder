import glob,os
def bin(fasta):
    cmds = ''
    faname = os.path.basename(fasta).split(
        '.fasta')[0]
    try:
        f1 = open('/media/user/T7/AnniZ/assembly/bin/%s.1.fa'%(faname),'r')
    except IOError:
        cmds = 'metabat2 -i %s -o /media/user/T7/AnniZ/assembly/bin/%s\n'%(fasta,faname)
    return cmds

def derepfas(donor):
    cmds = 'dRep dereplicate /media/user/T7/AnniZ/assembly/derep/%s --debug -g /media/user/T7/AnniZ/assembly/bin/%s*.fa\n'%(donor,donor)
    return cmds

alldonor = dict()
for files in glob.glob('/media/user/T7/AnniZ/assembly/fasta/*.fasta'):
    donor = os.path.basename(files)[:2]
    alldonor.setdefault(donor,[])
    alldonor[donor].append(files)

print(alldonor)
for donor in alldonor:
    if len(alldonor[donor]) > 1:
        f1 = open('subscripts/%s.sh'%(donor),'w')
        cmds = '#!/bin/bash\nsource ~/.bashrc\nsource activate meta\n'
        for fasta in alldonor[donor]:
            cmds += bin(fasta)
        cmds += derepfas(donor)
        f1.write(cmds)
        f1.close()    

f1 = open('allbash.sh','w')
f1.write('#!/bin/bash\n')
for bashfile in glob.glob('subscripts/*.sh'):
    #f1.write('nohup %s > %s.out &\n'%(bashfile,bashfile))
    f1.write('%s > %s.out\n'%(bashfile, bashfile))
f1.close()
