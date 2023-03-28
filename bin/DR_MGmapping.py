import glob,os
from Bio import SeqIO
SearchMG = False
input_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/DR/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/DR/MG'
MG_DR = '%s/MG_DR.fasta'%(output_folder)
WGS_DR = glob.glob('%s/*.DR.fasta'%(input_folder))

if SearchMG:
    #make MG DR db
    os.system('usearch -makeudb_usearch %s -output %s.udb'%(MG_DR,MG_DR))

    def run_usearch(MG_DR,fasta):
        outputfile = os.path.join(output_folder,os.path.split(fasta)[-1].replace('.fasta','.MG.txt'))
        os.system('usearch -usearch_global %s -db %s.udb -strand both -id 0.9 -blast6out %s -threads 40 -maxaccepts 0 -maxrejects 0'%(
            fasta,MG_DR,outputfile
        ))

    for fasta in WGS_DR:
        run_usearch(MG_DR, fasta)

os.system('find %s/ -type f -empty -delete'%(output_folder))
# sum result
def load_MG(MG_DR):
    allMGDR = set()
    for record in SeqIO.parse(MG_DR, 'fasta'):
        allMGDR.add(str(record.id))
    return allMGDR

def compareblast(result1,result2):
    if result2[1]<result1[1]:
        return result2
    return result1

def load_blast(blastfile):
    donor_species = os.path.split(blastfile)[-1].split('.')[0]
    species = donor_species.split('_')[1].replace('PaDi','PB')
    for lines in open(blastfile,'r'):
        lines_set = lines.split('\t')
        MGDR = lines_set[1]
        mismatch = lines_set[4]
        if MGDR not in allMGDRresult:
            allMGDRresult.setdefault(MGDR,[species,mismatch])
        else:
            allMGDRresult[MGDR] = compareblast(allMGDRresult[MGDR],[species,mismatch])

alloutput = glob.glob('%s/*DR.MG.txt'%(output_folder))
# load MG DR
allMGDR = load_MG(MG_DR)
# load blast result
allMGDRresult = dict()
for blastfile in alloutput:
    load_blast(blastfile)
# sum blast result
alloutput = ['DR\tspecies\tmismatch\n']
alloutputempty = []
for DR in allMGDR:
    if DR in allMGDRresult:
        alloutput.append('%s\t%s\t%s\n'%(DR,allMGDRresult[DR][0],allMGDRresult[DR][1]))
    else:
        alloutputempty.append('%s\t\t\n'%(DR))

f1 = open('%s/MG_WGS.mapping.txt'%(output_folder),'w')
f1.write(''.join(alloutput))
f1.write(''.join(alloutputempty))
f1.close()
