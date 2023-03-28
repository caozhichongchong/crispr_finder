import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/crispr_MG/anno/allspacers_withtimetag/'
time_meta = '/scratch/users/anniz44/scripts/1MG/virus/table_BN10_WGS_metadata_07052018.txt'
genome_name_meta = '/scratch/users/anniz44/genomes/donor_species/selected_species/file.change.name.new.txt'
clonal_pop_files = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/*.genome.cluster.txt'

def output_fasta(outputlist,outputfilename):
    f1 = open(outputfilename,'w')
    f1.write(''.join(outputlist))
    f1.close()

def load_clonal(genome_name_meta):
    clonal = dict()
    clonal_new = []
    for lines in open(genome_name_meta,'r'):
        if not lines.startswith('oldname'):
            lines_set = lines.split('\n')[0].split('\t')
            newname = lines_set[4].split('_')
            newname2 = '%s_%s_%s_%s'%(newname[0], newname[1],'_'.join(newname[3:]),newname[2])
            clonal.setdefault(lines_set[1],newname2)
            clonal_new.append('%s\t%s'%(newname2,lines))
    output_fasta(clonal_new, genome_name_meta + '.newname.txt')
    return clonal

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

# add clonal name
clonal = load_clonal(genome_name_meta + '.addtime.addclonal.txt')