import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import json

input_folder = '/scratch/users/anniz44/genomes/crispr_MG/CRISPRCasFinder/'
genome_foler = '/scratch/users/mit_alm/IBD_Evo/BL/Assembly_for_gene_flow/'
output_folder = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
filename_format = '_final'
crispr_level_cutoff = 2  # only level 2-4 + with cas
spacer_target = '/scratch/users/anniz44/genomes/crispr_MG/summary/BL.widespread.spacer.fasta'
species_target = 'BL'
extract_dis = 1000 # extract X bp
align_by_donor = True # cluster by 90% similarity

def load_target():
    SPall = []
    for record in SeqIO.parse(spacer_target, 'fasta'):
        SPall.append(str(record.seq))
    return SPall

def changepos(posold,posnew):
    return[min(posold[0],posnew[0]),max(posold[1],posnew[1])]

def parse_crispr(SPall,SPfile):
    allSPpos = dict()
    # only crisprcasfinder
    with open(SPfile) as f:
        data = json.load(f)
    for subseq in data['Sequences']:
        for crispr in subseq['Crisprs']:
            if crispr['Evidence_Level'] >= crispr_level_cutoff and subseq['Cas'] != []:
                contig_ID = '_'.join(crispr['Name'].split('_')[:-2])
                allSPpos.setdefault(contig_ID, [])
                allSPpos[contig_ID].append([])
                for spacer in crispr['Regions']:
                    if spacer['Type'] == 'Spacer':
                        record_seq = spacer['Sequence']
                        # complementary reverse + sort
                        record_seq2 = str(Seq(record_seq).reverse_complement())
                        record_list = [record_seq, record_seq2]
                        record_list.sort()
                        record_seq = record_list[0]
                        if record_seq in SPall:
                            pos1 = int(spacer['Start'])
                            pos2 = int(spacer['End'])
                            if allSPpos[contig_ID][-1] == []:
                                allSPpos[contig_ID][-1] = [pos1,pos2]
                            else:
                                allSPpos[contig_ID][-1] = changepos(allSPpos[contig_ID][-1],[pos1,pos2])
    genomename = os.path.split(SPfile)[-1].split('_final.scaffolds.fasta.result.json')[0]
    reffile = glob.glob('%s/%s/scaffolds.fasta'%(genome_foler,genomename))
    allSPoutput = []
    print(allSPpos)
    if allSPpos != {}:
        if reffile == []:
            print('no fasta for %s'%(SPfile))
        else:
            for record in SeqIO.parse(reffile[0], 'fasta'):
                contig_ID = str(record.id)
                contig_ID = '_'.join(contig_ID.split('_')[:-1])
                if allSPpos.get(contig_ID,[[]])!= [[]]:
                    contigseq = str(record.seq)
                    contigseqlen = len(contigseq)
                    contigposall = allSPpos[contig_ID]
                    for contigpos in contigposall:
                        if contigpos != []:
                            #print(contigseqlen,max(contigpos[0]-1 - extract_dis,0),min(contigpos[1] + extract_dis,contigseqlen))
                            allSPoutput.append('>%s_%s_%s_%s\n%s\n'%(genomename,contig_ID,contigpos[0],contigpos[1],
                                                                contigseq[max(contigpos[0]-1 - extract_dis,0):min(contigpos[1] + extract_dis,contigseqlen)]))
            f1.write(''.join(allSPoutput))

def fasta_separate_by_donor(input_file):
    HGTseqall = dict()
    for record in SeqIO.parse(input_file, 'fasta'):
        record_id = str(record.id)
        donor = record_id.split('_')[0]
        HGTseqall.setdefault(donor,'')
        HGTseqall[donor] += '>%s\n%s\n'%(record_id,str(record.seq))
    for donor in HGTseqall:
        outputfilename = '%s/multidonor/all.BL.widespreadspacer.neighbour.%s.fasta' % (output_folder,donor)
        f1 = open(outputfilename, 'w')
        f1.write(HGTseqall[donor])
        f1.close()
        #cutoff = 0.9
        #os.system(('usearch -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
        #               % (outputfilename, cutoff, outputfilename,
        #                  outputfilename, 8)))


f1 = open('%s/all.BL.widespreadspacer.neighbour.fasta' % (output_folder), 'w')
SPall = load_target()
allSPfiles = glob.glob('%s/*_%s*' % (input_folder, species_target))
for SPfile in allSPfiles:
    print('process %s' % (SPfile))
    parse_crispr(SPall, SPfile)
f1.close()
try:
    os.mkdir('%s/multidonor/' % (output_folder))
except IOError:
    pass
fasta_separate_by_donor('%s/all.BL.widespreadspacer.neighbour.fasta' % (output_folder))

# run mafft_BLHGT.sh
os.system('cat %s/multidonor/*.cluster.aa > %s/all.BL.widespreadspacer.neighbour.fasta'%(output_folder,output_folder))
#os.system('maffttree_run %s/all.BL.widespreadspacer.neighbour.fasta 8 fast'%(output_folder))