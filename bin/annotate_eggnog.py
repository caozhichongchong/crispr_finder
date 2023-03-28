import os
import glob
import copy
from Bio.Seq import Seq

database_metacyc = '/scratch/users/mit_alm/database/metacyc/genes.col'
database_eggnog = '/scratch/users/mit_alm/database/eggnog/2_annotations.tsv'
database_eggnog2 = '/scratch/users/mit_alm/database/eggnog/COG_catogary.txt'
database_kegg = '/scratch/users/anniz44/scripts/database/Kegg/ko_formated'
output_file = '/scratch/users/anniz44/genomes/crispr_MG/summary'
all_fasta_set = glob.glob(os.path.join(output_file, 'all.BL.widespreadspacer.neighbour.fasta.faa'))

# function annotation
def best_hit(Blast_output,small = 0):
    for gene_name in Blast_output:
        db_name,identity = Blast_output[gene_name]
        if len(identity) > 2:
            identity2 = copy.deepcopy(identity)
            if small == 1:
                identity2.sort()
            else:
                identity2.sort(reverse=True)
            top1 = identity2[0]
            # check the top 4 hits, keep COG annotation
            try:
                if not db_name[identity.index(top1)].startswith('COG'):
                    top1 = identity2[1]
                    if not db_name[identity.index(top1)].startswith('COG'):
                        top1 = identity2[2]
                        if not db_name[identity.index(top1)].startswith('COG'):
                            top1 = identity2[3]
                            if not db_name[identity.index(top1)].startswith('COG'):
                                top1 = identity2[0]
            except IndexError:
                pass
            Blast_output[gene_name]=[db_name[identity.index(top1)]#,
                                     #db_name[identity.index(top2)]
             ]
        else:
            Blast_output[gene_name]=db_name
    return Blast_output

def annotate_eggnog(blast_search):
    Anno_eggnog = dict()
    for lines in open(database_eggnog2,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        Anno_eggnog.setdefault(lines_set[0],lines.split('\n')[0])
    Blast_output = dict()
    DB_name = dict()
    DB_name_2 = dict()
    try:
        f1 = open(all_fasta + '.eggnog.all.txt','r')
    except IOError:
        os.system('cat %s > %s'%(' '.join(blast_search),all_fasta + '.eggnog.all.txt'))
    blast_search_file = all_fasta + '.eggnog.all.txt'
    for lines in open(blast_search_file, 'r'):
        if not lines.startswith('#'):
            db_name = ''
            identity = 0
            lines_set = lines.split('\n')[0].split(' ')
            gene_name = lines_set[0]
            for sub_line in lines_set[1:]:
                if sub_line != '' and sub_line != '-':
                    if db_name == '':
                        db_name = sub_line.split('.')[0]
                    elif identity == 0:
                        identity = float(sub_line)
                        break
            Blast_output.setdefault(gene_name, [[], []])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['', ''])
            DB_name_2.setdefault(db_name, [])
    for database in [database_eggnog]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[1]
                function_name = ''
                annotation_catogary = []
                for COG in lines_set[2]:
                    annotation_catogary.append(Anno_eggnog[COG])
                annotation_fun = lines_set[3]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
                    for annotation_catogary_sub in annotation_catogary:
                        DB_name_2[db_name].append('%s%s\t%s'%(function_name,annotation_fun,annotation_catogary_sub))
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name,DB_name_2]

def sum_kegg_eggnog(Blast_outputegg,DB_nameegg):
    Output = set()
    for gene_name in Blast_outputegg:
        for db_name in Blast_outputegg[gene_name]:
            if db_name in DB_nameegg:
                for funegg in DB_nameegg[db_name]:
                    Output.add('%s\t%s\t%s\n' % (
                        gene_name, db_name, funegg))
    f1 = open(all_fasta + '.all.eggnog.sum', 'w')
    f1.write(
        'gene_name\tEGGNOG\tannotation\tCOG\tCOG1\tCOG2\n' + ''.join(
            list(Output)))
    f1.close()

for all_fasta in all_fasta_set:
    all_fasta_folder, all_fasta_file = os.path.split(all_fasta)
    # sum up
    Blast_output3, DB_name3, DB_name3_2 = annotate_eggnog(glob.glob(all_fasta + '.eggnog.*.txt'))
    sum_kegg_eggnog(Blast_output3, DB_name3_2)
