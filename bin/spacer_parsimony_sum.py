# sum spacer_parsimony.py and set pgain and pnochang -> spacer_parsimony_sum.py
import os, glob
spacer_dir = '/scratch/users/anniz44/genomes/crispr_MG/summary/'
spacer_all_file = glob.glob('%s/*.pgain.txt'%(spacer_dir))
print(len(spacer_all_file))
spacer_all = dict()
for spacer_file in spacer_all_file:
    spacer_file_set = set()
    for lines in open(spacer_file,'r'):
        if not lines.startswith('Spacer'):
            spacer,pnochange,pgain,score1,score2 = lines.split('\n')[0].split('\t')
            pnochange = float(pnochange)
            pgain = float(pgain)
            pnochange_pgain = '%s\t%s'%(pnochange,pgain)
            if pnochange_pgain not in spacer_file_set:
                spacer_file_set.add(pnochange_pgain)
                spacer_all.setdefault(pnochange_pgain,0)
                spacer_all[pnochange_pgain] += float(score2)

alloutput = []
alloutput.append('Pnochange\tPgain\ttotal_score\n')
for pnochange_pgain in spacer_all:
    alloutput.append('%s\t%s\n'%(pnochange_pgain,spacer_all[pnochange_pgain]))

f1 = open('%s/all.Pgain.sum.txt'%(spacer_dir),'w')
f1.write(''.join(alloutput))
f1.close()
