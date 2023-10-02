# the start
#/scratch/users/anniz44/scripts/1MG/crispr
# run cirsprcasefinder runCRISPRCasFinder.py
cmds = 'python runCRISPRCasFinder.py'
# merge and summarize crispr spacers + DR merge_spacer_filter.py
cmds = 'python merge_spacer_filter.py'
# filter crisprcasfinder
# parsimony crispr
cmds = 'python spacer_parsimony.py'
cmds = 'python spacer_parsimony_sum.py '
cmds = 'python spacer_parsimony.py -pgain 1.30102999566 -pnochange 0.04575749056'
# annotate spacer
cmds = 'python annotate_spacer.py'
# structure of target spacers
cmds = 'python spacer_structure.py'
# phage seq from MG # not used
cmds = 'python MGphage.py'
# extract DR and SP of specific strains
cmds = 'python DR_extract.py'
# sum DR
cmds = 'python sumDR.py'
# map MG DR to WGS DR
cmds = 'python DR_MGmapping.py'
# remove contamination by DR
cmds = 'python DR_contamination.py'
# repeat parsimony crispr
cmds = 'python spacer_parsimony.py -pgain 1.30102999566 -pnochange 0.04575749056'


