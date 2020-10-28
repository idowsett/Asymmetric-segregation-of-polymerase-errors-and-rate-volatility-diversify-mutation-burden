import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

from matplotlib import collections  as mc
import pylab as pl


figsize = (16, 6)
fig, ax = plt.subplots(figsize=figsize)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.yticks(range (1,17))
plt.xlim(0, 1600000)
plt.ylim(0, 17)
ax.set_ylabel('Chromosome', fontname="Arial", fontsize=12)
ax.set_xlabel('Basepair Position', fontname="Arial", fontsize=12)
ax.set_title('Mutations Scored Over the Genome', fontname="Arial", fontsize=16)

# Blank Chromosomes Construction
lengths = [230218,813184,	316620,	1531933,576874,	270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066]
output = []
for i in range(0, 16):
    output.append([  (0, i+1), (lengths[i], i+1)    ])
c = "blue"
lc = mc.LineCollection([i for i in output], colors=c
, linewidths=2)
ax.add_collection(lc)


# Chromosome key lookup (For mutations and Masked genome construction)
dchrom = {'NC_001133' : 1, 'NC_001134' : 2, 'NC_001135' : 3, 'NC_001136' : 4, 'NC_001137' : 5 , 'NC_001138' : 6,
'NC_001139' : 7, 'NC_001140' : 8, 'NC_001141' : 9, 'NC_001142' : 10, 'NC_001143' : 11, 'NC_001144' : 12,
'NC_001145' : 13, 'NC_001146' : 14, 'NC_001147' : 15, 'NC_001148' : 16, 'NC_001224' : 'M'}



# Generate a means of looking up a position along a chromosome to see if it's masked ('N'), requires the referenced file: S288C-masked-genome.fasta
# This section builds the yeast masked genome positions ('N') as a list (N_list) in the form (Chromosome, position) of all masked positions
genome_dict, record_dict = {}, {}
for record in SeqIO.parse("S288C-masked-genome.fasta", "fasta"):
    record_dict[record.id] = record.seq
position = 0
N_list = []
for key,value in record_dict.items():
    position = 0
    for i in value:
        if i == 'N' and key != 'Mito':
            N_list.append((dchrom[key.split('|')[1]], position))
        position += 1
color = 'grey'
mark = '|'
for i in range (1,17):
    ax.scatter([x[1] for x in N_list if x[0] == i], ([i for x in N_list if x[0] == i]), s=None, c=color, marker=mark, label="Masked Base" if i ==1 else "", alpha=0.2,zorder=2)


# Mutation dataframe construction and graphing
def ExtractPosition(mutations):
    Position    =  int( mutations.split(':')[1])
    Chromosome  =   mutations.split('|')[1]
    return (Position)
# Generate a list of mutations:
mutation_file = 'FullMutationList.txt'
with open(mutation_file, 'r') as rr:
    mutation_lines  = [i.split('\t') for i in rr]
Chromosome      = [dchrom[i[3].split('|')[1]] for i in mutation_lines]
mutations = pd.DataFrame(mutation_lines, columns = ['Lineage', 'SubDivion', 'Colonies', 'Position', 'Segregant', 'DivisionID', 'Mutation', 'GenomicPosition'])
mutations['Position'] = mutations['Position'].apply(ExtractPosition)
mutations.insert( 4, 'chrom', Chromosome)
color = 'green'
mark = '|'
for i in range (1,17):
    ax.scatter([x.Position for x in mutations.itertuples() if x.chrom == i], ([i for x in mutations.itertuples() if x.chrom == i]), s=None, c=color, marker=mark, label="Mutation" if i ==1 else "", alpha=0.8,zorder=3)





# Optional Centromere addition:
cent_file = 'centromeres.txt'
with open(cent_file, 'r') as cent:
    cent_lines = [i.split('\t') for i in cent]
cents = [(int(i[2].split('..')[0])+int(i[2].split('..')[1]))/2 for i in cent_lines]
for i in range (1,17):
    ax.scatter(cents[i-1],i, s=20, c='red', marker='o', label="Centromere" if i ==1 else "", alpha=0.8, zorder=4)

ax.legend()
plt.show()
