# For Fig.3(h-j), Dowsett et al.
import glob
import pandas as pd
import numpy as np
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import poisson
from scipy.stats import nbinom
from scipy.stats import gamma
import random
import math
from statistics import mean
from statistics import stdev

# Introduction
"""The following Python 3 code simulates the extent to which asymmetric segregation of
mutations lead to overdispersion in human ultramutator cancers assuming that a single Poisson
process governs the mutator phenotype in all cell divisions. To illustrate this concept,
we extrapolated the haploid genome mutation rate observed in yeast to humans."""

# Output Files(s) folder:
from pathlib import Path
File_save_location  = Path.cwd()

# Defined functions
"""chr_mu_dict is a dictionary where the keys represent chromosomes in a diploid
cell and the values, their associated mutation rates. This assumes each chr
contributes proportionately to the genome-wide mutation rate."""
def mk_chr_mu_dict(mu,chr_lenlist,tot_bp,chr_list_a,chr_list_b):
    ch_mu = []
    for i in chr_lenlist:
        ch_mu.append(mu*(i/tot_bp))
    tuplesA=list(zip(chr_list_a,ch_mu))
    tuplesB=list(zip(chr_list_b,ch_mu))
    chr_mu_dict = {}
    for i in tuplesA:
        chr_mu_dict[i[0]]=i[1]
    for i in tuplesB:
        chr_mu_dict[i[0]]=i[1]
    return chr_mu_dict

# creates dict where keys are chrs, and values are arrays (length c) of fixed muts from a Poisson.
def chr_mut_p(chr_mu_dict,c):
    return {k:(poisson.rvs(mu=v, size=c)) for (k,v) in chr_mu_dict.items()}

# as above, except rate is doubled to represent actual rate of mispair formation.
def chr_mut_bp(chr_mu_dict,c):
    return {k:(poisson.rvs(mu=(v*2), size=c)) for (k,v) in chr_mu_dict.items()}

# multiplies value by randomly selected 0 or 1 to mimic Mendelian segregation.
def binomialize(x):
    return x * (random.randint(0,1))

# returns index of dispersion for a distribution
def index_dis(df):
    return ((df["sum"].std())**2)/df["sum"].mean()

# simulates c iterations of fixed mutations counts arising by a Poisson process.
def mk_df_p(chr_mu_dict,c):
    data_poisson_p = chr_mut_p(chr_mu_dict,c)
    df_p = pd.DataFrame.from_dict(data_poisson_p)
    df_p["sum"]=df_p.sum(axis=1)
    return df_p

# simulates c iterations of fixed mutations counts arising by a binomial-Poisson process.
def mk_df_bp(chr_mu_dict,c):
    data_poisson2_p = chr_mut_bp(chr_mu_dict,c)
    dfas_p = pd.DataFrame.from_dict(data_poisson2_p)
    df2_p = dfas_p.applymap(binomialize) #Mimics segregation of each chromosome.
    df2_p["sum"]=df2_p.sum(axis=1)
    return df2_p

n = 10000 # desired number of iterations to perform.

#Gathers the number of nucleotides on each chr in human genome (a and b for 2n cell)
chr_list_a = []
chr_list_b = []
chr_lenlist = []
for chrfile in glob.glob(r"G:\My Drive\Volatility File Share\HumanGenome\*.fa"):
    chr_rm = chrfile.split(".")[4] #gets the chr number from the fasta file name.
    chr_list_a.append(chr_rm + "a")
    chr_list_b.append(chr_rm + "b")
    with open(chrfile, 'r') as f:
        firstline = f.readline()
        count = int(firstline.split(":")[5]) #gets chr length from header of fasta file
    chr_lenlist.append(count)
tot_bp = sum(chr_lenlist)

# Extrapolating yeast mutation rates to human ultramutator cells.
yeast_1n_mu = 34.5 #mutations/divisions/genome
print("With a yeast haploid mutation rate of", yeast_1n_mu)
hu_mu = (tot_bp/11e7)*yeast_1n_mu
print("the comparable human haploid mutation rate is", hu_mu)

#Generates initial Poisson and binomial-Poisson distributions for comparisons.
mu=hu_mu
counter = 0
Dlist=[]
D2list=[]
chr_mu_dict = mk_chr_mu_dict(mu,chr_lenlist,tot_bp,chr_list_a,chr_list_b)
df_p = mk_df_p(chr_mu_dict,n)
df_bp = mk_df_bp(chr_mu_dict,n)
print("Index of Dispersion for Poisson, n=",n,",",index_dis(df_p))
print("Index of Dispersion for binomial-Poisson, n=",n,",",index_dis(df_bp))

# Simulates the change in mutation burden in individual tumor cell lineages over 10 divisions.
n = 1000
counter = 0
plist=[] #
bplist=[]
dfRT=pd.DataFrame() #To keep track of cumulative mutation burden during a Poisson process.
dfRT["Div"]=list(range(1,31,1))
df2RT=pd.DataFrame() #To keep track of cumulative mutation burden during a binomial-Poisson process.
df2RT["Div"]=list(range(1,31,1))
while counter < n:
    df_p_line = mk_df_p(chr_mu_dict,30)
    dfRT[counter+1]=df_p_line["sum"].cumsum()
    mutbur_p=df_p_line['sum'].sum()
    plist.append(mutbur_p)
    df_bp_line = mk_df_bp(chr_mu_dict,30)
    df2RT[counter+1]=df_bp_line["sum"].cumsum()
    mutbur_bp=df_bp_line['sum'].sum()
    bplist.append(mutbur_bp)
    counter += 1
print("The average mutation burden for Poisson after 30 div is",mean(plist),"+/-",stdev(plist))
print("The average mutation burden for binomial-Poisson after 30 div is",mean(bplist),"+/-",stdev(bplist))
print("Index of Dispersion for Poisson after 30 div, n=",n,",",(stdev(plist))**2/mean(plist))
print("Index of Dispersion for binomial-Poisson after 30 div, n=",n,",",(stdev(bplist))**2/mean(bplist))

# Introduces a row of zeros into each dataframe to represent the starting point for each lineage.
dfRT.loc[len(dfRT)]=0
df2RT.loc[len(dfRT)]=0
dfRTm=pd.melt(dfRT,id_vars=['Div'],value_vars=list(dfRT.columns)[1:],var_name='line', value_name='Total Mut')
df2RTm=pd.melt(df2RT,id_vars=['Div'],value_vars=list(df2RT.columns)[1:],var_name='line', value_name='Total Mut')

sns.set(style="white", palette="bright", color_codes=True, rc={"lines.linewidth": 1.0})
plt.rcParams['patch.linewidth'] = 0

# Plots Poisson vs binomial-Poisson distributions after 1 division.
fig1, ax = plt.subplots(figsize=(2.5, 2.5))

sns.distplot(df_p["sum"],
                  bins=list(range(900,3000,10)),
                  kde=False,
                  norm_hist=True,
                  color='gray')
sns.distplot(df_bp["sum"],
                  bins=list(range(900,3000,10)),
                  kde=False,
                  norm_hist=True,
                  color='orange')
ax.tick_params(bottom=True,left=True,labelsize=10)
ax.set(xlabel='Mutations/Division', ylabel='Density')

fig1.savefig(f"{File_save_location}/ASmut_hu.pdf",transparent = True,bbox_inches='tight')

# Plots Poisson vs binomial-Poisson distributions after 1 division.
fig2, ax = plt.subplots(figsize=(2.5, 2.5))

sns.distplot(plist,
                  bins=list(range(51000,63000,200)),
                  kde=False,
                  norm_hist=True,
                  color='gray')
sns.distplot(bplist,
                  bins=list(range(51000,63000,200)),
                  kde=False,
                  norm_hist=True,
                  color='orange')
ax.tick_params(bottom=True,left=True,labelsize=10)
ax.set(xlabel='Mutation Burden', ylabel='Density')

fig2.savefig(f"{File_save_location}/MB_hu30.pdf",transparent = True,bbox_inches='tight')

# Plots Poisson vs binomial-Poisson mutation accumulation over 10 divisions.
fig3, ax = plt.subplots(figsize=(2,2.5))
ax.tick_params(bottom=True,left=True,labelsize=10)
ax.set(xlim=(0, 10),ylim=(0,25000))
sns.lineplot(x="Div", y="Total Mut", hue='line', palette=sns.color_palette("Set1", n), legend=False,lw=0.5,data=df2RTm)
sns.lineplot(x="Div", y="Total Mut",legend=False, lw=1,color="k",data=dfRTm)

fig3.savefig(f"{File_save_location}/MB_huMutAc.pdf",transparent = True,bbox_inches='tight')
