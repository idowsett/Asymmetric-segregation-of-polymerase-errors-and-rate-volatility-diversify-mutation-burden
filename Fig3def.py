# For Fig.3d,e,f, Dowsett et al. (pol3-01/pol3-01 msh6/msh6 diploid cells)
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
"""The following Python 3 code simulates the extent to which mutator volatility
and asymmetric segregation of mutations lead to overdispersion  It first simulates
the contribution of asymmetric segregation alone, assuming a single Poisson process
governs the mutator phenotype in all cell divisions. The script then uses the parameters for
a negative binomial model of mismatches segregated to daughter (Dm) and mother (Mm)
to define a gamma distribution of lambda values governing the underlying mutation
rates of individual divisions. With these values, the script then
simulates the distribution of error counts in individual divisions."""

# Output Files(s) folder:
File_save_location = "C:/Users/Sett/Desktop/Temp Photos"

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

# Input parameters from glm.nb in R for the combined distribution of Dm and Mm counts.
munb = 4.927 # rate parameter for negative binomial (munb)
theta = 60.42 # shape parameter

# Conversions
muc = math.exp(munb) #mu corrected (muc) is munb exponentiated.
gv = (muc)**2/theta #variance of the Gamma distribution of lambda
gscale = gv/muc #scale paramater of the Gamma distribution of lambda.
hap_mu = muc/4 #rate of fixed mutations per haploid genome.
print("The haploid rate of fixed mutations is", hap_mu)

n = 10000 # number of iterations to perform.

#Counts the number of unmasked nucleotides on each chr in the masked haploid yeast genome
chr_list_a = [] #(list_a and list_b for two sets of chrs in 2n cell)
chr_list_b = []
chr_lenlist = []
for chrfile in glob.glob(r"G:\My Drive\Volatility File Share\YeastGenome\*_rm.*"):
    chr_rm = chrfile.split(".")[4] #gets the chr number from the fasta file name.
    chr_list_a.append(chr_rm + "a")
    chr_list_b.append(chr_rm + "b")
    with open(chrfile, 'r') as f:
        firstline = f.readline()
        count = int(firstline.split(":")[5]) #gets chr length from header of fasta file
        for line in f:  #counts "N"s in the fasta file and subtracts them from chr length.
            for i in line:
                if i == "N":
                    count = count - 1
    chr_lenlist.append(count)
tot_bp = sum(chr_lenlist)

#Generates Poisson and binomial-Poisson distributions.
mu=hap_mu
Dlist=[]
D2list=[]
chr_mu_dict = mk_chr_mu_dict(mu,chr_lenlist,tot_bp,chr_list_a,chr_list_b)
df_p = mk_df_p(chr_mu_dict,n)
df2_p = mk_df_bp(chr_mu_dict,n)
print("Index of Dispersion for Poisson, n=",n,",",index_dis(df_p))
print("Index of Dispersion for binomial-Poisson, n=",n,",",index_dis(df2_p))

#simulates variation in index of dispersion given a small sample size
counter = 0
sample = 200
while counter < n:
    df_ps = mk_df_p(chr_mu_dict,sample)
    df2_ps = mk_df_bp(chr_mu_dict,sample)
    Dlist.append(index_dis(df_ps))
    D2list.append(index_dis(df2_ps))
    counter += 1
print("The mean index of dispersion for Poisson, n=",sample,",", mean(Dlist),"+/-",stdev(Dlist))
print("The mean index of dispersion for binomial-Poisson, n=", sample,",",mean(D2list),"+/-",stdev(D2list))

#Generates negative binomial and binomial-negative binomial distributions.
"""The following returns a gamma distribution of rates of mismatch formation
and then divides each value by 4, since we're looking at the number of fixed mutations
per haploid genome. This distribution will be our input for n divisions,
each using a different lambda."""
data_gamma = (gamma.rvs(a=theta, scale=gscale, size=n))/4
df_nb = pd.DataFrame(columns=chr_list_a+chr_list_b)
dfas_nb = pd.DataFrame(columns=chr_list_a+chr_list_b)
for g in data_gamma:
    chr_mu_dict = mk_chr_mu_dict(g,chr_lenlist,tot_bp,chr_list_a,chr_list_b)
    data_poisson = chr_mut_p(chr_mu_dict, 1)
    data_poisson2 = chr_mut_bp(chr_mu_dict, 1)
    dfa = pd.DataFrame.from_dict(data_poisson)
    df_nb = df_nb.append(dfa, ignore_index=True)
    dfa2 = pd.DataFrame.from_dict(data_poisson2)
    dfas_nb = dfas_nb.append(dfa2, ignore_index=True)
df2_nb = dfas_nb.applymap(binomialize)
df_nb["sum"] = df_nb.sum(axis=1)
df2_nb["sum"] = df2_nb.sum(axis=1)
print("Index of Dispersion for negative binomial, n=",n,",",index_dis(df_nb))
print("Index of Dispersion for binomial-negative binomial, n=",n,",",index_dis(df2_nb))

sns.set(style="white", palette="bright", color_codes=True, rc={"lines.linewidth": 1.0})
plt.rcParams['patch.linewidth'] = 0

# Plots Poisson vs binomial-Poisson distributions.
fig1, ax = plt.subplots(figsize=(2.5, 2.5))
sns.distplot(df_p["sum"],
                  bins=list(range(20,120,1)),
                  kde=False,
                  norm_hist=True,
                  color='gray')
sns.distplot(df2_p["sum"],
                  bins=list(range(20,120,1)),
                  kde=False,
                  norm_hist=True,
                  color='orange')
ax.tick_params(bottom=True,left=True,labelsize=10)
ax.set(xlabel='Mutations/Division', ylabel='Density')

fig1.savefig(f"{File_save_location}/ASmut_ye2n.pdf",transparent = True,bbox_inches='tight')

#Plots the distributions of the indices of dispersion from Poisson and binomial-Poisson simulations.
fig2, ax = plt.subplots(figsize=(2.5, 2.5))

sns.distplot(Dlist,
                bins=[i for i in np.arange(0.25,5,0.02)],
                  kde=False, norm_hist=True,
                  color='gray')
ax = sns.distplot(D2list,
                  bins=[i for i in np.arange(0.25,5,0.02)],
                  kde=False, norm_hist=True,
                  color='orange')
ax.tick_params(bottom=True,left=True,labelsize=10)
ax.set(xlabel='Index of Dispersion', ylabel='Density')

fig2.savefig(f"{File_save_location}/IDmut_ye2n.pdf",transparent = True,bbox_inches='tight')

# Plots all distributions on one plot
fig4, ax = plt.subplots(figsize=(2.5, 2.5))
sns.distplot(df_p['sum'],
                  bins=list(range(20,120,2)),
             hist = False,
#                   kde=False,
                  norm_hist=True,
                  color='gray')
sns.distplot(df2_p['sum'],
                  bins=list(range(20,120,2)),
             hist = False,
#                   kde=False,
                  norm_hist=True,
                  color='slateblue')

sns.distplot(df_nb['sum'],
                  bins=list(range(20,120,2)),
             hist = False,
#                   kde=False,
                  norm_hist=True,
                  color='green')
sns.distplot(df2_nb['sum'],
                  bins=list(range(20,120,2)),
             hist = False,
#                   kde=False,
                  norm_hist=True,
                  color='salmon')
ax.tick_params(bottom=True,left=True,labelsize=10)
ax.set(xlim=(0, 150))
ax.set(xlabel='Mutations/Division', ylabel='Density')

fig4.savefig(f"{File_save_location}/CombinedModelalt_ye2n.pdf",transparent = True,bbox_inches='tight')
