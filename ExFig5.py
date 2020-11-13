# For Extended Data Fig.5, Dowsett et al.
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
"""The following Python 3 code simulates the extent to which the mutation counts in
pairs of segregant groups (Da vs Db and Ma vs Mb) are correlated."""

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

# simulates c iterations of fixed mutations counts arising and segregating to two cells.
def mk_df_bp(chr_mu_dict,c):
    data_poisson2_p = chr_mut_bp(chr_mu_dict,c)
    dfas_p = pd.DataFrame.from_dict(data_poisson2_p)
    df_binomialize(dfas_p) #produces two dfs of mutation counts from reciprocol segregants
    df2_cell1["sum"]=df2_cell1.sum(axis=1)
    add_chr_totals(df2_cell1)
    df2_cell2["sum"]=df2_cell2.sum(axis=1)
    add_chr_totals(df2_cell2)

def df_binomialize(df_input):
    df_input_bintem = (df_input * 0) + 1 #makes a dataframe of 1's w/dimensions of df_input
    df_input_binom = df_input_bintem.applymap(binomialize) #Randomly changes 1's to 0 for cell1
#     print(df_input_binom)
    df_input_binom_inv=abs(df_input_binom - 1) #Inverse of dfas_p_binom for cell2
#     print(df_input_binom_inv)
    df2_cell1 = df_input.mul(df_input_binom) #Multiplies mutations in cell1 by 1 or 0 to mimic segregation.
    df2_cell2 = df_input.mul(df_input_binom_inv) #Mimics mutations in cell2 the inverse to mimic segregation.
    df2_cell1["sum"] = df2_cell1.sum(axis=1)
    add_chr_totals(df2_cell1)
    df2_cell2["sum"] = df2_cell2.sum(axis=1)
    add_chr_totals(df2_cell2)
    return df2_cell1,df2_cell2

def add_chr_totals(df_cell): #adds muts on a and b sets of chromosomes to get observed mut/chr.
    df_cell["I"] = df_cell["Ia"] + df_cell["Ib"]
    df_cell["II"] = df_cell["IIa"] + df_cell["IIb"]
    df_cell["III"] = df_cell["IIIa"] + df_cell["IIIb"]
    df_cell["IV"] = df_cell["IVa"] + df_cell["IVb"]
    df_cell["V"] = df_cell["Va"] + df_cell["Vb"]
    df_cell["VI"] = df_cell["VIa"] + df_cell["VIb"]
    df_cell["VII"] = df_cell["VIIa"] + df_cell["VIIb"]
    df_cell["VIII"] = df_cell["VIIIa"] + df_cell["VIIIb"]
    df_cell["IX"] = df_cell["IXa"] + df_cell["IXb"]
    df_cell["X"] = df_cell["Xa"] + df_cell["Xb"]
    df_cell["XI"] = df_cell["XIa"] + df_cell["XIb"]
    df_cell["XII"] = df_cell["XIIa"] + df_cell["XIIb"]
    df_cell["XIII"] = df_cell["XIIIa"] + df_cell["XIIIb"]
    df_cell["XIV"] = df_cell["XIVa"] + df_cell["XIVb"]
    df_cell["XV"] = df_cell["XVa"] + df_cell["XVb"]
    df_cell["XVI"] = df_cell["XVIa"] + df_cell["XVIb"]


# Input parameters from glm.nb in R for the combined distribution of Dm and Mm counts.
munb = 4.927 # rate parameter for negative binomial (munb)
theta = 60.42 # shape parameter

# Conversions
muc = math.exp(munb) #mu corrected (muc) is munb exponentiated.
gv = (muc)**2/theta #variance of the Gamma distribution of lambda
gscale = gv/muc #scale paramater of the Gamma distribution of lambda.
hap_mu = muc/4 #rate of fixed mutations per haploid genome.
print("The haploid rate of fixed mutations is", hap_mu)

n = 100 # number of iterations to perform.

#Counts the number of unmasked nucleotides on each chr in the masked haploid yeast genome
chr_list_a = [] #(list_a and list_b for two sets of chrs in 2n cell)
chr_list_b = []
chr_list = []
chr_lenlist = []
for chrfile in glob.glob(r"G:\My Drive\Volatility File Share\YeastGenome\*_rm.*"):
    chr_rm = chrfile.split(".")[4] #gets the chr number from the fasta file name.
    chr_list_a.append(chr_rm + "a")
    chr_list_b.append(chr_rm + "b")
    chr_list.append(chr_rm)
    with open(chrfile, 'r') as f:
        firstline = f.readline()
        count = int(firstline.split(":")[5]) #gets chr length from header of fasta file
        for line in f:  #counts "N"s in the fasta file and subtracts them from chr length.
            for i in line:
                if i == "N":
                    count = count - 1
    chr_lenlist.append(count)
tot_bp = sum(chr_lenlist)

#Models the occurrance and segregation of fixed mutations in pairs of segregating cells.

"""The following returns a gamma distribution of rates of mismatch formation
and then divides each value by 4, since we're looking at the number of fixed mutations
per haploid genome. This distribution will be our input for n divisions,
each using a different lambda. It then takes the dataframe of mismatch counts and generates
the reciprocol mutation burdens expected in the resulting pairs of segregant cells."""

data_gamma = (gamma.rvs(a=theta, scale=gscale, size=n))/4
df_nb = pd.DataFrame(columns=chr_list_a+chr_list_b)
dfas_nbD = pd.DataFrame(columns=chr_list_a+chr_list_b) #D stands for Daughter cell
dfas_nbM = pd.DataFrame(columns=chr_list_a+chr_list_b) #M stands for Mother cell
for g in data_gamma: #with the same lambda we're going to call mutations twice and compare
    chr_mu_dictD = mk_chr_mu_dict(g,chr_lenlist,tot_bp,chr_list_a,chr_list_b)
    chr_mu_dictM = mk_chr_mu_dict(g,chr_lenlist,tot_bp,chr_list_a,chr_list_b)
    data_poisson2D = chr_mut_bp(chr_mu_dictD, 1)
    data_poisson2M = chr_mut_bp(chr_mu_dictM, 1)
    dfa2D = pd.DataFrame.from_dict(data_poisson2D)
    dfa2M = pd.DataFrame.from_dict(data_poisson2M)
    dfas_nbD = dfas_nbD.append(dfa2D, ignore_index=True)
    dfas_nbM = dfas_nbM.append(dfa2M, ignore_index=True)
df2_bnb_cell1D,df2_bnb_cell2D = df_binomialize(dfas_nbD)
df2_bnb_cell1M,df2_bnb_cell2M = df_binomialize(dfas_nbM)
df2totalbnbDsums=pd.DataFrame()
df2totalbnbMsums=pd.DataFrame()
df2totalbnbDsums["Segregant 1"]=df2_bnb_cell1D["sum"]
df2totalbnbDsums["Segregant 2"]=df2_bnb_cell2D["sum"]
df2totalbnbMsums["Segregant 1"]=df2_bnb_cell1M["sum"]
df2totalbnbMsums["Segregant 2"]=df2_bnb_cell2M["sum"]
df2totalbnbDsums

"""Creates data frame of actual data (Da v Db and Ma v Mb) for comparison."""
df_data = pd.read_csv("G:\My Drive\Volatility File Share\Data\MutInSegGroups.csv")

"""Plots the linear regressions of the bnb model and the actual data."""
fig1, ax = plt.subplots(figsize=(2,2))
sns.regplot(x="Segregant 1", y="Segregant 2", data=df2totalbnbDsums, marker='o', color="silver", scatter_kws={'s':7},line_kws={'color':"peru"})
sns.regplot(x="Db", y="Da", data=df_data, color="green",marker='o', scatter_kws={'s':11})
sns.regplot(x="Ma", y="Mb", data=df_data, color="blue",marker='o', scatter_kws={'s':11})
ax.set(ylim=(0, 150))
ax.set(xlim=(0, 150))
ax.set(xlabel='Segregant 1', ylabel='Segregant 2')
ax.tick_params(bottom=True,left=True, labelsize=8)
fig1.savefig(f"{File_save_location}/lineRegSegsnb100v2.pdf",transparent = True,bbox_inches='tight')
