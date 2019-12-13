Data Collection Pipeline:
Data collection performed using a Mac Mojave 10.14.5 environment.

A custom sequence analysis bash script called eex_yeast_pileline utilizes established tools including BWA (0.7.17), 
Samblaster (v.0.1.24), Picard-tools (2.2.2) , Samtools (1.8), GATK (4.0.6.0), Varscan (2.3.9), and a Python 2 script 
(Variant_deSNPer2013a.py) written to remove single nucleotide polymorphisms in our strain background from the variant calls.  
A Python 2 script (JLSlineage_caller.py) calls shared variants in our lineages.




Data Analysis:
Figures and extended figures generated using code are provided as standalone Python 3.6.5 or R v3.5.3 scripts with included statistical tests.
All figures and extended figure scripts run using a Windows 10 environment, these are named for the figure number.
Some editing for clarity of text/colors performed for generation of final figures provided in full manuscript.

Required 3rd party Python 3.6.5 packages (and all dependencies), available through PIP:
latest versions of Scipy, pandas, plotly, matplotlib, seaborn.
May be run through a standard Windows 10 command line (e.g. Powershell). 

An R script performs mixture modeling and AIC scoring for Fig1 (Fig1b.R). 
Required R v3.5.3 libraries (and all dependencies) for file Fig1b.R:
countreg, data.table, dplyr, flexmix, ggplot2
