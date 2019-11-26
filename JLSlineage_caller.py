#!/usr/bin/env python


from collections import defaultdict
from argparse import ArgumentParser
import pandas as pd


def main():

    parser = ArgumentParser()
    parser.add_argument("--lineage_info", action='store', dest="lineage_info", help='csv file consisting of the strain file name followed by the index in the lineage', required=True)
    parser.add_argument("--min_depth", action='store', dest="min_depth", type=int, help="Minimum depth to evaluate site", default=18)
    o = parser.parse_args()

    strain_dict = defaultdict(lambda: 0)
    scored_sites = defaultdict(lambda: 0)
    var_dict = defaultdict(lambda: [])

    lineage_info_file = open(o.lineage_info, 'r')

    for line in lineage_info_file:
        # creates a list of a line from the line of the CSV
        splitline = line.strip('\r\n').split(',')
        # creates a dictionary entry where key is the rep_event and the value is the above list
        strain_dict.update({int(splitline[0]): splitline[1:]})

    lineage_info_file.close()

    site_count_file = open("site_counts.txt", 'w')

    sub_lineage_list = []

    for rep_event in strain_dict.keys():
        print rep_event
        # For the first two strains it goes through the pileup files and counts and scores the available sites that have sufficient coverage
        for sub_lineage in strain_dict[rep_event][:2]:
            sub_lineage_list.append(sub_lineage)
            nuc_ctr = 0
            pileup_file = open(sub_lineage + '_pe.mpileup', 'r')

            for line in pileup_file:

                splitline = line.strip('\n').split('\t')

                if int(splitline[3]) > o.min_depth:
                    nuc_ctr += 1
                    scored_sites[splitline[0] + ":" + str(splitline[1])] += 1

            pileup_file.close()

            site_count_file.write(sub_lineage + '\t' + str(nuc_ctr) + '\n')

    scored_site_file = open("scored_sites.txt", 'w')
    scored_site_ctr = 0
    #creates a shared set for counting the shared sites.
    shared_set=set()
    for site in scored_sites.keys():
        #From the scored_sites dictionary that kept track of sites with at least 18-fold coverage in one strain, we now assess how many times we counted it and compare it to the
        #the number of strains in our sublineage list.  For example, if we have 12 strains in our isolate_list and we observed it 12 times than we can add the site to the shared site set.
        if scored_sites[site] == len(sub_lineage_list):
 	    shared_set.add(site)
            scored_site_ctr += 1
            scored_site_file.write(site + '\n')

    site_count_file.write("# Scored Sites:\t" + str(scored_site_ctr))

    scored_site_file.close()
    site_count_file.close()

    total_var_set = set()
    total_var_dict = defaultdict(lambda:set())
    for rep_event in strain_dict.keys():

        temp_list = []
        for sub_lineage in strain_dict[rep_event]:

            var_file = open(sub_lineage + ".noSNP.var", 'r')
            var_set = set()

            for line in var_file:
                splitline = line.strip('\n').split('\t')
                # all vars have the required coverage.  This asks if they are in the shared_set.  If they are then we score it.
                if splitline[0] + ":" + splitline[1] in shared_set and splitline[0] != "Chrom":
                    var_set.add(splitline[0] + ":" + splitline[1])
                    total_var_set.add(splitline[0] + ":" + splitline[1])

                    if splitline[1] not in total_var_dict[splitline[0]]:
                        total_var_dict[splitline[0]].add(int(splitline[1]))

            temp_list.append(var_set)

        strain_dict[rep_event] = temp_list

    pos_file = open("position_file.txt", 'w')
    for chrom in sorted(total_var_dict.keys()):

        pos_file.write(chrom + '\n')
        x = sorted(list(total_var_dict[chrom]))

        for pos in x:
            pos_file.write(str(pos) + '\n')

    pos_file.close()


    rep_event_keys = sorted(strain_dict.keys())

    df_master = pd.DataFrame()
    total_var_list = sorted(list(total_var_set))
    for x in rep_event_keys:

        dgd2_name = 'd' + str(x) + '_gd' + str(x) + '-2'
        gd1ggd1_name = 'gd' + str(x) + '-1_ggd' + str(x) + '-1'
        dgd1_name = 'd' + str(x) + '_gd' + str(x) + '-1'
        ddnplus1_name = 'd' + str(x) + '_d' + str(x + 1)

        df = pd.DataFrame({dgd2_name : pd.Series(sorted(list(strain_dict[x][0].intersection(strain_dict[x][2]))),
                                   index=sorted(list(strain_dict[x][0].intersection(strain_dict[x][2])))),
                                   gd1ggd1_name : pd.Series(sorted(list(strain_dict[x][1].intersection(strain_dict[x][3]))),
                                  index=sorted(list(strain_dict[x][1].intersection(strain_dict[x][3])))),
                                   dgd1_name: pd.Series(sorted(list(strain_dict[x][0].intersection(strain_dict[x][1]))),
                                    index=sorted(list(strain_dict[x][0].intersection(strain_dict[x][1])))),
                                    ddnplus1_name : pd.Series(sorted(list(strain_dict[x][0].intersection(strain_dict[x][4]))),
	                                index=sorted(list(strain_dict[x][0].intersection(strain_dict[x][4]))))}).reindex(total_var_list)

        #df.to_csv(str(x) + '_test.csv')
        df_master = pd.concat([df_master,df], axis=1)
        #print df_master
    df_master.to_csv("test.csv")

if __name__ == "__main__":
    main()
