import os
import glob
import sys

import collections
from collections import defaultdict

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, Gapped
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import argparse

parse = argparse.ArgumentParser()

parse.add_argument("--path",type=str, help="path to vespasian output",required=True)

parse.add_argument("--AA_fasta",type=str,help="path to and name of the original aligned AA fasta file that analysis was run on",required=True)
#parse.add_argument("--species",type=str,help="name of species in alignment file under analysis")
parse.add_argument("--species_list",type=str,help="path to and name of branches file that vespasian analysis was run on",required=True)

parse.add_argument("--b_range", type=int, default=20,help="number of sites either side of selected site to include in blanks analysis default:20") 
parse.add_argument("--c_range", type=int, default=5,help="number of sites either side of selected site to include in conservation analysis default:5") 
parse.add_argument("--gaps_cutoff", type=float, default=0.75,help="filter for number of gaps allowed around selected site default:0.75") 
parse.add_argument("--dist", type=float, default=0.1,help="filter for blank count either side of selected site for default:0.1") 

args = parse.parse_args()


path = args.path
aligned_fasta = args.AA_fasta
species  = args.species
blanks_site_range = args.b_range
conserv_site_range  = args.c_range
gaps_cutoff  = args.gaps_cutoff
similarity_cutoff  = args.dist


def passed_final_output(dict1, dict2, dict3):
    '''Take return of passed sites function and export a CSV of outputs'''
    combine_dict = {**dict1 , **dict2, **dict3}
    df = pd.DataFrame.from_dict(combine_dict, orient = 'index')
    return df.to_csv('passed_positively_selected_sites.tsv', sep='\t', index_label = 'Gene_family', header = False)


def failed_final_output(dict1, dict2, dict3):
    '''Take return of failed sites function and export a CSV of outputs'''
    combine_dict = {**dict1 , **dict2, **dict3}
    df = pd.DataFrame.from_dict(combine_dict, orient = 'index')
    return df.to_csv('failed_positively_selected_sites.tsv', sep='\t', index_label = 'Gene_family', header = False)
    

def parse_species(branches):
    '''Function to get list of multiple branch labelled species (if 
    that option was used in vespasian analysis)'''
    branches_file = args.species_list
    if isinstance(branches_file, str):
        species_label = branches_file
    else:
        with open(branches_file) as f:
            for line in f:
                lines = line.split(':')
                species_label = lines[1].strip()
                species_label = species_label.strip('[').strip(']').replace(', ', ',').split(',')

    return(species_label)


def report_parse(path):
    posSel_site_dict = dict()
    original_align_fasta = aligned_fasta#requires path and file name for original AA alignemnt file

    seq_len = []
    for record in SeqIO.parse(original_align_fasta, 'fasta'):
        ID = record.id
        seq = str(record.seq)
        posSel_site_dict[ID] = seq
        seq_len.append(len(seq))
    seq_len = seq_len[0]

    AA_position_sites = []
    vespasian_summary = path + 'summary.tsv'
    with open(vespasian_summary) as f:
        for line in f:
            lines = line.split('\t')
            fam = lines[0]
            tree = lines[1]
            model = lines[2]
            pos_sel_sites = lines[6]
            neb_sites = lines[7]
            beb_sites = lines[8]
            if pos_sel_sites == 'yes':
                beb_sites_split = beb_sites.split(' ')
                for site in beb_sites_split:
                    locus = site.split(':')[0]
                    position = locus[:-1]
                    AA = locus[-1]
                    AA_position_sites.append(position)

    posSel_sites = []
    for pos in range(seq_len):
        pos = pos + 1
        if str(pos) in AA_position_sites:
            posSel_sites.append('N')
        else:
            posSel_sites.append('-')
    posSel_sites = ''.join(posSel_sites)
    posSel_site_dict['PS_Sites|'] = posSel_sites

    return posSel_site_dict


def retrieve_selected_sites(path):
    '''Use alignment outputted from codeml_reader to give 
    co-ordinates of positively selected sites'''

    selected_sites = defaultdict(list)
    
    posSel_site_dict = report_parse(path)

    selected_sequence = []
    positively_selected_sites = []
    for ID, seq in posSel_site_dict.items():
        if 'PS_Sites|' in ID:
            for character in seq:
                selected_sequence.append(character)

    for position, char in enumerate(selected_sequence):
        if char != '-':
            positively_selected_sites.append(position)

    family = path.rpartition('/')[0].rpartition('/')[2]

    selected_sites[family] = positively_selected_sites
    sel_site_count = [len(site) for key, site in selected_sites.items()]
    print(sel_site_count[0], 'positively selected sites')

    return selected_sites


def unique_selected_sites(path, sites, species):
    '''From all sites under positive selection find only those that have amino acid differences
    unique to the species specified'''
    
    dictionary = defaultdict(list)
    fail_dictionary = defaultdict(list)

    for family, positions in sites.items():
        pos_sites_count = len(positions)

        posSel_site_dict = report_parse(path)

        seq_dict = {} 
        for ID, seq in posSel_site_dict.items():
            if ID[:8] != 'PS_Sites':
                ID = ID.split('|')[0]
                seq_dict[ID] = seq

        #Write list of positively selected sites in each species, see if selected species is unique and all other conserved
        passed_sites = []
        for index in positions: #iterate through all positively selected sites
            index = int(index)
            species_character = []
            other_species = []
            
            for key, value in seq_dict.items():#if species' amino acid is not a gap add all species to dictionary 
                if (key in species) and value[index] == '-' :
                    pass
                if (key in species) and value[index] != '-' :
                    species_character.append(value[index])
                if (key not in species):
                    other_species.append(value[index])
            if '-' in species_character:
                species_character.remove('-')
            if '-' in other_species:
                other_species.remove('-')
            if (len(set(species_character)) == 1 ) and (len(set(other_species)) == 1) and (x not in species_character for x in other_species) :
                passed_sites.append(str(index + 1))#Increase index number by one to account for python array

        dictionary[family] = passed_sites
        passed_sites_count = len(passed_sites)

        failed_unique_sites = []#Save +ve selected sites that fail the unique sites filter
        for site in positions:
            site = site + 1
            site = str(site)
            if site not in passed_sites:
                failed_unique_sites.append(site)
        fail_dictionary[family] = failed_unique_sites

    print(pos_sites_count - passed_sites_count, 'positively selected sites not unique')

    return dictionary, fail_dictionary

def pull_sites(selected_sites, site_range):
    '''From an output file containing gene_families and selected site, pull out selected sites
    and defined flanking regions'''
    
    gene_dict = defaultdict(dict)
    for gene_id, sites in selected_sites.items():
        sites_dict = {}
        for site in sites:
            surrounding_sites = []
            range_max = int(site) + int(site_range) + 1 # zero indexed range
            range_min = int(site) - int(site_range) - 1 # zero indexed range
            [{surrounding_sites.append(i)} for i in range(int(site) + 1 , range_max, 1)] #don't include original site
            [{surrounding_sites.append(i)} for i in range(int(site) - 1 , range_min, -1)]
            surrounding_sites.sort()
            sites_dict[site] = surrounding_sites
        gene_dict[gene_id] = sites_dict

    return gene_dict

def count_blanks(path, gene_dict, gaps_cutoff, species):
    '''Iterate through sites with unique amino acid substitutions and check defined flanking regions filtering 
    any sites that contain more gaps than threshold'''

    blank_passed_sites = defaultdict(list)
    blank_failed_sites = defaultdict(list)

    for family in gene_dict.keys():
        original_dict = report_parse(path)

        seq_dict = {}
        for ID, seq in original_dict.items():
            ID = ID.split('|')[0]
            if ID in species:
                seq_dict[ID] = seq

        for gene, sel_site in gene_dict.items():
            sel_site_count = len(sel_site)
            if gene == family:
                sites = sel_site

        if isinstance(species, str):#If analysis a single species or multiple species
            passed = []
            failed = []
            for sel_site, flank_sites in sites.items():
                py_index = [x-1 for x in flank_sites] #account for zero index
                region = []
                for index in py_index:
                    for k, v in seq_dict.items():
                        if index < len(v):
                            region.append(v[index])

                if region.count('-') < (gaps_cutoff*len(py_index)):
                    passed.append(sel_site)
                else:
                    failed.append(sel_site)
                    
            passed_sites_count = len(passed)
            blank_passed_sites[family] = passed
            blank_failed_sites[family] = failed

        else:
            passed = []
            failed = []
            for sel_site, flank_sites in sites.items():
                py_index = [x-1 for x in flank_sites] #account for zero index
                species_region = {}
                for k, v in seq_dict.items():
                    region = []
                    for index in py_index:
                        if index < len(v):
                            region.append(v[index])
                        species_region[k] = region 

                sp_passed = {}
                for sp, region in species_region.items():
                    passed_sites = []
                    if region.count('-') < (gaps_cutoff*len(py_index)):
                        passed_sites.append(sel_site)
                    else:
                        continue
                    sp_passed[sp] = passed_sites

                if len(sp_passed) >= 1:#If at least one species from your list passes gap filter, add the +ve site to the dictionary
                    sites_passed = sp_passed.values()
                    sites_list = []
                    for sites in sites_passed:
                        sites_list.append(sites[0])
                    passed_site = sites_list[0]
                    passed.append(passed_site)
                    passed_sites_count = len(passed)
                    blank_passed_sites[family] = passed
                else:
                    failed.append(sel_site)
                    blank_failed_sites[family] = failed


    print(sel_site_count - passed_sites_count, 'sites dropped due to too many gaps')
    return blank_passed_sites, blank_failed_sites


def conservation_analysis(path, gene_dict, similarity_cutoff, species):

    conserv_passed_sites = defaultdict(list)
    conserv_rejected_sites = defaultdict(list)  
    
    for family in gene_dict.keys():
        passed_sites = []
        rejected_sites = []

        original_dict = report_parse(path)
       
        #Grab all the sites in the flanks of the selected sites 
        seq_dict = {}
        sp_list = []
        for ID, seq in original_dict.items():
            if ID[:8] != 'PS_Sites':
                ID = ID.split('|')[0]
                sp_list.append(ID)
                seq_dict[ID] = seq

        for gene, sel_site in gene_dict.items():
            if gene == family:
                sites = sel_site

        #Take the flanked sites and combine them into a matrix
        for sel_site, flank_sites in sites.items():
            py_index = [x-1 for x in flank_sites] #account for zero index
            region = {}
            for k, v in seq_dict.items():
                region_seqs = []
                for index in py_index:
                    if index < len(v):
                        region_seqs.append(v[index])
                region[k]=(region_seqs)
            df = pd.DataFrame.from_dict(region, orient="index")

            all_chrs = df.iloc[1:, :].sum(0)
            consensus = []
            for i in range(0, len(all_chrs)):
                unique_chrs = "/".join(set(all_chrs[i])) 
                most_freq_chr = collections.Counter(all_chrs[i]).most_common(1)[0]
                consensus.append(most_freq_chr[0]) #output most common amino acid

            df.loc[len(df)] = consensus

            species_seq  = (df.iloc[:, :].sum(1)[species])

            consensus_seq = (df.iloc[[-1]].sum(1))
            for con_ID, con_seq in consensus_seq.items():
                con_seq = con_seq


            with open('temp.txt', 'w') as interim_file:   
                if isinstance(species, str):#If analysis a single species or multiple species
                    interim_file.write('>' + species + '\n' + species_seq + '\n' + '>consenus' + '\n' + con_seq + '\n')
                else:
                    for sp_ID, sp_seq in species_seq.items():
                        interim_file.write('>' + sp_ID + '\n' + sp_seq + '\n')
                    interim_file.write('>consenus' + '\n' + con_seq + '\n')

            aln = AlignIO.read('temp.txt', 'fasta')
            os.remove('temp.txt')
        
            calculator = DistanceCalculator('blosum62')
            dm = calculator.get_distance(aln)


            if isinstance(species, str):#If analysis a single species or multiple species
                if dm[0][1] < similarity_cutoff:
                    passed_sites.append(sel_site)
                else:
                    rejected_sites.append(sel_site)
            else:
                species_count = len(species)
                dist_list = []
                for pos in range(species_count):
                    dist_list.append(dm[pos][-1])


                if any(dist < similarity_cutoff for dist in dist_list):#Conservation is passed if any species seq passes filter
                    passed_sites.append(sel_site)
                else:
                    rejected_sites.append(sel_site)

##Currently, sites pass if at least one species satisfies filter, use below code if you want an average across species instead
#                avg_dist_species = sum(dist_list)/len(dist_list)#Conservation is passed if the average distance score passes filter
#                if avg_dist_species < similarity_cutoff:
#                    passed_sites.append(sel_site)
#                else:
#                    rejected_sites.append(sel_site)

        conserv_passed_sites[family] = passed_sites
        conserv_rejected_sites[family] = rejected_sites

    reject_count = [len(site) for key, site in conserv_rejected_sites.items()]
    print(reject_count[0], 'sites failed conservation filter')
    
    return conserv_passed_sites, conserv_rejected_sites


species_list  = parse_species(args.species_list)

sites = retrieve_selected_sites(path)

#If just unique sites function is run - need to create some flag specific to this function
#unique_selected_sites(path, sites, species_list)


#Write results to file for analysis using all 3 functions
if (blanks_site_range is not None) and (conserv_site_range is not None) and (gaps_cutoff is not None) and (similarity_cutoff is not None):
    gene_dict_blank = pull_sites(sites, blanks_site_range)
    gene_dict_conserv = pull_sites(sites, conserv_site_range)

    unique_sites_passed = unique_selected_sites(path, sites, species_list)
    blank_count_pass = count_blanks(path, gene_dict_blank, gaps_cutoff, species_list)
    conserv_score_pass = conservation_analysis(path, gene_dict_conserv, similarity_cutoff, species_list)
    passed_final_output(unique_sites_passed[0], blank_count_pass[0], conserv_score_pass[0])
    failed_final_output(unique_sites_passed[1], blank_count_pass[1], conserv_score_pass[1])



