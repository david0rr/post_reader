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



'''Order: 
        [1] get dictionary of all (working directory, species) >>> compulsory function *
        [2] get dictionary of unique sites (dictionary, species) >>> optional function *
        [3] pull blank sites (dictionary, blank_range) >>> optional_function
        [4] analyse blank sites (previous_dict, gap_cutoff, species) >>> compulsory function if function [3] performed
        [5] pull conserved sites (dictionary, conserv_site_range) >>> optional_function
        [6] analyse conserv_sites (previous_dict, conserv_site_range, species) >>> compulsory function if function [5] performed
        
        [7] Find unique outputs in final dict 
        [8] Convert to CSV 
'''

def function_set(path, species, options, unique_sub = False, blanks = False, conservation = False):
    '''Take input from wrapper and determine which functions to run'''
    
    blanks_site_range = int(options['b_range'])
    conserv_site_range = int(options['c_range'])
    gaps_cutoff = float(options['gaps_cutoff'])
    minimum_distance = float(options['dist'])
    
    print ('Flanked gaps range size: ' + str(blanks_site_range) + ' sites')
    print ('Minimum proportion of gaps: ' + str(gaps_cutoff))
    print ('Conserved region analysis size: ' + str(conserv_site_range) + ' sites')
    print ('Minimum distance between consensus and species sequence: ' + str(minimum_distance))
    
    sites_index = retrieve_selected_sites(path)
    if unique_sub == True:
        print ('Running unique substitution identification...')
        passed_unique_selected_sites = unique_selected_sites(path, sites_index, species)
        if blanks == True:
            print ('Running sites with flanking gaps removal...')
            blank_count = pull_sites(passed_unique_selected_sites, blanks_site_range)
            blank_passed_sites = count_blanks(path, blank_count, gaps_cutoff, species)
            if conservation == True:
                print ('Running sites within unconserved regions removal...')
                conserv_count = pull_sites(blank_passed_sites, conserv_site_range)
                conserved_passed_sites = conservation_analysis(path, conserv_count, minimum_distance, species)
                return final_output(conserved_passed_sites) #UCB
            else:
                return final_output(blank_passed_sites) #UB
        elif conservation == True:
            conserv_count = pull_sites(passed_unique_selected_sites, conserv_site_range)
            conserved_passed_sites = conservation_analysis(path, conserv_count, minimum_distance, species)
            return final_output(conserved_passed_sites) #UC
        else:
            return final_output(passed_unique_selected_sites) #U
    elif blanks == True:
        print ('Running sites with flanking gaps removal...')
        blank_count = pull_sites(sites_index, blanks_site_range)
        blank_passed_sites = count_blanks(path, blank_count, gaps_cutoff, species)
        if conservation == True:
            print ('Running sites within unconserved regions removal...')
            conserv_count = pull_sites(blank_passed_sites, conserv_site_range)
            conserved_passed_sites = conservation_analysis(path, conserv_count, minimum_distance, species)
            return final_output(conserved_passed_sites) #CB            
        else:
            return final_output(blank_passed_sites) #B
    elif conservation == True:
        print ('Running sites within unconserved regions removal...')
        conserv_count = pull_sites(sites_index, conserv_site_range)
        conserved_passed_sites = conservation_analysis(path, conserv_count, minimum_distance, species)
        return final_output(conserved_passed_sites) #C
    else:
        return 'ERROR: NO FUNCTIONS PROVIDED'



def final_output(dictionary):
    '''Take return of wrapper and export a CSV of outputs'''
    final_output = {}
    for k, v in dictionary.items():
        if v != '':
            final_output[k] = v
    df = pd.DataFrame.from_dict(final_output, orient = 'index')
    return df.to_csv('filtered_positively_selected_sites.tsv', sep='\t', index_label = 'Gene_family', header = ['Positively_selected_sites'])
    

def retrieve_selected_sites(path):
    '''Use alignment outputted from codeml_reader to give 
    co-ordinates of positively selected sites'''
    
    
    alignment_files = glob.glob(path + '*/PosSites*modelA.fasta')
    selected_sites = {}
    
    for alignment in alignment_files:
        
        selected_sequence = []
        positively_selected_sites = []
        for record in SeqIO.parse(alignment, 'fasta'):
            if 'PS_Characters|' in record.id: 
                for character in str(record.seq):
                    selected_sequence.append(character)
        for position, char in enumerate(selected_sequence):
            if char != '-':
                positively_selected_sites.append(position)
        #print (positively_selected_sites)
        aa_positively_selected_sites = positively_selected_sites[0::3] # Transform from nuc positions to amino acid
        myInt = 3
        aa_positively_selected_sites[:] = [int(x / myInt) for x in aa_positively_selected_sites]
        family = alignment.rpartition('/')[0].rpartition('/')[2]
        non_list = ', '.join(map(str, aa_positively_selected_sites))
        selected_sites[family] = non_list   
    return selected_sites



def unique_selected_sites(path, sites, species):
    '''From all sites under positive selection find only those that have amino acid differences
    unique to the species specified'''
    
    dictionary = {}
    
    for family, positions in sites.items():
        alignment = (path + family + '/PosSites_' + family + '_' + species + '_modelA.fasta')
          
        records = list(SeqIO.parse(alignment, "fasta")) 
        seq_dict = {} #dictionary of all sequences in amino acid form, seperated by species
        for i in range(1,len(records)):#skip selection data output sequences from VESPA output
            if 'PS_Characters' not in records[i].id:
                id = records[i].id.partition('|')[0]
                seq = str(records[i].seq)
                aa_seq = Seq(seq, Gapped(generic_dna, "-")).translate()
                seq_dict[id] = str(aa_seq)
        
        #Write list of positively selected sites in each species, see if Saiga is unique and all other conserved
        passed_sites = []
        list_positions = positions.split(",")
        #print (list_positions)
        for index in list_positions: #iterate through all positively selected sites
            index = int(index)
            species_character = []
            other_species = []
            for key, value in seq_dict.items():#if species' amino acid is not a gap add all species to dictionary 
                if key == species and value[index] == '-':
                    pass
                if key == species and value[index] != '-':
                    species_character.append(value[index])
                else:
                    other_species.append(value[index])
            if all(x == other_species[0] and x not in species_character or x == '-' or x == 'X' for x in other_species):
                passed_sites.append(str(index + 1))#Increase index number by one to account for python array

        passed_sites_unlist = ', '.join(passed_sites)
        dictionary[family] = passed_sites_unlist

        print (sites)
        print (passed_sites_unlist)
    return dictionary
  
    

def pull_sites(selected_sites, site_range):
    '''From an output file containing gene_families and selected site, pull out selected sites
    and defined flanking regions'''
    
    gene_dict = {}
    for gene_id, sites in selected_sites.items():
        sites_dict = {}
        for site in sites.split():
            site = site.replace(',','')
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
    blank_passed_sites = {}
    for family in gene_dict.keys():
        file = path + family + '/PosSites_' + family + '_' + species + '_modelA.fasta'  
        seq_dict = {}
        for record in SeqIO.parse(file, "fasta"):
            #print (species + '|')
            if record.id.startswith(species + '|'): #species input is case sensitive
                id = record.id.partition('|')[0]
                seq = str(record.seq)
                aa_seq = Seq(seq, Gapped(generic_dna, "-")).translate()
                seq_dict[id] = str(aa_seq)
  
        for gene, sel_site in gene_dict.items():
            if gene == family:
                sites = sel_site
        passed = []
        for sel_site, flank_sites in sites.items():
            py_index = [x-1 for x in flank_sites] #account for zero index
            region = []
            for index in py_index:
                for k, v in seq_dict.items():
                    if index < len(v):
                        region.append(v[index])

            if region.count('-') < (gaps_cutoff*len(py_index)):
                passed.append(sel_site)
#                print ('pass')
#            else:
#                print ('too gappy') #option to send unused sites to report
        passed_unlist = ', '.join(passed)
        blank_passed_sites[family] = passed_unlist
    return blank_passed_sites

def conservation_analysis(path, gene_dict, similarity_cutoff, species):
    
    conserv_passed_sites = {}
    conserv_rejected_sites = {}    
    
    for family in gene_dict.keys():
        passed_sites = []
        rejected_sites = []

        file = path + family + '/PosSites_' + family + '_' + species + '_modelA.fasta'  
       
        #Grab all the sites in the flanks of the selected sites 
        records = list(SeqIO.parse(file, "fasta"))
        seq_dict = {}
    
        for i in range(2,len(records)):
            id = records[i].id.partition('|')[0]
            seq = str(records[i].seq)
            aa_seq = Seq(seq, Gapped(generic_dna, "-")).translate()
            seq_dict[id] = str(aa_seq)
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
                #consensus.append(unique_chrs) #output all different amino acids present
            df.loc[len(df)] = consensus
            #print (df)
            #print (consensus)

            species_seq  = (df.iloc[:, :].sum(1)[species])
            consensus_seq = (df.iloc[:, :].sum(1)[-1])
            fasta_seq = '>' + species + '\n' + species_seq + '\n' + '>consenus' + '\n' + consensus_seq + '\n'
        
            with open('temp.txt', 'w') as interim_file:
                interim_file.write(fasta_seq)

            aln = AlignIO.read('temp.txt', 'fasta')
            os.remove('temp.txt')
        
            calculator = DistanceCalculator('blosum62')
            dm = calculator.get_distance(aln)

            if dm[0][1] < similarity_cutoff:
                passed_sites.append(sel_site)
            else:
                rejected_sites.append(sel_site)
        
        conserv_passed_sites[family] = ', '.join(passed_sites)
        conserv_rejected_sites[family] = ', '.join(rejected_sites)

    return conserv_passed_sites


def wrapper(path, species, functions, variables):
    '''Take command line input and deliver options for function set'''
    path = sys.argv[1]
    species  = sys.argv[2]
    functions = sys.argv[3]
    variables = sys.argv[4:]

    unique_sub = False
    blanks = False
    conservation = False    
    
    options = {}

    options['b_range'] = 20
    options['c_range'] = 5
    options['gaps_cutoff'] = 0.75 
    options['dist'] = 0.1
    
        
    if 'U' in functions:
        unique_sub = True
    if 'B' in functions:
        blanks = True
    if 'C' in functions:
        conservation = True
       

    for variable in variables:
        if 'b_range' in variable:
            options['b_range'] = variable.rpartition('=')[2]
        if 'c_range' in variable:
            options['c_range'] = variable.rpartition('=')[2]
        if 'gaps_cutoff' in variable:
            options['gaps_cutoff'] = variable.rpartition('=')[2]
        if 'dist' in variable:
            options['dist'] = variable.rpartition('=')[2]

    return function_set(path, species, options, unique_sub, blanks, conservation)
    



## EXAMPLE
## post_reader.py /Users/david_orr/Documents/Saiga/codeml_outputs... saiga BCV b_range=10 c_range=10 gaps_cutoff=0.8 dist=0.3

wrapper(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
