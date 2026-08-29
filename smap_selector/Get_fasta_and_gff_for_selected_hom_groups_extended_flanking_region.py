#!usr/bin/python3

#===============================================================================
# Description
#===============================================================================

# Sabine VAN GLABEKE March 2021

"""
This script creates input fasta and gff files for script to design amplicons
Required input:
- Genome data in gff format (tab-delimited)
- Genome data in fasta format
- Gene family information for the (coding) genes, separated per gene family type (tab-delimited)
- Species to extract, corresponding with species indicated in gene family info file
- Requested border region extending both sides of the gene in bp (default 0 or enter a positive value)
- List with selected homology groups or list with selected genes
"""

#===============================================================================
# Import modules
#===============================================================================

import os, sys, argparse
from datetime import datetime
import gffutils
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import re

#===============================================================================
# Parse arguments
#===============================================================================

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Get sequences of genes in selected homology groups in fasta format and gff format. List with selected homology groups or list with selected genes can be given as input')

# Add positional arguments (mandatory)
parser.add_argument('gff',
                    help = 'Genome data in gff format (tab-delimited)')

parser.add_argument('fasta',
                    help = 'Genome data in fasta format')

parser.add_argument('gene_family_info',
                    help = 'Gene family information for the (coding) genes, separated per gene family type (tab-delimited)')

parser.add_argument('species',
                    help = 'Species to extract, corresponding with species indicated in gene family info file')

# Add optional arguments
parser.add_argument('-r','--region',
                    default = 0,
                    type = int,
                    help = 'Flanking region extending both sides of the gene in bp (default: 0 or enter a positive value)')

parser.add_argument('-f', '--hom_groups',
                    default = '',
                    type = str,
                    help = 'List with selected homology groups (default: empty and given list with selected genes is used)')

parser.add_argument('-g', '--genes',
                    default = '',
                    type = str,
                    help = 'List with selected genes (default: empty and given list with selected homology groups is used)')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

#===============================================================================
# Functions
#===============================================================================

def print_date ():
    """
    Print the current date and time to stderr.
    """
    sys.stderr.write('----------------\n')
    sys.stderr.write('{}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M')))
    sys.stderr.write('----------------\n\n')
    return

def create_hom_group_list_from_gene_list(genes, gene_family_info):
    """
    Create an unique homology group list from the given gene list.
    """
    genes_in = open(genes).read().split()
    gene_family_info_dataframe = pd.read_csv(gene_family_info, sep = '\t', header = None, low_memory=False, comment='#', names=["hom_group", "species", "gene"])
    hom_groups_selected_dataframe = gene_family_info_dataframe[(gene_family_info_dataframe['gene'].isin(genes_in))]
    genes_not_present = set(genes_in) - set(hom_groups_selected_dataframe['gene'])
    genes_not_present_length = len(genes_not_present)
    genes_not_present_path = os.path.splitext(genes)[0] + '_not_present.txt'
    genes_not_present_out = open(genes_not_present_path, 'w')
    for gene_not_present in genes_not_present:
        genes_not_present_out.write(gene_not_present + "\n")
    genes_not_present_out.close()
    hom_groups_present_path = os.path.splitext(genes)[0] + '_selected_homology_groups.txt'
    hom_groups_present_out = open(hom_groups_present_path, 'w')
    for gene in set(hom_groups_selected_dataframe['hom_group']):
        hom_groups_present_out.write(gene + "\n")
    hom_groups_present_out.close()
    sys.stderr.write('* {} gene(s) are not present in {} and are listed in {}\n\n'.format(genes_not_present_length, gene_family_info, genes_not_present_path))
    return hom_groups_present_path

def create_gene_list_per_hom_group(hom_list, gene_family_info, species):
    """
    Create an unique gene list from homology groups.
    """
    hom_list_in = open(hom_list).read().split()
    gene_family_info_dataframe = pd.read_csv(gene_family_info, sep = '\t', header = None, low_memory=False, comment='#', names=["hom_group", "species", "gene"])
    hom_list_not_present_path = os.path.splitext(hom_list)[0] + '_not_present.txt'
    hom_list_not_present_out = open(hom_list_not_present_path, 'w')
    hom_group_dataframe_per_species = gene_family_info_dataframe[(gene_family_info_dataframe['hom_group'].isin(hom_list_in)) & (gene_family_info_dataframe['species'] == species)]
    gene_list_per_hom_group = hom_group_dataframe_per_species.groupby('hom_group')['gene'].apply(list).reset_index(name='hom_group_list')
    hom_list_present = gene_list_per_hom_group['hom_group'].tolist()
    hom_list_not_present = set(hom_list_in) - set(hom_list_present)
    hom_list_not_present_length = len(hom_list_not_present)
    for hom_group in hom_list_not_present:
        hom_list_not_present_out.write(hom_group + "\n")
    hom_list_not_present_out.close()
    sys.stderr.write('* {} homology group(s) are not present in {} for species "{}" and are listed in {}\n\n'.format(hom_list_not_present_length, gene_family_info, species, hom_list_not_present_path))
    return gene_list_per_hom_group

def create_gff_db(gff):
    """
    Create a database from the gff file.
    """
    gff_db = gffutils.create_db(gff, ':memory:', merge_strategy="create_unique", keep_order=True, force=True)
    return gff_db

def get_sequences_of_genes_in_selected_hom_groups_from_genome_fasta_using_gff_file(gff, fasta, hom_list, gene_family_info, species, region):
    """
    Get sequences of genes in selected homology groups in fasta format.
    """
    gene_list_per_hom_group = create_gene_list_per_hom_group(hom_list, gene_family_info, species)
    hom_list_present = gene_list_per_hom_group['hom_group'].tolist()
    type='gene'
    if region < 0:
        raise ValueError('Extended flanking region must be an positive value!')
        sys.exit()
    elif region > 0:
        gff_extended = create_gff_file_with_extended_region(gff, region, fasta)
        gff_db = create_gff_db(gff_extended)
    else:
        gff_db = create_gff_db(gff)
    for hom_group in hom_list_present:
        fasta_out_path = os.path.split(fasta)[0] + '/' + str(hom_group) + '_' + str(species) + '_selected_genes_extended_' + str(region) + '_bp.fasta'
        fasta_out = open(fasta_out_path, 'w')
        gene_list = gene_list_per_hom_group['hom_group_list'][gene_list_per_hom_group['hom_group'] == hom_group].item()
        for line in gff_db.features_of_type(type):
            if line.id in gene_list:
                line_seq = line.sequence(fasta)
                line_seq = Seq(line_seq)
                fasta_out.write('>' + str(line.id) + '\n' + str(line_seq) + '\n')
        fasta_out.close()
    return

def create_corresponding_gff_file(gff, hom_list, gene_family_info, species, region, fasta):
    """
    Get features of genes in selected homology groups in gff format.
    """
    gene_list_per_hom_group = create_gene_list_per_hom_group(hom_list, gene_family_info, species)
    hom_list_present = gene_list_per_hom_group['hom_group'].tolist()
    if region < 0:
        raise ValueError('Extended flanking region must be an positive value!')
        sys.exit()
    elif region > 0:
        gff_extended = create_gff_file_with_extended_region(gff, region, fasta)
        gff_db = create_gff_db(gff_extended)
    else:
        gff_db = create_gff_db(gff)
    for hom_group in hom_list_present:
        gff_in = open(gff)
        gff_out_path = os.path.split(gff)[0] + '/' + str(hom_group) + '_' + str(species) + '_selected_genes_extended_' + str(region) + '_bp.gff'
        gff_out = open(gff_out_path, 'w')
        gene_list = gene_list_per_hom_group['hom_group_list'][gene_list_per_hom_group['hom_group'] == hom_group].item()
        for line in gff_in:
            if line.startswith('#') or line.startswith('\n'):
                continue
            seqid, source, type, start, end, score, strand, phase, attribute = line.strip().split('\t')
            for gene in gene_list:
                if re.search(str(gene)+'[^0-9]+', attribute):
                    gene_ID = gene
                    if strand == '+':
                        new_start = int(start) - int(gff_db[gene_ID].start) + 1
                        new_end = int(end) - int(gff_db[gene_ID].start) + 1
                        gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene,source,type,new_start,new_end,score,strand,phase,attribute))
                    elif strand == '-':
                        new_start = abs(int(end) - int(gff_db[gene_ID].end)) + 1
                        new_end = abs(int(start) - int(gff_db[gene_ID].end)) + 1
                        gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t+\t{}\t{}\n".format(gene,source,type,new_start,new_end,score,phase,attribute))
                    else:
                        continue
        gff_out.close()
        gff_in.close()
    return

def create_gff_file_with_extended_region(gff, region, fasta):
    """
    Create a gff file with extended region.
    """
    gff_in = open(gff)
    gff_out_path = os.path.splitext(gff)[0] + '_extended_' + str(region) + '_bp.gff'
    gff_out = open(gff_out_path, 'w')
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    for line in gff_in:
        if line.startswith('#') or line.startswith('\n'):
            gff_out.write(line)
        else:
            seqid, source, type, start, end, score, strand, phase, attribute = line.strip().split('\t')
            if int(start) - int(region) < 0:
                new_start = 1
            else:
                new_start = int(start) - int(region)
            if (int(end) + int(region)) > len(fasta_dict[seqid].seq):
                new_end = len(fasta_dict[seqid].seq)
            else:
                new_end = int(end) + int(region)
            gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seqid,source,type,new_start,new_end,score,strand,phase,attribute))
    gff_in.close()
    gff_out.close()
    return gff_out_path

#===============================================================================
# Script
#===============================================================================

if __name__ == '__main__':
    
    print_date()
    
    if args['hom_groups'] != '':
        sys.stderr.write('* Get sequences of genes in selected homology groups in fasta format...\n\n')
        get_sequences_of_genes_in_selected_hom_groups_from_genome_fasta_using_gff_file(args['gff'], args['fasta'], args['hom_groups'], args['gene_family_info'], args['species'], args['region'])
        sys.stderr.write('* Get features of genes in selected homology groups in gff format...\n\n')
        create_corresponding_gff_file(args['gff'], args['hom_groups'], args['gene_family_info'], args['species'], args['region'], args['fasta'])
        sys.stderr.write('* Finished\n\n')
    elif args['genes'] != '':
        sys.stderr.write('* Get sequences of genes in homology groups in fasta format...\n\n')
        hom_list_unique = create_hom_group_list_from_gene_list(args['genes'], args['gene_family_info'])
        get_sequences_of_genes_in_selected_hom_groups_from_genome_fasta_using_gff_file(args['gff'], args['fasta'], hom_list_unique, args['gene_family_info'], args['species'], args['region'])
        sys.stderr.write('* Get features of genes in homology groups in gff format...\n\n')
        create_corresponding_gff_file(args['gff'], hom_list_unique, args['gene_family_info'], args['species'], args['region'], args['fasta'])
        sys.stderr.write('* Finished\n\n')
    else:
        sys.stderr.write('* Missing argument hom_groups or genes to run the script...\n\n')
        sys.exit()
    print_date()
