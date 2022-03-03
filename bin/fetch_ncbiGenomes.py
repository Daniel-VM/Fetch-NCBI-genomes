#!/usr/bin/env python3

# ========================================================
#      Fetching NCBI complete genomes
# ========================================================
# File name: fetch_ncbiGenomes.py
# Author: Daniel Valle Millares
# Date created: 03/03/2022
# Date last modified: 03/03/2022
# Python Version: 3.8.8
# 
# Description: 
#       
#       Scans the NCBI-nucleotide database to find the complete genomes of an organism.
#       To this end, it performs cross-search between NCBI databases, 
#       download complete genomes, and outputs a report of the process. 
#       The input can be a single string defining the name of the specie (common or scientific name)
#       you look for or a configuration file (tab separated format, see below) containing 
#       multiple species. An extensive information (usage, arguments ...) can be found in README.md.

# Load modules
from pickle import FALSE
from Bio import Entrez
from numpy import NaN
import pandas as pd
import argparse
import collections
import sys
import re
import os

# parsing arguments
parser = argparse.ArgumentParser(description="Fetch NCBI complete genomes...")
parser.add_argument('-i', '--input', type = str,  help = 'Specie name between quotes: --input "Bovine herpesvirus 1"')
parser.add_argument('-o', '--outdir', type = str, help = 'Directory to place the results')
requiredNamed = parser.add_argument_group('required arguments')
parser.add_argument('-c', '--config', type = str, help = 'Configuration file in tab separated format. Column1: specie three-code letter, Column2: specie complete name')
requiredNamed.add_argument('-e', '--email', required = True, type = str, help = 'User email')
args = vars(parser.parse_args())

# Set email
Entrez.email = args['email']

# This function performs cross-search between NCBI databases and retrieves the NCBI-nucleotide recird record of complete genome.
def get_ncbi_records(taxon):
    
    '''Getn NCBI taxonomy'''
    # From specie name to tax id (source: https://harryincupboard.blog/3)
    if not re.match(r'\d+', taxon):
        # Get taxonomy ID using Entrez
        taxon2 = '"' + taxon + '"'
        handle = Entrez.esearch(
            db='taxonomy', term=taxon2, rettype='gb', retmode='text')
        record = Entrez.read(handle, validate=False)
        handle.close()
        # If there's no result
        if not record['IdList']:
            sys.exit(
                '[ERROR] The taxon "{}" you provided is invalid. '
                'Please check NCBI Taxonomy'.format(taxon))
        tax_id = record['IdList']
    else:
        tax_id = taxon

    # Now connect NCBI again and perform a cross search by using the taxonomiy identifier. It returns the record_id for the complete genome.
    search_term = 'txid' + str(tax_id[0]) + '[Organism] AND refseq[filter]' 
    handle2 = Entrez.esearch(db="nucleotide", term = search_term , retmax="200", rettype="fasta", retmode="text")
    record2 = Entrez.read(handle2, validate=False)
    
    # Check whether ncbi search succed
    if not record2['IdList']:
        id_record = NaN
        iscomplete = False
    else:
        id_record = record2['IdList']
        iscomplete = True
    handle2.close()

    dict = {
            'taxon_name': taxon,
            'taxon_id': tax_id,
            'records': id_record,
            'complete_genome': iscomplete
            }
    return(dict)

# Verify user's input/output
if args['input'] is not None and args['config'] is not None:
    sys.exit(
        '[ERROR] More than one input type has been provided: -i "{}" and -c "{}".\n'
        .format(args['input'], args['config'])
    )

if args['outdir'] is None:
    args['outdir'] = os.getcwd()
outFile = os.path.join(args['outdir'], 'genome.fa')

# Here is where the process starts...

# Catching user's input type:
if args['input'] is not None:
    # Search for NCBI records & download complete genome. 
    try:
        ncbi_record = get_ncbi_records(args['input'])
    except SystemExit:
        sys.exit("[ERROR] No taxonomy ID was found for '{}'. \nAborting..."
                 .format(str(args['input']))
        )

    with open(outFile, 'w') as out_fasta:
            handle3 = Entrez.efetch(db="nucleotide", id = ncbi_record['records'], retmax="200", rettype="fasta", retmode="text")         
            out_fasta.write(handle3.read())
            handle3.close()

# Searching for each specie that has been defined in the user's config file
if args['config'] is not None:
    # Report of this process will be keep into var(data)
    data = collections.defaultdict(list)

    # Read user's config file
    with open(args['config'],'r') as config_file:
        codes = config_file.readlines()

        with open(outFile, 'w') as out_fasta:
            for code in codes:
                tax_acroname = code.strip().split('\t')[0]
                tax_name = code.strip().split('\t')[1]

                # Search for NCBI records and download NCBI genomes. to find/fetch specie's genome by using the complete name of the specie. If it fails, try to search again but using the specie name in three-letter code format.
                try:
                    try:
                        ncbi_record = get_ncbi_records(tax_name) 
                    except SystemExit:
                        print("[ERROR] No taxonomy ID was found for '{}'. Trying now with '{}'".format(str(tax_name), str(tax_acroname)))
                        ncbi_record = get_ncbi_records(tax_acroname) 
                except SystemExit:
                    print("Still not found. Droping {} ".format(str(tax_name)))
                    ncbi_record = {
                        'taxon_name': tax_name,
                        'taxon_id': [NaN],
                        'records': [NaN],
                        'complete_genome': False
                        }            

                # Keep the output of each ENREZ query
                data['taxon_name'].append(ncbi_record['taxon_name'])
                data['taxon_id'].append(ncbi_record['taxon_id'])
                data['records'].append(ncbi_record['records'])
                data['complete_genome'].append(ncbi_record['complete_genome'])

                # Now connect for the last time to NCBI and download organism complete genome based on NCBI's nucleotideDB record
                if ncbi_record['complete_genome'] == False:
                    continue
                handle3 = Entrez.efetch(db="nucleotide", id = ncbi_record['records'], retmax="200", rettype="fasta", retmode="text")
                out_fasta.write(handle3.read())
            
            # Close Entrex connection    
            handle3.close()

    # Transform summary dict to dataframe and save it
    df = pd.DataFrame(data)
    df.to_csv(os.path.join(args['outdir'], 'report_ncbiGenomes.csv'), index=False)
    print(df)

# DONE