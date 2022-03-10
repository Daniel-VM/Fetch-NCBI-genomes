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
from ncbi_utils import RefGenome_records
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
parser.add_argument('-c', '--config', type = str, help = 'Input is single column file containing the specie common or scientifc name.')
parser.add_argument('-ac', '--accession_ids', type = str,  help = 'Input is a file containing a single column of accession identifiers.')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-e', '--email', required = True, type = str, help = 'User email')
args = vars(parser.parse_args())

# Set email
Entrez.email = args['email']

# ===============================================
#   Validate User's Input / Output
# ===============================================

if args['input'] is not None and args['config'] is not None:
    sys.exit(
        '[ERROR] More than one input type has been provided: -i "{}" and -c "{}".\n'
        .format(args['input'], args['config'])
    )

if args['outdir'] is None:
    args['outdir'] = os.getcwd()
outFile = os.path.join(args['outdir'], 'genome.fa')

# ===============================================
#   Search for the reference genome - NCBI RefSeq
# ===============================================
# Save the summary of the process in var(data)
data = collections.defaultdict(list)

# Looking for the refernece genome o f a single specie: 
if args['input'] is not None:
    try:
        ncbi_record = RefGenome_records(taxon = args['input'], email = args['email'])
    except SystemExit:
        sys.exit("[ERROR] No taxonomy ID was found for '{}'. \nAborting..."
                 .format(str(args['input']))
        )
    with open(outFile, 'w') as out_fasta:
            nucleotide_handle = Entrez.efetch(db="nucleotide", id = ncbi_record['records'], retmax="200", rettype="fasta", retmode="text")
            out_fasta.write(nucleotide_handle.read())
            nucleotide_handle.close()

# Looking for refernece genomes of many species by common or scientific name ( User's input is a config file with the scientific names of the species). 
if args['config'] is not None:
    with open(outFile, 'w') as out_fasta:
        with open(args['config'],'r') as config_file:
            codes = config_file.readlines()
            for code in codes:
                tax_name = code.strip()
                print('Handling "{}"...'.format(str(tax_name)))
    
                # Get NCBI records of specie's reference genomes
                try:
                    ncbi_record = RefGenome_records( taxon = tax_name, email = args['email']) 
                except SystemExit:
                    print("Still not found. Droping {} ".format(str(tax_name)))
                    ncbi_record = {
                        'taxon_name': tax_name,
                        'taxon_id': [NaN],
                        'records': [NaN],
                        'complete_genome': False
                        }    

                # If not reference genome is available 
                if ncbi_record['complete_genome'] == False:
                    data['taxon_name'].append(ncbi_record['taxon_name'])
                    data['taxon_id'].append(ncbi_record['taxon_id'])
                    data['accession'].append(NaN)
                    data['accession_title'].append(NaN)
                    data['complete_genome'].append(ncbi_record['complete_genome'])
                    continue

                # Retrieve complete genome info and fasta sequence
                handle_fasta = Entrez.efetch(db="nucleotide", id = ncbi_record['records'], retmax="200", rettype="fasta", retmode="text")
                genome_sequence = handle_fasta.read()
                out_fasta.write(genome_sequence)
                handle_fasta.close()

                # Get genome accession identifiers
                accession = re.search(">(.*?) ",genome_sequence).group(1)
                try:
                    accession_title =  re.search(
                                                ">{} (.*?),".format(str(accession)), genome_sequence
                                                ).group(1)
                except:
                    accession_title =  re.search(
                                                ">{} (.*?) complete".format(str(accession)), genome_sequence
                                                ).group(1)

                # Save process summary
                data['taxon_name'].append(ncbi_record['taxon_name'])
                data['taxon_id'].append(ncbi_record['taxon_id'])
                data['accession'].append(accession)
                data['accession_title'].append(accession_title)
                data['complete_genome'].append(ncbi_record['complete_genome'])


# Looking for refernece genomes of many species by NCBI-nucliotide accession number ( User's input is a config file with accession numbers). 
if args['accession_ids'] is not None:
    with open(outFile, "a+") as out_fasta:
        with open(args['accession_ids'], 'r') as config_file:
            codes = config_file.readlines()       
            for code in codes:
                accession_id = code.strip()
                print('Handling "{}"...'.format(str(accession_id)))

                # Check whether the accession number has been previously loaded to the output
                if accession_id in data['accession']:
                    print('The accession "{}" has been previously processed... Proceeding with the next item.'
                            .format(str(accession_id)))
                    continue

                # Retrieve genome sequence
                handle_fasta = Entrez.efetch(db="nucleotide", id = accession_id, retmax="200", rettype="fasta", retmode="text")
                out_fasta.write(handle_fasta.read())
                handle_fasta.close()
                
                # Retrieve complet report of accession_id from NCBI-nucleotide database
                try:
                    esummary_handle = Entrez.esummary(db = "nucleotide", id = accession_id, report = 'full')
                    esummary_record = Entrez.read(esummary_handle)
                    esummary_handle.close()
                except:
                    print('[ERROR] There is a problem with {}'.format(str(accession_id)))
                    data['taxon_name'].append('Accession Error')
                    data['taxon_id'].append(NaN)
                    data['accession_title'].append(NaN)
                    data['accession'].append(accession_id)
                    data['complete_genome'].append(False)
                    esummary_handle.close()
                    continue

                # Get species taxon name
                taxonomy_search_term = int(esummary_record[0]['TaxId'])
                try:              
                    taxonomy_handle = Entrez.esummary(db='taxonomy', id = taxonomy_search_term, rettype='gb', retmode='text')
                    taxonomy_record = Entrez.read(taxonomy_handle, validate=False)
                    tax_name = taxonomy_record[0]['ScientificName']
                    taxonomy_handle.close()
                except:
                    tax_name = '[ERROR] Taxid error'

                # Gather variables that will be printed in the summary
                tax_id = int(esummary_record[0]['TaxId'])
                record_title =  esummary_record[0]['Title']
                ncbi_accession = esummary_record[0]['AccessionVersion']

                # Save process summary
                data['taxon_name'].append( tax_name )
                data['taxon_id'].append( tax_id )
                data['accession_title'].append( record_title )
                data['accession'].append( ncbi_accession )
                if 'complete genome' in record_title:
                    data['complete_genome'].append(True)
                else:
                    data['complete_genome'].append(False)
                
# Transform summary dict to dataframe and save it
df = pd.DataFrame(data)
df.to_csv(os.path.join(args['outdir'], 'report_ncbiGenomes.csv'), index=False)
print(df)