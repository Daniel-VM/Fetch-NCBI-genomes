#!/usr/bin/env python3

# ===============================================================
#   Set of functions to deal with NCBI data
# ===============================================================
# File name: ncbi_utils.py
# Author: Daniel Valle Millares
# Date created: 10/03/2022
# Date last modified: 10/03/2022
# Python Version: 3.8.8
# 
# Description: 
#       
#       These funtions ease to obtain and manage NCBI data, like 
#       genes, CDS, partial sequences, complete genomes...
#       
#       They employ bio.Entrez module to carry out these task by performing 
#       cross search between the NCBI databases. In addition, they provide 
#       a summary of the process.

# Performs cross-search between NCBI databases (NCBI-taxonomy & NCBI-nucleotide) and retrieves the NCBI-nucleotide record of complete genome/s given a target specie.
def RefGenome_records(taxon, email):
    '''
    Uses the scientific or common name [taxon] of an specie and looks for its NCBI taxonomy ID
    and performs cross-search between NCBI-taxonomy and NCBI-nucleotide database to retrieve the 
    record that will be used afterwards to get the complete genome of the taxon. 
    '''    
    # Load modules
    from pickle import FALSE
    from Bio import Entrez
    from numpy import NaN
    import re
    import sys
    
    Entrez.email = email
    
    # From taxon name to tax id (source: https://harryincupboard.blog/3)
    if not re.match(r'\d+', taxon):
        # Get taxonomy ID using Entrez
        taxon2 = '"' + taxon + '"'
        handle_taxonomy = Entrez.esearch(
            db='taxonomy', term=taxon2, rettype='gb', retmode='text')
        record_taxonomy = Entrez.read(handle_taxonomy, validate=False)
        handle_taxonomy.close()
        # If there's no result
        if not record_taxonomy['IdList']:
            sys.exit(
                '[ERROR] The taxon "{}" you provided is invalid. '
                'Please check NCBI Taxonomy'.format(taxon))
        tax_id = record_taxonomy['IdList']
    else:
        tax_id = taxon

    # Now connect to NCBI again and perform a cross search by using the taxonomy identifier and retrieve the record_id of the complete genome.
    search_term = 'txid' + str(tax_id[0]) + '[Organism] AND refseq[filter] AND complete genome[title]'
    handle_nucleotide = Entrez.esearch(db="nucleotide", term = search_term , retmax="200", rettype="fasta", retmode="text")
    record_nucleotide = Entrez.read(handle_nucleotide, validate=False)

    
    # Check whether ncbi search succed
    if not record_nucleotide['IdList']:
        id_record = NaN
        iscomplete = False
    else:
        id_record = record_nucleotide['IdList']
        iscomplete = True
    handle_nucleotide.close()

    dict = {
            'taxon_name': taxon,
            'taxon_id': tax_id,
            'records': id_record,
            'complete_genome': iscomplete
            }
    return(dict)