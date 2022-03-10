# Fetch-NCBI-genomes

**NCBI-utils** provides a set of functions to perform cross-search between NCBI databases and retrieve information and/or sequences for target species. So far the **fetch_ncbiGenomes.py** retrieves the reference sequence (NCBI-Refseq) for a single specie or multiple species based on their common or scientific name or a list of accession numbers [See NCBI-nucleotide](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/).

## INPUT
The input can be a single string defining the name of the specie (common or scientific name) or a configuration file indicating the scientifc name or the sequence accession numbers for the target specie. See [input templates](https://github.com/Daniel-VM/ncbi-utils/tree/main/conf)

*CONFIGURATION FILE*
```console
File format:
                Bovine herpesvirus 1
                Bovine herpesvirus 5
                BK polyomavirus
                Bovine leukemia virus
                
                ======= or =========
                
                NC_055230.1
                NC_055330.1
                NC_055331.1
                NC_055332.1
```

## PARAMETERS

`--input` / `-i`:  String indicaring the name of the target specie: --input "Bovine herpesvirus 1". 

`--config` / `-c`: Configuration file containing the names of the species.

`--accession` / `-ac`: Configuration file containing the sequence's accession numbers.

`--outdir` / `-o`: Directory to place the results.

`--email` / `-e`:  User's email (mandatory).

#### `--help / -h`: Help page.
```console
python  fetch_ncbiGenomes.py --help
```

## Usage and output

Depending on whether you search the genome of a single or multiple organisms two different command line inputs: 

1. **Single target organism**

The typical command :
```console
python  fetch_ncbiGenomes.py -i "Bovine herpesvirus 1" \
                             -e example@domain.com \
                             -o path_outDir/

```

>  genome.fa: Fasta file containig the DNA sequences of the target specie.

1. **Multiple target organisms:**
The typical command:
```console
python  fetch_ncbiGenomes.py -c conf/species.config \
                             -e example@domain.com \
                             -o path_outDir/

```

> genome.fa: multifasta file.
> report_ncbiGenomes.csv: Dataframe summarizing the process.


:warning: Note: The Bio.Entrez module in python requieres an email to use this resource. See [Bio.Entrez](https://biopython.org/docs/1.76/api/Bio.Entrez.html).
