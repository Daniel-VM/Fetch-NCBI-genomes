# Fetch-NCBI-genomes

**fetch_ncbiGenomes.py** searches the NCBI database for the complete genome/s of an organism. To this end, it performs cross-search at NCBI, identifies the target genome and downloads it. If multiple organisms are searched, then this script provides a report summarizing the process. 

## INPUT
The input can be a single string defining the name of the specie (common or scientific name) or a configuration file (tab separated format, see below) containing several species according to this schema (See [conf/species.config]()):

*CONFIGURATION FILE*
```console
File format:    [acronym name]\t[common or scientific name ]

Example Config_species.txt :
                bhv1	Bovine herpesvirus 1
                bhv5	Bovine herpesvirus 5
                bkv     BK polyomavirus
                blv     Bovine leukemia virus
```

## PARAMETERS

#### `--input` / `-i`:  A single organism name: --input "Bovine herpesvirus 1". 
#### `--config` / `-c`: Configuration file containing multiple speces.
#### `--outdir` / `-o`: Directory to place the results
#### `--email` / `-e`:  User's email

#### `--help / -h`: Help page
```console
python  download_viralGenomes.py --help
```

## Usage and output

Depending on whether you search the genome of a single or multiple organisms (by using a configuraton file) two different escenarios are possible: 


1. **Single target organism**

The typical command :
```console
python  download_viralGenomes.py -i "Bovine herpesvirus 1" \
                                 -e uuuu@ddddd.xx \
                                 -o path_outDir/

```
>  genomes.fa: Fasta file containig the DNA sequences of the target specie

1. **Multiple target organisms:**

The typical command:
```console
python  download_viralGenomes.py -i "Bovine herpesvirus 1" \
                                 -e uuuu@ddddd.xx \
                                 -o path_outDir/

```

> genomes.fa: multifasta file.
> report_ncbiGenomes.csv: Dataframe summarizing the process.


:warning: Note: The Bio.Entrez module in python requieres an email to use this resource. See [Bio.Entrez](https://biopython.org/docs/1.76/api/Bio.Entrez.html).
