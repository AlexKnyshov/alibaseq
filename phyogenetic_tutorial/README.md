# ALiBaSeq tutorial - phylogenetic dataset

### Table of Contents
* [Dependencies](https://github.com/AlexKnyshov/alibaseq/phyogenetic_tutorial#dependencies)  
* [Description](https://github.com/AlexKnyshov/alibaseq/phyogenetic_tutorial#description)  
* [Tutorial](https://github.com/AlexKnyshov/alibaseq/phyogenetic_tutorial#tutorial)  

## Dependencies
Python 2 or 3, Biopython, NCBI BLAST

## Description

Four UCE loci of the Hemiptera 2.7K version 1 kit are being retrieved from 3 samples: a WGS sample (based on *Rhodnius prolixus* genome), an RNA-seq sample (based on *Lygus hesperus* transcriptome), and a hybrid capture sample (based on *Rhodnius robustus* UCE capture). Retrieved sequences are to be appended to the bait sequences used for search, which were derived from the fourth genome (based on *Cimex lectularius* genome). Resulting FASTA files can be used for downstream steps (alignment, filtering / trimming, phylogenetic reconstruction).

## Tutorial

We assume user's working directory is `alibaseq/phyogenetic_tutorial` and consequently, `reference`, `assemblies`, and `baits` folders are in the current directory, while scripts used are in the folder above the working directory. BLAST and Python executables are assumed to be in path. BLAST multithreading is set to 4.

First, create a BLAST database for the reference genome:
```
makeblastdb -in reference/reference_genome.fasta -dbtype nucl -parse_seqids
```

Next, search the baits (derived from the reference genome) against reference genome (needed for RBH check)
```
mkdir refblast_results
bash ../blast_wrapper.sh baits/ reference/reference_genome.fasta 1e-10 blastn 4 n
mv *.blast refblast_results
```

Create a list of assemblies to blast (needed for use with provided blast and reciprocal search wrapper scripts):
```

ls assemblies/ > assembly_list.txt

```

Create BLAST databases for sample assemblies:
```
for f in assemblies/*.fasta
do
	makeblastdb -in $f -dbtype nucl -parse_seqids
done
```

Search for baits in sample assemblies:
```
mkdir blast_results
bash ../blast_wrapper.sh baits/ assemblies/ 1e-10 dc-megablast 4 n assembly_list.txt
mv *.blast blast_results
```

Search contigs found in baits-vs-sample search against the reference genome:
```
bash ../reciprocal_search.sh blast_results assemblies \
reference/reference_genome.fasta dc-megablast 4 n ../reciprocal_get_contigs.py \
assembly_list.txt
```

Run alibaseq (Python 2 `python` executable is assumed): 
```
python ../alibaseq.py -x a -f M -b blast_results -t assemblies -q baits -s tutorial_logs -o tutorial_results \
-e 1e-10 --is --amalgamate-hits -r blast_results -R refblast_results/reference_genome.fasta.blast
```
If Python 3 version is needed, replace `../alibaseq.py` with `../alibaseqPy3.py`, if Python 3 executable is `python3`, use that instead of `python`.

Final output is in the `tutorial_results` folder.