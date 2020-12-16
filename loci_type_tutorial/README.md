# ALiBaSeq tutorial - different types of loci
Based on *R. prolixus* genome

## Dependencies
Python 2 or 3, Biopython, NCBI BLAST

## How to run the test

Navigate to the test folder.
If running on a local PC, simply execute this command:
```
bash run_test.sh
```

If running on a cluster system with software available through `module load` commands, edit the `run_test.sh` to add the appropriate `module load` command, or execute the following command (assuming the BLAST module is available under `ncbi-blast`):
```
echo "module load ncbi-blast" | cat - run_test.sh | bash
```

## Contents of the test

Both baits and the sample represent the sequences from the *R. prolixus* genome. While it is possible to simply run the blastn for the test, discontinuous megablast is used to highlight paralog presence for some of the loci.
* First, the sample file is unzipped, the blast database is created
* Then, blast is run using `baits.fas` file as the query against the database
* Then, ALiBaSeq is run using 'extract single best match' option; results have `best` suffix or prefix
* Then, ALiBaSeq is run using 'extract all matches' option; results have `all` suffix or prefix
Despite the `baits` folder is not needed for the search (it could be used for search instead of the single bait file via the blast helper script), in order to demostrate ALiBaSeq's primary function to compile a phylogenetic dataset and to ease the comparison between baits and results, the test results will be output one per file with the corresponding bait prepended. Output `FASTA` files can be fed to an multiple sequence aligner of the choice to show aligned sequences.

## Information on test baits

We included few baits for loci of variable complexity:
* SingleExon - a bait from a single exon of a gene without close paralogs
* MultiExon - a bait derived from several spliced exons of a gene without close paralogs
* MultiContig - a bait derived from several spliced exons of a gene without close paralogs; in the sample the gene is represented by two contigs, broken up in one of the intronic regions of the gene
* ParalogsOnDifContig - a bait derived from several spliced exons of a gene with one similar paralog, located on a different contig
* ParalogsOnSameContig - a bait derived from several spliced exons of a gene with several similar paralogs, located sequentially on the same contig