# absx
Alignment-Based Sequence eXtraction

## Dependencies
Python 2.7, Biopython

## Installation
Clone the repository like this:
```
https://github.com/AlexKnyshov/absx.git
```

## Workflow
### Basic extraction
Requirements:
- bait sequences in files, one file per locus (can contain multiple sequences)
- target files (typically, assembled contigs) in fasta format
- list of file names of target files (in case there is more than one target)

#### Search
Seach using tblastx:
```
bash blast_wrapper.sh ./query_folder/ ./assembly_folder/ 1e-05 tblastx 4 n ./list_of_files_to_seach_against.txt
```
After search is complete, the output (*.blast) files can be placed in a separate folder
```
mkdir blast_results
mv *.blast blast_results
```
#### Sequence extraction
Basic extraction can be done like this:
```
python absx.py -x b -f S -b ./blast_results/ -c 1 -e 1e-05 --is --ac tdna-tdna -m b-i-e -t ./assembly_folder/
```

### Reciprocal search
In case a reciprocal check is needed, 
####
The search can be done as follows:
```
bash reciprocal_search.sh ./blast_results/ ./assembly_folder/ ./Reference_assembly.fas tblastx 4 n ./reciprocal_get_contigs.py ./list_of_files_to_seach_against.txt
mkdir reciprocal_results
mv ./blast_results/*_reciprocal.blast reciprocal_results
```
In order to assess the best matches of queries to the reference assembly, it also needs to be searched against
#### Sequence extraction
Sequence extraction with consideration of reciprocal best hit can be done as follows:
```
python absx.py -x b -f S -b ./blast_results/ -c 1 -e 1e-05 --is --ac tdna-tdna -m b-i-e -t ./assembly_folder/ -r ./reciprocal_results -R Reference_assembly.fas.blast
```
### Hmmer usage
coming soon...

### Other features and parameter description
option `-f` specifies whether a single alignment table, or multiple tables are used. In the latter, only path to the folder needs be specified.

option `-b` specifies the path to alignment table file or a folder with such tables. When multiple tables are used, match between the table and the assembly is done based on cutting out `.blast` extension from the table and using the resulting name to search for the assembly file.

option `-x` specifies the extraction type: `n` extracts the whole contig, `s` extracts only single best hit, `a` extracts all hit regions and joins them together, `b`extracts region between two outmost hit regions.
                        
option `-t` specifies the path to the target file (when `-f S`) or to the folder with target files. Filename matching is done as described for `-b` option.

option `-q` specifies the path to the folder with query file(s) to which extracted results are to be appended

option `--om` specifies the way sequences are output. When set to `query`, sequences are grouped into per query files, when set to `target`, they are grouped into per target files, and when set to `combined`, all sequences from all tables are extracted into a single file.
                        
option `-e` specifies the evalue cutoff (default: 0.01). Nothing will be considered below this cutoff as it filters out initial alignment table parsing.

option `-c` specifies the max number of (super)contigs to extract. By default it is set to 0, which causes all supercontigs to be extracted. If a single (best matching) contig needs to be extracted, so `-c` should be set to 1. 

option `--fl` specifies flanks on each side in bp. This option is only available when `-x` is set to `s` or `a`

option `--lr` specifies local single best match check. When set to `range`, each region of the target contig (after joining multiple hits) is allowed to be matched to only one query. When set to `actual`, individual hits are checked for the same condition prior to being joined together. Can be switched off by setting `none`.
                     
option `--is` turns on contig stiching

option `--bt` specifies the alignment table type (only for forward searches; reciprocal and reference tables are always parsed as `blast`). `blast` is a standard blast table, `hmmer22` is a --domtblout table of hmmer, `hmmer18` is a protein --tblout table of hmmer, `hmmer15` is a dna --tblout table of hmmer. 

option `--ac` specifies the alignment table type. `dna-dna` is default, and has no special effects, except is not allowed to run with `--bt hmmer18` as latter is a protein table. In `tdna-tdna`, `tdna-aa` and `aa-tdna` overlapping hits are checked for frameshift before joining. When set to `tdna-aa` or `aa-tdna`, coordinates of blast tables are modified accordingly to convert between AA and NT values. In `tdna-aa` and `aa-aa` output sequence translation is not allow.

option `--translate` turns on sequence translation (for `-x s` or `-x a`), only works when appropriate `--ac`

option `--hit-ovlp` specifies max allowed hit overlap on query, in bp (default: 5). If two hits overlap more than this amount they are considered to be indeed overlapping.

option `--ctg-ovlp` specifies max allowed contig overlap on query, in bp (default: 1). If two contigs overlap more than this amount they are considered to be indeed overlapping.

option `--recip-ovlp` specifies max allowed hit/contig overlap on query for reciprocal check, in bp (default: 10). If two hits/contigs overlap more than this amount, they are considered to be sufficiently overlapping to pick only one best out of the two.

option `-r` specifies path to the reciprocal search output table file or the folder with such files. In case of multiple files, match is done based on suffix `_reciprocal.blast`.

option `-R` specifies path to the query search against the reference assembly.
                        
option `-m` specifies the order of metrics to be used for discriminating between hits/contigs (e - evalue, b - bitscore, i - identity) (default: e-b-i)

option `--rescale-metric` rescales the metric value by the length of the match

option `--ref-hs` turns on hit sticher on the reciprocal table (slow) (default: False). Typically reciprocal table is much larger, and takes considerable amount of time to parse. The alternative, default, approach is for a given target region to simply pick the best hit to reference that is located in the same region.

