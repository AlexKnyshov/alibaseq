echo "start"
if [[ $1 == "" ]]
then
	echo "########################################"
	echo "arguments"
	echo "\$1 folder with blast / hmm files"
	echo "\$2 path to assemblies"
	echo "\$3 path to reference database"
	echo "\$4 command: tblastx, blastn, etc."
	echo "\$5 cpu: 4"
	echo "\$6 cluster: y / n"
	echo "\$7 path to reciprocal_get_contigs.py"
	echo "\$8 list (for multiple blast / hmm files)"
	echo "########################################"
else 
	if [[ $6 == y ]]
	then
		echo "loading modules"
		module load ncbi-blast
		source activate python2std #biocluster specific option. loading biopython basically
	else 
		echo "no modules loaded"
	fi
	if [ -z "$8" ]
	then
	    echo "single sample option selected"
		echo "extract contigs"
		cut -f2 $1 | sort | uniq > contigs_to_extract.txt
		python $7 contigs_to_extract.txt $2
		$4 -db $3 -query extracted_contigs.fas -out $1_reciprocal.blast -outfmt 6 -num_threads $5
		rm contigs_to_extract.txt extracted_contigs.fas
	else
	  	echo "multiple sample option selected, number of samples: $(cat $8 | wc -l)"
	  	while read sample
	  	do
	  		echo "processing sample" $sample 
	  		echo "get contig names"
	  		cut -f2 $1/$sample".blast" | sort | uniq > contigs_to_extract.txt
	  		echo "get contigs"
	  		python $7 contigs_to_extract.txt $2/$sample
	  		echo "run blast"
	  		$4 -db $3 -query extracted_contigs.fas -out $1/$sample".blast_reciprocal.blast" -outfmt 6 -num_threads $5
	  		rm contigs_to_extract.txt extracted_contigs.fas
		done < $8
	fi
fi
echo "done"