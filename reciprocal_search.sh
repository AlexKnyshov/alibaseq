echo "start"
if [[ $1 == "" ]]
then
	echo "########################################"
	echo "arguments"
	echo "\$1 folder with blast / hmm / lastz files"
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
	else 
		echo "no modules loaded"
	fi
	if [[ $4 == dc-megablast ]]
	then
		prog="blastn -task dc-megablast"
	else
		prog=$4
	fi

	pycmd=""
	if ! command -v python &> /dev/null
	then
		if ! command -v python3 &> /dev/null
		then
			echo "ERROR: python could not be found"
		else
			pyV=$( python3 -V 2>&1 | cut -f2 -d" " | cut -f1 -d.)
			if [[ $pyV == 3 ]]
			then
				echo "python 3 found, executable python3"
				pycmd="python3"
			else
				echo "ERROR: python3 executable has version 2"
			fi
		fi
	else
		pyV=$( python -V 2>&1 | cut -f2 -d" " | cut -f1 -d.)
		if [[ $pyV == 2 ]]
		then
			echo "python 2 found, executable python"
			pycmd="python"
		elif [[ $pyV == 3 ]]
		then
			echo "python 3 found, executable python"
			pycmd="python"
		else
			echo "ERROR: python version check"
		fi
	fi




	if [ -z "$8" ]
	then
	    echo "single sample option selected"
		echo "extract contigs"
		if [[ ${1##*.} == "hmmer" ]]
		then
			cut -f1 $1 | sort | uniq > contigs_to_extract.txt
		else
			cut -f2 $1 | sort | uniq > contigs_to_extract.txt
		fi
		# cut -f2 $1 | sort | uniq > contigs_to_extract.txt
		$pycmd $7 contigs_to_extract.txt $2
		$prog -db $3 -query extracted_contigs.fas -out $1"_reciprocal.blast" -outfmt 6 -num_threads $5 -evalue 0.0001
		rm contigs_to_extract.txt extracted_contigs.fas
	else
	  	echo "multiple sample option selected, number of samples: $(cat $8 | wc -l)"
	  	while read sample
	  	do
	  		echo "processing sample" $sample 
	  		echo "get contig names"
	  		if [[ -f $1/$sample".hmmer" ]]
	  		then
    			cut -f1 $1/$sample".hmmer" | sort | uniq > contigs_to_extract.txt
    			outname=$1/$sample".hmmer_reciprocal.blast"
    		elif [[ -f $1/$sample".lastz" ]]
    		then
    			cut -f2 $1/$sample".lastz" | sort | uniq > contigs_to_extract.txt
    			outname=$1/$sample".lastz_reciprocal.blast"
    		else
    			cut -f2 $1/$sample".blast" | sort | uniq > contigs_to_extract.txt
    			outname=$1/$sample".blast_reciprocal.blast"
    		fi
	  		# cut -f2 $1/$sample".blast" | sort | uniq > contigs_to_extract.txt
	  		echo "get contigs"
	  		$pycmd $7 contigs_to_extract.txt $2/$sample
	  		echo "run blast"
	  		$prog -db $3 -query extracted_contigs.fas -out $outname -outfmt 6 -num_threads $5 -evalue 0.0001
	  		rm contigs_to_extract.txt extracted_contigs.fas
		done < $8
	fi
fi
echo "done"