echo "test"

echo "unzip the assembly"
gunzip sample.fas.gz

echo "make blast database"
makeblastdb -in sample.fas -dbtype nucl -parse_seqids

echo "run blastn"
blastn -task dc-megablast -query baits.fas -db sample.fas -outfmt 6 -out sample.fas.blast

echo "run alibaseq"
alibaseqCommand1="-b sample.fas.blast -t sample.fas -f S -x a -c 1 -o best_output -s best --is -q ./baits/"
alibaseqCommand2="-b sample.fas.blast -t sample.fas -f S -x a -c 0 -o all_output -s all --is -q ./baits/"
if ! command -v python &> /dev/null
then
	if ! command -v python3 &> /dev/null
	then
		echo "python could not be found"
	else
		pyV=$( python3 -V 2>&1 | cut -f2 -d" " | cut -f1 -d.)
		if [[ $pyV == 3 ]]
		then
			echo "python 3 found, executable python3"
			python3 ../alibaseqPy3.py $alibaseqCommand1
			python3 ../alibaseqPy3.py $alibaseqCommand2
		else
			echo "no correct python command found"
		fi
	fi
else
	pyV=$( python -V 2>&1 | cut -f2 -d" " | cut -f1 -d.)
	if [[ $pyV == 2 ]]
	then
		echo "python 2 found, executable python"
		python ../alibaseq.py $alibaseqCommand1
		python ../alibaseq.py $alibaseqCommand2
	elif [[ $pyV == 3 ]]
	then
		echo "python 3 found, executable python"
		python ../alibaseqPy3.py $alibaseqCommand1
		python ../alibaseqPy3.py $alibaseqCommand2
	else
		echo "no correct python command found"
	fi
fi

echo "cleanup and gzip"
rm sample.fas.n*
gzip sample.fas