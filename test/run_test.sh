echo "test"

echo "unzip the assembly"
gunzip sample.fas.gz

echo "make blast database"
makeblastdb -in sample.fas -dbtype nucl -parse_seqids

echo "run blastn"
blastn -query baits.fas -db sample.fas -outfmt 6 -out sample.fas.blast

echo "run alibaseq"
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
			python3 ../alibaseqPy3.py -b sample.fas.blast -t sample.fas -f S -x a -c 1 -o test_output -s test
		else
			echo "no correct python command found"
		fi
	fi
else
	pyV=$( python -V 2>&1 | cut -f2 -d" " | cut -f1 -d.)
	if [[ $pyV == 2 ]]
	then
		echo "python 2 found, executable python"
		python ../alibaseq.py -b sample.fas.blast -t sample.fas -f S -x a -c 1 -o test_output -s test
	elif [[ $pyV == 3 ]]
	then
		echo "python 3 found, executable python"
		python ../alibaseqPy3.py -b sample.fas.blast -t sample.fas -f S -x a -c 1 -o test_output -s test
	else
		echo "no correct python command found"
	fi
fi

echo "cleanup and gzip"
rm sample.fas.n*
gzip sample.fas